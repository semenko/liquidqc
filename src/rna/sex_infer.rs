//! Sex inference from chrY:autosome ratio + XIST + RPS4Y1.
//!
//! Pure finalize-time computation against [`crate::rna::chrom_metrics`] and
//! the per-gene counts. Always emits `predicted_sex`; the
//! `sex_swap_warning` qc_flag was intentionally dropped (no required user
//! metadata in v1) — downstream consumers compare against their own metadata.

use indexmap::IndexMap;
use serde::Serialize;

use crate::gtf::Gene;
use crate::rna::chrom_metrics::ChromMetricsResult;
use crate::rna::dupradar::counting::GeneCounts;

/// Y-chromosome / autosome read fraction above which we consider the chrY
/// signal "strong". Tuned for the typical male whole-blood cfRNA range
/// (~0.005) while sitting well above the female background (~0.0001 from
/// pseudoautosomal-region reads). Conservative: lower thresholds inflate
/// false-male calls in low-Y females.
const Y_AUTOSOME_RATIO_STRONG: f64 = 1e-3;

/// XIST read fraction (`xist_reads / fc_assigned`) above which the XIST
/// signal is "strong". Females typically express XIST at ~10⁻⁴ of total
/// assigned reads in cfRNA; males have orders-of-magnitude less. Tuned to
/// the female regime (10⁻⁴) — the earlier 10⁻⁵ value sat in male leakage
/// territory and forced "ambiguous" calls on male blood samples with even
/// a few residual XIST reads.
const XIST_FRACTION_STRONG: f64 = 1e-4;

/// RPS4Y1 read fraction above which the male signal from this gene alone
/// is "strong". RPS4Y1 is a Y-chromosome ribosomal protein gene that is
/// expressed in males and absent in females; its read count is an
/// orthogonal male signal independent of total chrY mapping density (which
/// can pick up reads on PAR / pseudoautosomal regions in females).
const RPS4Y1_FRACTION_STRONG: f64 = 1e-5;

/// Minimum total mapped reads required to attempt a confident call. Below
/// this we report `unknown` regardless of ratios because counting noise
/// dominates at low depth.
const MIN_MAPPED_FOR_CALL: u64 = 1_000;

/// Schema-stable predicted sex value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum PredictedSex {
    Male,
    Female,
    Ambiguous,
    Unknown,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct SexInferenceResult {
    pub chr_y_reads: u64,
    pub autosome_reads: u64,
    /// `chr_y_reads / autosome_reads`. 0.0 when autosome_reads is 0.
    pub chr_y_autosome_read_ratio: f64,
    pub xist_reads: u64,
    pub rps4y1_reads: u64,
    /// `xist_reads / denominator`. Denominator is `featurecounts.assigned`
    /// (or 0).
    pub xist_fraction: f64,
    /// `rps4y1_reads / denominator`.
    pub rps4y1_fraction: f64,
    pub denominator: u64,
    pub predicted_sex: PredictedSex,
}

/// Compute the sex-inference result. Returns `predicted_sex = Unknown` and
/// zeroed metrics when the inputs are insufficient (no chrom-metrics, no
/// genes resolved by symbol, or `mapped < MIN_MAPPED_FOR_CALL`).
pub fn compute(
    chrom_metrics: Option<&ChromMetricsResult>,
    genes: &IndexMap<String, Gene>,
    gene_counts: &IndexMap<String, GeneCounts>,
    fc_assigned: u64,
) -> SexInferenceResult {
    let (chr_y_reads, autosome_reads, ratio) = match chrom_metrics {
        Some(cm) => derive_chr_y_autosome(cm),
        None => (0, 0, 0.0),
    };

    let symbol_reads = read_counts_for_symbols(genes, gene_counts, &["XIST", "RPS4Y1"]);
    let xist_reads = symbol_reads[0];
    let rps4y1_reads = symbol_reads[1];

    let xist_fraction = crate::rna::safe_fraction(xist_reads, fc_assigned);
    let rps4y1_fraction = crate::rna::safe_fraction(rps4y1_reads, fc_assigned);

    let predicted_sex = predict(autosome_reads, ratio, xist_fraction, rps4y1_fraction);

    SexInferenceResult {
        chr_y_reads,
        autosome_reads,
        chr_y_autosome_read_ratio: ratio,
        xist_reads,
        rps4y1_reads,
        xist_fraction,
        rps4y1_fraction,
        denominator: fc_assigned,
        predicted_sex,
    }
}

fn derive_chr_y_autosome(cm: &ChromMetricsResult) -> (u64, u64, f64) {
    let mut chr_y: u64 = 0;
    let mut autosome: u64 = 0;
    for (name, entry) in &cm.contigs {
        if is_chr_y(name) {
            chr_y = chr_y.saturating_add(entry.mapped);
        } else if is_autosome(name) {
            autosome = autosome.saturating_add(entry.mapped);
        }
    }
    let ratio = if autosome == 0 {
        0.0
    } else {
        chr_y as f64 / autosome as f64
    };
    (chr_y, autosome, ratio)
}

fn is_chr_y(name: &str) -> bool {
    let trimmed = name
        .strip_prefix("chr")
        .or_else(|| name.strip_prefix("CHR"));
    let core = trimmed.unwrap_or(name);
    core.eq_ignore_ascii_case("Y")
}

fn is_autosome(name: &str) -> bool {
    let trimmed = name
        .strip_prefix("chr")
        .or_else(|| name.strip_prefix("CHR"));
    let core = trimmed.unwrap_or(name);
    if core.is_empty() {
        return false;
    }
    // Autosomes are 1..=22 (human). Match exactly so we don't grab
    // alt/scaffold/random contigs that contain digits.
    core.chars().all(|c| c.is_ascii_digit()) && {
        let n: u32 = core.parse().unwrap_or(0);
        (1..=22).contains(&n)
    }
}

/// Single-pass resolver for several gene symbols at once. Walks `genes` once,
/// summing `fc_reads` into a fixed-shape array of totals (one per requested
/// symbol). Avoids the O(N · S) full-rescan of an earlier per-symbol loop.
fn read_counts_for_symbols<const N: usize>(
    genes: &IndexMap<String, Gene>,
    gene_counts: &IndexMap<String, GeneCounts>,
    symbols: &[&str; N],
) -> [u64; N] {
    let mut totals = [0u64; N];
    for (gid, gene) in genes {
        let Some(g_sym) = gene.attributes.get("gene_name") else {
            continue;
        };
        let Some(idx) = symbols.iter().position(|s| g_sym == *s) else {
            continue;
        };
        if let Some(c) = gene_counts.get(gid) {
            totals[idx] = totals[idx].saturating_add(c.fc_reads);
        }
    }
    totals
}

fn predict(
    autosome_reads: u64,
    y_ratio: f64,
    xist_fraction: f64,
    rps4y1_fraction: f64,
) -> PredictedSex {
    if autosome_reads < MIN_MAPPED_FOR_CALL {
        return PredictedSex::Unknown;
    }
    // Two orthogonal male signals: chrY/autosome mapping density (sensitive
    // but pollutable by PAR mismaps) and RPS4Y1 expression (specific to
    // male chrY transcription). Either alone suffices for a "male" call.
    let male_strong = y_ratio > Y_AUTOSOME_RATIO_STRONG || rps4y1_fraction > RPS4Y1_FRACTION_STRONG;
    let female_strong = xist_fraction > XIST_FRACTION_STRONG;
    match (male_strong, female_strong) {
        (true, false) => PredictedSex::Male,
        (false, true) => PredictedSex::Female,
        (true, true) => PredictedSex::Ambiguous,
        (false, false) => PredictedSex::Unknown,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::chrom_metrics::{ChromMetricsEntry, ChromMetricsResult};
    use std::collections::{BTreeMap, HashMap};

    fn cm(entries: &[(&str, u64)]) -> ChromMetricsResult {
        let mut contigs = BTreeMap::new();
        let mut total_mapped: u64 = 0;
        for (n, mapped) in entries {
            total_mapped += mapped;
            contigs.insert(
                n.to_string(),
                ChromMetricsEntry {
                    contig_length: 1_000_000,
                    total_records: *mapped,
                    mapped: *mapped,
                    duplicates: 0,
                    read_fraction: 0.0,
                    proper_pair_count: 0,
                    proper_pair_insert_size_mean: 0.0,
                },
            );
        }
        ChromMetricsResult {
            total_mapped,
            contigs,
        }
    }

    fn gene(symbol: &str) -> Gene {
        let mut attrs = HashMap::new();
        attrs.insert("gene_name".to_string(), symbol.to_string());
        Gene {
            gene_id: format!("ensg_{symbol}"),
            chrom: "chr1".to_string(),
            start: 0,
            end: 100,
            strand: '+',
            exons: vec![],
            effective_length: 0,
            attributes: attrs,
            transcripts: vec![],
        }
    }

    fn gene_counts_for(reads: &[(&str, u64)]) -> IndexMap<String, GeneCounts> {
        let mut m: IndexMap<String, GeneCounts> = IndexMap::new();
        for (gid, n) in reads {
            m.insert(
                (*gid).to_string(),
                GeneCounts {
                    fc_reads: *n,
                    ..Default::default()
                },
            );
        }
        m
    }

    #[test]
    fn chr_y_match_handles_chr_prefix() {
        assert!(is_chr_y("chrY"));
        assert!(is_chr_y("CHRY"));
        assert!(is_chr_y("Y"));
        assert!(!is_chr_y("chrY_KI270740v1_random"));
        assert!(!is_chr_y("chrYY")); // double-Y is not Y
    }

    #[test]
    fn autosome_match_only_human_1_to_22() {
        for i in 1..=22 {
            assert!(is_autosome(&format!("chr{}", i)));
            assert!(is_autosome(&format!("{}", i)));
        }
        assert!(!is_autosome("chrX"));
        assert!(!is_autosome("chrY"));
        assert!(!is_autosome("chrM"));
        assert!(!is_autosome("chr23"));
        assert!(!is_autosome("chr1_random"));
    }

    #[test]
    fn predicts_male_with_strong_y_and_weak_xist() {
        let metrics = cm(&[("chr1", 10_000), ("chrY", 50)]);
        let mut genes = IndexMap::new();
        genes.insert("ensg_RPS4Y1".to_string(), gene("RPS4Y1"));
        let counts = gene_counts_for(&[("ensg_RPS4Y1", 100)]);
        let r = compute(Some(&metrics), &genes, &counts, 10_050);
        assert!(r.chr_y_autosome_read_ratio > Y_AUTOSOME_RATIO_STRONG);
        assert_eq!(r.predicted_sex, PredictedSex::Male);
    }

    #[test]
    fn predicts_female_with_weak_y_and_strong_xist() {
        let metrics = cm(&[("chr1", 10_000), ("chrY", 1)]);
        let mut genes = IndexMap::new();
        genes.insert("ensg_XIST".to_string(), gene("XIST"));
        let counts = gene_counts_for(&[("ensg_XIST", 50)]);
        let r = compute(Some(&metrics), &genes, &counts, 10_001);
        assert!(r.chr_y_autosome_read_ratio < Y_AUTOSOME_RATIO_STRONG);
        assert!(r.xist_fraction > XIST_FRACTION_STRONG);
        assert_eq!(r.predicted_sex, PredictedSex::Female);
    }

    #[test]
    fn predicts_ambiguous_with_both_signals_strong() {
        let metrics = cm(&[("chr1", 10_000), ("chrY", 100)]);
        let mut genes = IndexMap::new();
        genes.insert("ensg_XIST".to_string(), gene("XIST"));
        let counts = gene_counts_for(&[("ensg_XIST", 50)]);
        let r = compute(Some(&metrics), &genes, &counts, 10_100);
        assert_eq!(r.predicted_sex, PredictedSex::Ambiguous);
    }

    #[test]
    fn predicts_unknown_at_low_depth() {
        // Strong Y ratio numerically, but only 5 autosome reads → unknown.
        let metrics = cm(&[("chr1", 5), ("chrY", 5)]);
        let genes: IndexMap<String, Gene> = IndexMap::new();
        let counts: IndexMap<String, GeneCounts> = IndexMap::new();
        let r = compute(Some(&metrics), &genes, &counts, 10);
        assert_eq!(r.predicted_sex, PredictedSex::Unknown);
    }

    #[test]
    fn unknown_when_chrom_metrics_missing() {
        let genes: IndexMap<String, Gene> = IndexMap::new();
        let counts: IndexMap<String, GeneCounts> = IndexMap::new();
        let r = compute(None, &genes, &counts, 0);
        assert_eq!(r.predicted_sex, PredictedSex::Unknown);
    }
}
