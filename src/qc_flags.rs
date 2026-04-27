//! QC-flag rule engine for the v1 envelope.
//!
//! Pure module: takes references, returns `Vec<String>`. Wired into
//! `process_single_bam()` after every per-tool result is finalized but before
//! envelope construction. Thresholds are intentionally conservative; they are
//! also the user-visible interpretive layer of the envelope, so the rules are
//! kept simple and documented inline.
//!
//! Two flags defined by the schema are not emitted in Phase 1:
//! - [`QcFlag::CfdnaContaminationSuspected`] — needs fragmentomics signal
//!   (Phase 2+).
//! - [`QcFlag::SexSwapWarning`] — needs user-supplied sex metadata.
//!
//! Their enum variants are kept so the wire format is stable.

use crate::cli::Strandedness;
use crate::rna::dupradar::counting::CountResult;
use crate::rna::dupradar::fitting::FitResult;
use crate::rna::preseq::PreseqResult;
use crate::rna::rseqc::bam_stat::BamStatResult;

/// Schema-defined QC-flag identifiers.
///
/// The string forms (returned by [`QcFlag::as_str`]) match the `qc_flags`
/// enum in `schema/v1/liquidqc.schema.json`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QcFlag {
    SingleEndNoTlen,
    LowPairedFraction,
    LowMappingRate,
    /// Phase 2 — requires fragmentomics signal; not emitted in Phase 1.
    #[allow(dead_code)]
    CfdnaContaminationSuspected,
    HighRrnaFraction,
    TagmentationEndmotifBias,
    UnstrandedEndbiasUnreliable,
    ReadLengthCapsLongFragments,
    /// Phase 4+ — requires user-supplied sex metadata; not emitted in Phase 1.
    #[allow(dead_code)]
    SexSwapWarning,
    LowComplexityLibrary,
    HighDuplicationNonbiological,
}

impl QcFlag {
    pub fn as_str(self) -> &'static str {
        match self {
            QcFlag::SingleEndNoTlen => "single_end_no_tlen",
            QcFlag::LowPairedFraction => "low_paired_fraction",
            QcFlag::LowMappingRate => "low_mapping_rate",
            QcFlag::CfdnaContaminationSuspected => "cfdna_contamination_suspected",
            QcFlag::HighRrnaFraction => "high_rrna_fraction",
            QcFlag::TagmentationEndmotifBias => "tagmentation_endmotif_bias",
            QcFlag::UnstrandedEndbiasUnreliable => "unstranded_endbias_unreliable",
            QcFlag::ReadLengthCapsLongFragments => "read_length_caps_long_fragments",
            QcFlag::SexSwapWarning => "sex_swap_warning",
            QcFlag::LowComplexityLibrary => "low_complexity_library",
            QcFlag::HighDuplicationNonbiological => "high_duplication_nonbiological",
        }
    }
}

// Phase 1 thresholds. Tunable; see CHANGELOG when changing.
const LOW_PAIRED_FRACTION_THRESHOLD: f64 = 0.7;
const LOW_MAPPING_RATE_THRESHOLD: f64 = 0.7;
const HIGH_RRNA_FRACTION_THRESHOLD: f64 = 0.4;
const READ_LENGTH_CAPS_LONG_FRAGMENTS_MAX: u64 = 75;
const HIGH_DUPLICATION_INTERCEPT_THRESHOLD: f64 = 0.5;
const HIGH_DUPLICATION_MIN_GENES: u64 = 1000;
/// Slope of preseq's last extrapolated segment below this is "flat".
const LOW_COMPLEXITY_PRESEQ_SLOPE_THRESHOLD: f64 = 0.2;
const LOW_COMPLEXITY_FALLBACK_DUP_FRACTION: f64 = 0.95;
const LOW_COMPLEXITY_FALLBACK_MAX_BIOTYPE_ASSIGNED: u64 = 1_000_000;

/// All accumulator references needed to evaluate Phase 1 rules.
pub struct QcContext<'a> {
    pub paired_end: bool,
    pub strandedness: Strandedness,
    pub library_prep: &'a str,
    pub bam_stat: Option<&'a BamStatResult>,
    pub count_result: Option<&'a CountResult>,
    pub preseq: Option<&'a PreseqResult>,
    pub dupradar_fit: Option<&'a FitResult>,
    pub dupradar_genes_with_reads: u64,
    pub featurecounts_biotype_rrna: u64,
}

/// Evaluate Phase 1 rules and return matching flag names in schema order.
pub fn evaluate(ctx: &QcContext<'_>) -> Vec<String> {
    let mut flags: Vec<QcFlag> = Vec::new();

    if !ctx.paired_end {
        flags.push(QcFlag::SingleEndNoTlen);
    }

    if let Some(stat) = ctx.bam_stat {
        if ctx.paired_end && stat.mapped > 0 {
            let frac = stat.proper_pairs as f64 / stat.mapped as f64;
            if frac < LOW_PAIRED_FRACTION_THRESHOLD {
                flags.push(QcFlag::LowPairedFraction);
            }
        }
        if stat.total_records > 0 {
            let frac = stat.mapped as f64 / stat.total_records as f64;
            if frac < LOW_MAPPING_RATE_THRESHOLD {
                flags.push(QcFlag::LowMappingRate);
            }
        }
        if stat.max_len > 0 && stat.max_len <= READ_LENGTH_CAPS_LONG_FRAGMENTS_MAX {
            flags.push(QcFlag::ReadLengthCapsLongFragments);
        }
    }

    if let Some(c) = ctx.count_result {
        if c.fc_biotype_assigned > 0 {
            let frac = ctx.featurecounts_biotype_rrna as f64 / c.fc_biotype_assigned as f64;
            if frac > HIGH_RRNA_FRACTION_THRESHOLD {
                flags.push(QcFlag::HighRrnaFraction);
            }
        }
    }

    if ctx.library_prep.to_ascii_lowercase().contains("tagment") {
        flags.push(QcFlag::TagmentationEndmotifBias);
    }

    if matches!(ctx.strandedness, Strandedness::Unstranded) {
        flags.push(QcFlag::UnstrandedEndbiasUnreliable);
    }

    if low_complexity_library(ctx) {
        flags.push(QcFlag::LowComplexityLibrary);
    }

    if let Some(fit) = ctx.dupradar_fit {
        if fit.intercept > HIGH_DUPLICATION_INTERCEPT_THRESHOLD
            && ctx.dupradar_genes_with_reads >= HIGH_DUPLICATION_MIN_GENES
        {
            flags.push(QcFlag::HighDuplicationNonbiological);
        }
    }

    flags.into_iter().map(|f| f.as_str().to_string()).collect()
}

/// Preseq-based low-complexity detection with a duplication-fraction fallback
/// for runs where preseq was disabled.
fn low_complexity_library(ctx: &QcContext<'_>) -> bool {
    if let Some(p) = ctx.preseq {
        if let Some(slope) = preseq_tail_slope(p) {
            return slope < LOW_COMPLEXITY_PRESEQ_SLOPE_THRESHOLD;
        }
    }
    if let (Some(stat), Some(c)) = (ctx.bam_stat, ctx.count_result) {
        if stat.mapped > 0 {
            let dup_frac = stat.duplicates as f64 / stat.mapped as f64;
            if dup_frac > LOW_COMPLEXITY_FALLBACK_DUP_FRACTION
                && c.fc_biotype_assigned < LOW_COMPLEXITY_FALLBACK_MAX_BIOTYPE_ASSIGNED
            {
                return true;
            }
        }
    }
    false
}

/// Slope of preseq's last extrapolation segment: ΔE[distinct] / Δreads.
///
/// A flat tail (slope near 0) means the library has saturated and additional
/// sequencing won't recover new molecules.
fn preseq_tail_slope(p: &PreseqResult) -> Option<f64> {
    if p.curve.len() < 2 {
        return None;
    }
    let n = p.curve.len();
    let (x1, y1, _, _) = p.curve[n - 2];
    let (x2, y2, _, _) = p.curve[n - 1];
    let dx = x2 - x1;
    if dx <= 0.0 {
        return None;
    }
    Some((y2 - y1) / dx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::dupradar::counting::{CountResult, GeneCounts};
    use indexmap::IndexMap;

    fn empty_count_result() -> CountResult {
        CountResult {
            gene_counts: IndexMap::<String, GeneCounts>::new(),
            stat_total_reads: 0,
            stat_assigned: 0,
            stat_ambiguous: 0,
            stat_no_features: 0,
            stat_total_fragments: 0,
            stat_total_mapped: 0,
            stat_total_dup: 0,
            stat_singleton_unmapped_mates: 0,
            fc_assigned: 0,
            fc_ambiguous: 0,
            fc_no_features: 0,
            fc_multimapping: 0,
            fc_unmapped: 0,
            fc_singleton: 0,
            fc_chimera: 0,
            biotype_reads: vec![],
            biotype_names: vec![],
            fc_biotype_assigned: 0,
            fc_biotype_ambiguous: 0,
            fc_biotype_no_features: 0,
            rseqc: None,
            qualimap: None,
            fragmentomics: None,
        }
    }

    fn base_ctx<'a>(library_prep: &'a str) -> QcContext<'a> {
        QcContext {
            paired_end: true,
            strandedness: Strandedness::Forward,
            library_prep,
            bam_stat: None,
            count_result: None,
            preseq: None,
            dupradar_fit: None,
            dupradar_genes_with_reads: 0,
            featurecounts_biotype_rrna: 0,
        }
    }

    #[test]
    fn paired_end_with_clean_inputs_emits_no_flags() {
        let flags = evaluate(&base_ctx("neb_next_ultra_ii_directional"));
        assert!(flags.is_empty(), "unexpected flags: {flags:?}");
    }

    #[test]
    fn single_end_emits_single_end_no_tlen() {
        let mut ctx = base_ctx("neb_next");
        ctx.paired_end = false;
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"single_end_no_tlen".to_string()));
    }

    #[test]
    fn low_paired_fraction_fires_below_threshold() {
        let mut stat = BamStatResult::default();
        stat.mapped = 1000;
        stat.proper_pairs = 100; // 10% — well below threshold
        stat.total_records = 1100;
        let mut ctx = base_ctx("neb_next");
        ctx.bam_stat = Some(&stat);
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"low_paired_fraction".to_string()));
    }

    #[test]
    fn low_mapping_rate_fires_below_threshold() {
        let mut stat = BamStatResult::default();
        stat.total_records = 1000;
        stat.mapped = 100; // 10% mapping rate
        stat.proper_pairs = 100;
        let mut ctx = base_ctx("neb_next");
        ctx.bam_stat = Some(&stat);
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"low_mapping_rate".to_string()));
    }

    #[test]
    fn high_rrna_fraction_fires() {
        let mut count = empty_count_result();
        count.fc_biotype_assigned = 1000;
        let mut ctx = base_ctx("neb_next");
        ctx.count_result = Some(&count);
        ctx.featurecounts_biotype_rrna = 600; // 60% rRNA
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"high_rrna_fraction".to_string()));
    }

    #[test]
    fn tagmentation_substring_fires() {
        let ctx = base_ctx("illumina_rna_prep_with_enrichment_tagmentation");
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"tagmentation_endmotif_bias".to_string()));
    }

    #[test]
    fn tagmentation_substring_case_insensitive() {
        let ctx = base_ctx("Some_TAGMENT_thing");
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"tagmentation_endmotif_bias".to_string()));
    }

    #[test]
    fn unstranded_fires() {
        let mut ctx = base_ctx("neb_next");
        ctx.strandedness = Strandedness::Unstranded;
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"unstranded_endbias_unreliable".to_string()));
    }

    #[test]
    fn short_reads_fire_read_length_caps() {
        let mut stat = BamStatResult::default();
        stat.max_len = 75; // boundary value: still fires
        let mut ctx = base_ctx("neb_next");
        ctx.bam_stat = Some(&stat);
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"read_length_caps_long_fragments".to_string()));
    }

    #[test]
    fn long_reads_do_not_fire_read_length_caps() {
        let mut stat = BamStatResult::default();
        stat.max_len = 150;
        let mut ctx = base_ctx("neb_next");
        ctx.bam_stat = Some(&stat);
        let flags = evaluate(&ctx);
        assert!(!flags.contains(&"read_length_caps_long_fragments".to_string()));
    }

    #[test]
    fn high_dup_intercept_fires() {
        let fit = FitResult {
            beta0: 0.0,
            beta1: 0.0,
            intercept: 0.7,
            slope: 0.5,
        };
        let mut ctx = base_ctx("neb_next");
        ctx.dupradar_fit = Some(&fit);
        ctx.dupradar_genes_with_reads = 5000;
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"high_duplication_nonbiological".to_string()));
    }

    #[test]
    fn high_dup_intercept_with_few_genes_does_not_fire() {
        let fit = FitResult {
            beta0: 0.0,
            beta1: 0.0,
            intercept: 0.7,
            slope: 0.5,
        };
        let mut ctx = base_ctx("neb_next");
        ctx.dupradar_fit = Some(&fit);
        ctx.dupradar_genes_with_reads = 50; // below MIN_GENES
        let flags = evaluate(&ctx);
        assert!(!flags.contains(&"high_duplication_nonbiological".to_string()));
    }

    #[test]
    fn low_complexity_preseq_flat_curve_fires() {
        let preseq = PreseqResult {
            curve: vec![(1.0, 1.0, 0.9, 1.1), (1_000_000.0, 1.0001, 0.99, 1.001)],
        };
        let mut ctx = base_ctx("neb_next");
        ctx.preseq = Some(&preseq);
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"low_complexity_library".to_string()));
    }

    #[test]
    fn low_complexity_fallback_high_dup_low_assigned() {
        let mut stat = BamStatResult::default();
        stat.mapped = 1000;
        stat.duplicates = 970; // 97% dup
        let mut count = empty_count_result();
        count.fc_biotype_assigned = 50_000; // < threshold
        let mut ctx = base_ctx("neb_next");
        ctx.bam_stat = Some(&stat);
        ctx.count_result = Some(&count);
        let flags = evaluate(&ctx);
        assert!(flags.contains(&"low_complexity_library".to_string()));
    }
}
