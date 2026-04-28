//! Reference-genome fingerprinting from BAM `@SQ` records.
//!
//! Identifies which canonical assembly a BAM was aligned to by matching
//! a small set of `(chromosome, length)` discriminators against the BAM
//! header. Output feeds `--gtf` cache resolution and the `fetch-references`
//! subcommand.

use anyhow::{Context, Result};
use rust_htslib::bam::{Read as BamRead, Reader};

/// Canonical reference assembly.
///
/// Add a new variant by extending [`KnownGenome::all`] and the discriminator
/// table in [`KnownGenome::discriminators`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KnownGenome {
    /// GRCh38 / hg38 (includes patches that don't change canonical lengths).
    Grch38,
    /// GRCh37 / hg19.
    Grch37,
    /// T2T-CHM13 v2.0.
    T2tChm13V2,
}

impl KnownGenome {
    /// All known genomes the fingerprinter recognises.
    pub fn all() -> &'static [KnownGenome] {
        &[
            KnownGenome::Grch38,
            KnownGenome::Grch37,
            KnownGenome::T2tChm13V2,
        ]
    }

    /// Stable short name used as the cache filename stem and for the
    /// envelope. Must be filesystem-safe.
    pub fn cache_name(&self) -> &'static str {
        match self {
            KnownGenome::Grch38 => "GRCh38",
            KnownGenome::Grch37 => "GRCh37",
            KnownGenome::T2tChm13V2 => "T2T-CHM13v2",
        }
    }

    /// Human-readable label.
    pub fn label(&self) -> &'static str {
        match self {
            KnownGenome::Grch38 => "GRCh38 (hg38)",
            KnownGenome::Grch37 => "GRCh37 (hg19)",
            KnownGenome::T2tChm13V2 => "T2T-CHM13 v2.0",
        }
    }

    /// Map a user-supplied genome name (case-insensitive) to a variant.
    /// Accepts canonical cache names plus common aliases.
    pub fn from_user_string(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "grch38" | "hg38" => Some(KnownGenome::Grch38),
            "grch37" | "hg19" => Some(KnownGenome::Grch37),
            "t2t-chm13v2" | "t2t" | "chm13" | "t2t-chm13" => Some(KnownGenome::T2tChm13V2),
            _ => None,
        }
    }

    fn discriminators(&self) -> &'static [(&'static str, u64)] {
        match self {
            KnownGenome::Grch38 => &[
                ("1", 248_956_422),
                ("2", 242_193_529),
                ("3", 198_295_559),
                ("4", 190_214_555),
                ("X", 156_040_895),
                ("Y", 57_227_415),
            ],
            KnownGenome::Grch37 => &[
                ("1", 249_250_621),
                ("2", 243_199_373),
                ("3", 198_022_430),
                ("4", 191_154_276),
                ("X", 155_270_560),
                ("Y", 59_373_566),
            ],
            KnownGenome::T2tChm13V2 => &[
                ("1", 248_387_328),
                ("2", 242_696_752),
                ("3", 201_105_948),
                ("4", 193_574_945),
                ("X", 154_259_566),
                ("Y", 62_460_029),
            ],
        }
    }
}

/// Floor at 3 — no two human references share three (chrom, length)
/// discriminators, and 3 tolerates BAMs that drop chrY or chr4.
const MIN_DISCRIMINATOR_HITS: usize = 3;

/// `chr`-prefix style of the BAM header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChromStyle {
    /// `chr1`, `chr2`, … (UCSC / GENCODE).
    Chr,
    /// `1`, `2`, … (Ensembl / NCBI).
    NoChr,
    /// Mixed or no canonical autosomes — fingerprint failed.
    Unknown,
}

impl std::fmt::Display for ChromStyle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChromStyle::Chr => f.write_str("chr-prefixed"),
            ChromStyle::NoChr => f.write_str("no-prefix"),
            ChromStyle::Unknown => f.write_str("unknown"),
        }
    }
}

/// Result of fingerprinting a BAM header.
#[derive(Debug, Clone)]
pub struct GenomeFingerprint {
    /// Identified assembly, or `None` if no known genome matched.
    pub genome: Option<KnownGenome>,
    /// Detected chromosome-naming style.
    pub chr_style: ChromStyle,
    /// Number of `@SQ` records inspected.
    pub n_records: usize,
    /// Discriminator hit count for the winning genome (0 when nothing
    /// matched).
    pub winning_hit_count: usize,
}

/// Read a BAM/SAM/CRAM header and return its `(chrom_name, length)` list.
/// Header-only — no record decoding.
pub fn read_bam_refs(path: &str) -> Result<Vec<(String, u64)>> {
    let reader = Reader::from_path(path)
        .with_context(|| format!("Cannot open alignment file '{}' for header read", path))?;
    let header = reader.header();
    Ok((0..header.target_count())
        .map(|tid| {
            let name = String::from_utf8_lossy(header.tid2name(tid)).into_owned();
            let len = header.target_len(tid).unwrap_or(0);
            (name, len)
        })
        .collect())
}

/// Open a BAM/SAM/CRAM and fingerprint its header.
pub fn fingerprint_bam(path: &str) -> Result<GenomeFingerprint> {
    Ok(fingerprint_records(&read_bam_refs(path)?))
}

/// Fingerprint a list of `(chrom_name, length)` tuples — pure function,
/// used by both [`fingerprint_bam`] and tests.
pub fn fingerprint_records(records: &[(String, u64)]) -> GenomeFingerprint {
    use std::collections::HashSet;
    let observed: HashSet<(String, u64)> = records
        .iter()
        .map(|(name, len)| (strip_chr_prefix(name).to_string(), *len))
        .collect();

    let (best_genome, best_hits) = KnownGenome::all()
        .iter()
        .map(|g| {
            let n = g
                .discriminators()
                .iter()
                .filter(|(c, l)| observed.contains(&((*c).to_string(), *l)))
                .count();
            (*g, n)
        })
        .max_by_key(|(_, n)| *n)
        .unwrap_or((KnownGenome::Grch38, 0));

    let genome = (best_hits >= MIN_DISCRIMINATOR_HITS).then_some(best_genome);

    GenomeFingerprint {
        genome,
        chr_style: detect_chr_style(records),
        n_records: records.len(),
        winning_hit_count: best_hits,
    }
}

/// Return the chromosome name with a leading `chr` (or `Chr`/`CHR`) stripped.
fn strip_chr_prefix(name: &str) -> &str {
    if name.len() >= 3 {
        let prefix = &name[..3];
        if prefix.eq_ignore_ascii_case("chr") {
            return &name[3..];
        }
    }
    name
}

fn detect_chr_style(records: &[(String, u64)]) -> ChromStyle {
    let mut has_chr = false;
    let mut has_nochr = false;
    for (name, _) in records {
        let core = strip_chr_prefix(name);
        if !is_canonical_chrom(core) {
            continue;
        }
        if name.starts_with("chr") || name.starts_with("Chr") || name.starts_with("CHR") {
            has_chr = true;
        } else {
            has_nochr = true;
        }
    }
    match (has_chr, has_nochr) {
        (true, false) => ChromStyle::Chr,
        (false, true) => ChromStyle::NoChr,
        _ => ChromStyle::Unknown,
    }
}

/// True if `name` is a canonical human chromosome name (no `chr` prefix):
/// `1..=22`, `X`, `Y`, `M`, or `MT`.
fn is_canonical_chrom(name: &str) -> bool {
    matches!(name, "X" | "Y" | "M" | "MT")
        || name
            .parse::<u32>()
            .map(|n| (1..=22).contains(&n))
            .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, len: u64) -> (String, u64) {
        (name.to_string(), len)
    }

    #[test]
    fn fingerprint_grch38_chr_prefix() {
        let recs = vec![
            rec("chr1", 248_956_422),
            rec("chr2", 242_193_529),
            rec("chr3", 198_295_559),
            rec("chrX", 156_040_895),
            rec("chrY", 57_227_415),
            rec("chrM", 16_569),
        ];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, Some(KnownGenome::Grch38));
        assert_eq!(fp.chr_style, ChromStyle::Chr);
        assert!(fp.winning_hit_count >= MIN_DISCRIMINATOR_HITS);
    }

    #[test]
    fn fingerprint_grch37_no_prefix() {
        let recs = vec![
            rec("1", 249_250_621),
            rec("2", 243_199_373),
            rec("3", 198_022_430),
            rec("X", 155_270_560),
        ];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, Some(KnownGenome::Grch37));
        assert_eq!(fp.chr_style, ChromStyle::NoChr);
    }

    #[test]
    fn fingerprint_t2t_chm13() {
        let recs = vec![
            rec("chr1", 248_387_328),
            rec("chr2", 242_696_752),
            rec("chr3", 201_105_948),
        ];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, Some(KnownGenome::T2tChm13V2));
    }

    #[test]
    fn fingerprint_synthetic_returns_none() {
        // The repo's test.bam — synthetic chr1/chr2 of 20000 bp.
        let recs = vec![rec("chr1", 20_000), rec("chr2", 20_000)];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, None);
        assert_eq!(fp.chr_style, ChromStyle::Chr);
    }

    #[test]
    fn fingerprint_partial_match_below_floor_returns_none() {
        // Only chr1 length matches GRCh38; need at least 3 hits to declare.
        let recs = vec![rec("chr1", 248_956_422), rec("chr2", 1), rec("chr3", 1)];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, None);
    }

    #[test]
    fn fingerprint_with_decoys_still_matches() {
        let recs = vec![
            rec("chr1", 248_956_422),
            rec("chr2", 242_193_529),
            rec("chr3", 198_295_559),
            rec("chrEBV", 171_823),
            rec("chr1_KI270706v1_random", 175_055),
        ];
        let fp = fingerprint_records(&recs);
        assert_eq!(fp.genome, Some(KnownGenome::Grch38));
    }

    #[test]
    fn strip_chr_prefix_handles_case_variants() {
        assert_eq!(strip_chr_prefix("chr1"), "1");
        assert_eq!(strip_chr_prefix("Chr1"), "1");
        assert_eq!(strip_chr_prefix("CHR1"), "1");
        assert_eq!(strip_chr_prefix("1"), "1");
        assert_eq!(strip_chr_prefix("hr"), "hr");
        assert_eq!(strip_chr_prefix("MT"), "MT");
    }

    #[test]
    fn cache_names_are_unique_and_path_safe() {
        use std::collections::HashSet;
        let names: HashSet<&str> = KnownGenome::all().iter().map(|g| g.cache_name()).collect();
        assert_eq!(names.len(), KnownGenome::all().len());
        for g in KnownGenome::all() {
            let n = g.cache_name();
            assert!(!n.contains('/'));
            assert!(!n.contains(' '));
            assert!(!n.is_empty());
        }
    }

    #[test]
    fn from_user_string_resolves_canonical_and_aliases() {
        assert_eq!(
            KnownGenome::from_user_string("GRCh38"),
            Some(KnownGenome::Grch38)
        );
        assert_eq!(
            KnownGenome::from_user_string("hg38"),
            Some(KnownGenome::Grch38)
        );
        assert_eq!(
            KnownGenome::from_user_string("HG19"),
            Some(KnownGenome::Grch37)
        );
        assert_eq!(
            KnownGenome::from_user_string("  t2t  "),
            Some(KnownGenome::T2tChm13V2)
        );
        assert_eq!(KnownGenome::from_user_string("nope"), None);
    }
}
