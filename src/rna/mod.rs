//! RNA-Seq quality control and analysis modules.
//!
//! Contains dupRadar duplication rate analysis, featureCounts-compatible output,
//! RSeQC tool reimplementations, the Phase 2 Tier 1 fragmentomics accumulators
//! (end-motifs, soft-clip k-mers, fragment-size bins, periodicity FFT), the
//! Phase 3 Tier 2 per-gene sparse table, and the Phase 4 Tier 3 add-ons
//! (per-chromosome metrics, cycle quality, splice-site dinucleotides, gene
//! classes, marker panels, sex inference, saturation curve, SNP fingerprint).

pub mod bam_flags;
pub mod chrom_metrics;
pub mod cpp_rng;
pub mod cycle_quality;
pub mod dupradar;
pub mod featurecounts;
pub mod fragmentomics;
pub mod gene_class;
pub mod panels;
pub mod per_gene;
pub mod preseq;
pub mod qualimap;
pub mod rseqc;
pub mod saturation;
pub mod sex_infer;
pub mod snp_fingerprint;
pub mod splice_dinuc;

/// `numerator / denominator` as f64, with `denominator == 0` mapped to `0.0`.
///
/// The "fraction with zero-safe denominator" pattern recurs at every
/// finalize site (gene-class, sex-inference, fragment-size, periodicity bands,
/// per-gene rates, panel aggregation, snp depth, saturation, etc.). Pulling
/// it into one place keeps the divide-by-zero contract uniform across the
/// envelope.
pub fn safe_fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn safe_fraction_zero_denominator_is_zero() {
        assert_eq!(safe_fraction(7, 0), 0.0);
        assert_eq!(safe_fraction(0, 0), 0.0);
    }

    #[test]
    fn safe_fraction_basic_division() {
        assert!((safe_fraction(1, 2) - 0.5).abs() < 1e-12);
        assert_eq!(safe_fraction(0, 100), 0.0);
        assert_eq!(safe_fraction(100, 100), 1.0);
    }
}
