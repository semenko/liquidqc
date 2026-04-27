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
