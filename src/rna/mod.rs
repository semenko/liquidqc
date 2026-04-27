//! RNA-Seq quality control and analysis modules.
//!
//! Contains dupRadar duplication rate analysis, featureCounts-compatible output,
//! RSeQC tool reimplementations, and the liquidqc Phase 2 Tier 1 fragmentomics
//! accumulators (end-motifs, soft-clip k-mers, fragment-size bins, periodicity FFT).

pub mod bam_flags;
pub mod cpp_rng;
pub mod dupradar;
pub mod featurecounts;
pub mod fragmentomics;
pub mod per_gene;
pub mod preseq;
pub mod qualimap;
pub mod rseqc;
