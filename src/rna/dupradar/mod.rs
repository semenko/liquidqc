//! dupRadar analysis modules.
//!
//! Duplication rate analysis for RNA-Seq data, including BAM read counting,
//! duplication matrix construction, logistic regression fitting, and plot generation.

pub mod counting;
pub mod dupmatrix;
pub mod fitting;
pub mod plots;
