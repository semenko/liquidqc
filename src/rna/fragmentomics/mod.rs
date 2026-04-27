//! Phase 2 Tier 1 fragmentomics accumulators.
//!
//! Hooks into the existing per-chromosome BAM dispatcher in
//! [`crate::rna::dupradar::counting::count_chromosome`]. Each worker owns a
//! [`FragmentomicsAccumulators`] container alongside its
//! [`crate::rna::rseqc::accumulators::RseqcAccumulators`]; per-record
//! dispatch is a single extra method call per record. Workers are merged at
//! the same point as the RSeQC accumulators.
//!
//! See module-level doc comments on each submodule for the per-feature
//! algorithm and conventions:
//! - [`fragment_size`] — `|TLEN|` histogram + 8-bin output + summary fractions.
//! - [`end_motifs`] — 5'/3' 4-mers from the reference FASTA. Requires a
//!   `--fasta` path; otherwise this submodule's accumulator is `None`.
//! - [`soft_clips`] — read-anchored soft-clip 4-mers + per-end clip rates.
//! - [`periodicity`] — pure FFT over the merged TLEN histogram (no per-read state).
//!
//! All accumulators in the per-record hot path are infallible. The only
//! fallible step is opening the FASTA reader at construction time.

pub mod common;
pub mod end_motifs;
pub mod fragment_size;
pub mod periodicity;
pub mod soft_clips;

use anyhow::Result;
use rust_htslib::bam::record::Record;
use std::path::Path;

pub use end_motifs::{EndMotifAccum, EndMotifLookup, EndMotifResult};
pub use fragment_size::{FragmentSizeAccum, FragmentSizeBins, FragmentSizeResult};
pub use periodicity::{compute as compute_periodicity, PeriodicityResult};
pub use soft_clips::{SoftClipAccum, SoftClipKmerByLen, SoftClipResult};

/// Run-time inputs for the four Phase 2 Tier 1 accumulators.
///
/// Each accumulator is automatically inhibited when its prerequisite is
/// missing: end-motifs require `--fasta`; end-motifs and fragment-size
/// require paired-end. Soft-clips run on any input.
#[derive(Debug, Clone, Copy)]
pub struct FragmentomicsConfig {
    pub paired_end: bool,
    pub mapq_cut: u8,
}

/// Per-worker container. Each chromosome worker constructs one of these
/// alongside its [`crate::rna::rseqc::accumulators::RseqcAccumulators`] and
/// dispatches into it from the same record-iteration loop.
#[derive(Debug)]
pub struct FragmentomicsAccumulators {
    pub end_motifs: Option<EndMotifAccum>,
    pub soft_clips: Option<SoftClipAccum>,
    pub fragment_size: Option<FragmentSizeAccum>,
    cfg: FragmentomicsConfig,
}

impl FragmentomicsAccumulators {
    /// Construct a per-worker container.
    ///
    /// `fasta_path` must be `Some(path)` when the user passed `--fasta`,
    /// else `None`. End-motif extraction is silently disabled when
    /// `fasta_path` is `None` (the dispatcher emits the
    /// `end_motifs_skipped_no_fasta` qc_flag at envelope-build time).
    ///
    /// Single-end input (`!cfg.paired_end`) disables the two paired-only
    /// accumulators (end-motifs, fragment-size). Soft-clips still run.
    pub fn new(cfg: &FragmentomicsConfig, fasta_path: Option<&Path>) -> Result<Self> {
        let end_motifs = match (cfg.paired_end, fasta_path) {
            (true, Some(p)) => Some(EndMotifAccum::new(p)?),
            _ => None,
        };
        let fragment_size = if cfg.paired_end {
            Some(FragmentSizeAccum::new())
        } else {
            None
        };
        Ok(Self {
            end_motifs,
            soft_clips: Some(SoftClipAccum::new()),
            fragment_size,
            cfg: *cfg,
        })
    }

    /// Per-record hot path. `chrom` is the contig name as it appears in the
    /// BAM header; required by [`EndMotifAccum`] for FASTA lookups.
    pub fn process_read(&mut self, record: &Record, chrom: &str) {
        if let Some(em) = self.end_motifs.as_mut() {
            em.process_read(record, chrom, self.cfg.mapq_cut);
        }
        if let Some(sc) = self.soft_clips.as_mut() {
            sc.process_read(record, self.cfg.mapq_cut);
        }
        if let Some(fs) = self.fragment_size.as_mut() {
            fs.process_read(record, self.cfg.mapq_cut);
        }
    }

    /// Merge a peer worker into self. Per-feature merges are additive.
    pub fn merge(&mut self, other: FragmentomicsAccumulators) {
        match (self.end_motifs.as_mut(), other.end_motifs) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.end_motifs = Some(b),
            _ => {}
        }
        match (self.soft_clips.as_mut(), other.soft_clips) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.soft_clips = Some(b),
            _ => {}
        }
        match (self.fragment_size.as_mut(), other.fragment_size) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.fragment_size = Some(b),
            _ => {}
        }
    }

    #[allow(dead_code)]
    pub fn config(&self) -> &FragmentomicsConfig {
        &self.cfg
    }
}
