//! Per-record sample-level accumulators that ride the single BAM pass.
//!
//! Hooks into the per-chromosome BAM dispatcher in
//! [`crate::rna::dupradar::counting::count_chromosome`]. Each worker owns one
//! [`FragmentomicsAccumulators`] alongside its
//! [`crate::rna::rseqc::accumulators::RseqcAccumulators`]; per-record dispatch
//! is a single extra method call per record. Workers are merged at the same
//! point as the RSeQC accumulators.
//!
//! Originally introduced for the Phase 2 Tier 1 fragmentomics features
//! ([`end_motifs`], [`soft_clips`], [`fragment_size`], [`periodicity`]). The
//! container also carries the Phase 4 Tier 3 single-pass accumulators
//! (per-chromosome metrics in [`crate::rna::chrom_metrics`] and per-cycle base
//! quality in [`crate::rna::cycle_quality`]) so they share the same dispatch
//! path; their submodules live outside this `fragmentomics` namespace because
//! they aren't fragmentomics features.
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

use crate::rna::chrom_metrics::ChromMetricsAccum;
use crate::rna::cycle_quality::CycleQualityAccum;

/// Run-time inputs for the per-record sample-level accumulators.
///
/// Each accumulator is automatically inhibited when its prerequisite is
/// missing: end-motifs require `--fasta`; end-motifs and fragment-size
/// require paired-end. Soft-clips, per-chromosome metrics, and per-cycle
/// quality all run on any input.
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
    /// Phase 4 Tier 3: per-chromosome read counts + proper-pair |TLEN| mean.
    pub chrom_metrics: Option<ChromMetricsAccum>,
    /// Phase 4 Tier 3: per-cycle base-quality drop-off (R1 / R2).
    pub cycle_quality: Option<CycleQualityAccum>,
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
    /// accumulators (end-motifs, fragment-size). Soft-clips, per-chromosome
    /// metrics, and per-cycle quality still run.
    ///
    /// `num_tids` is the BAM header reference count, used to size the
    /// per-chromosome counter vector inside `ChromMetricsAccum`. The header
    /// is only available inside `count_reads` (the caller resolves it from
    /// the open `bam::Reader`), so this can't be stashed in the config.
    pub fn new(
        cfg: &FragmentomicsConfig,
        fasta_path: Option<&Path>,
        num_tids: usize,
    ) -> Result<Self> {
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
            chrom_metrics: Some(ChromMetricsAccum::new(num_tids)),
            cycle_quality: Some(CycleQualityAccum::new()),
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
        if let Some(cm) = self.chrom_metrics.as_mut() {
            cm.process_read(record);
        }
        if let Some(cq) = self.cycle_quality.as_mut() {
            cq.process_read(record);
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
        match (self.chrom_metrics.as_mut(), other.chrom_metrics) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.chrom_metrics = Some(b),
            _ => {}
        }
        match (self.cycle_quality.as_mut(), other.cycle_quality) {
            (Some(a), Some(b)) => a.merge(b),
            (None, Some(b)) => self.cycle_quality = Some(b),
            _ => {}
        }
    }

    #[allow(dead_code)]
    pub fn config(&self) -> &FragmentomicsConfig {
        &self.cfg
    }
}
