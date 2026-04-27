//! Per-cycle base-quality drop-off, separately for read 1 and read 2.
//!
//! BAM stores SEQ/QUAL in reference orientation, so reverse-strand reads are
//! walked back-to-front to recover the sequencing-cycle index.

use rust_htslib::bam::record::Record;
use serde::Serialize;

use crate::rna::bam_flags::{BAM_FREAD2, BAM_FREVERSE};
use crate::rna::fragmentomics::common::is_primary_mapped;

/// Maximum number of sequencing cycles tracked. Reads longer than this only
/// contribute their first `MAX_CYCLES` bases (in sequencing order).
pub const MAX_CYCLES: usize = 200;

#[derive(Debug, Clone)]
struct CycleStats {
    sum: [u64; MAX_CYCLES],
    count: [u64; MAX_CYCLES],
}

impl CycleStats {
    fn new() -> Self {
        Self {
            sum: [0; MAX_CYCLES],
            count: [0; MAX_CYCLES],
        }
    }

    fn merge(&mut self, other: CycleStats) {
        for (a, b) in self.sum.iter_mut().zip(other.sum.iter()) {
            *a = a.saturating_add(*b);
        }
        for (a, b) in self.count.iter_mut().zip(other.count.iter()) {
            *a = a.saturating_add(*b);
        }
    }
}

/// Per-worker accumulator. Two `CycleStats` arrays — one for read 1 (or single-
/// end), one for read 2.
#[derive(Debug, Clone)]
pub struct CycleQualityAccum {
    r1: CycleStats,
    r2: CycleStats,
}

impl Default for CycleQualityAccum {
    fn default() -> Self {
        Self::new()
    }
}

impl CycleQualityAccum {
    pub fn new() -> Self {
        Self {
            r1: CycleStats::new(),
            r2: CycleStats::new(),
        }
    }

    pub fn process_read(&mut self, record: &Record) {
        if !is_primary_mapped(record, 0) {
            return;
        }
        let quals = record.qual();
        if quals.is_empty() {
            return;
        }
        let qlen = quals.len();
        let flags = record.flags();
        let is_reverse = flags & BAM_FREVERSE != 0;
        let is_r2 = flags & BAM_FREAD2 != 0;
        let stats = if is_r2 { &mut self.r2 } else { &mut self.r1 };
        let take = qlen.min(MAX_CYCLES);
        if is_reverse {
            // Reverse-strand reads have qual stored in reference orientation,
            // so cycle 0 (first off the sequencer) is the LAST element.
            for (cyc, q) in quals.iter().rev().take(take).enumerate() {
                stats.sum[cyc] = stats.sum[cyc].saturating_add(*q as u64);
                stats.count[cyc] = stats.count[cyc].saturating_add(1);
            }
        } else {
            for (cyc, q) in quals.iter().take(take).enumerate() {
                stats.sum[cyc] = stats.sum[cyc].saturating_add(*q as u64);
                stats.count[cyc] = stats.count[cyc].saturating_add(1);
            }
        }
    }

    pub fn merge(&mut self, other: CycleQualityAccum) {
        self.r1.merge(other.r1);
        self.r2.merge(other.r2);
    }

    /// Cycles with zero observations are truncated from the tail.
    pub fn into_result(self, paired_end: bool) -> CycleQualityResult {
        let r1 = means_truncated(&self.r1);
        let r2_block = if paired_end {
            Some(means_truncated(&self.r2))
        } else {
            None
        };
        CycleQualityResult {
            read1_mean_quality: r1,
            read2_mean_quality: r2_block,
            max_cycles_tracked: MAX_CYCLES as u32,
        }
    }
}

fn means_truncated(stats: &CycleStats) -> Vec<f64> {
    // Find the longest prefix with at least one observation.
    let mut last = 0usize;
    for (i, c) in stats.count.iter().enumerate() {
        if *c > 0 {
            last = i + 1;
        }
    }
    let mut out = Vec::with_capacity(last);
    for i in 0..last {
        let c = stats.count[i];
        let m = if c > 0 {
            stats.sum[i] as f64 / c as f64
        } else {
            0.0
        };
        out.push(m);
    }
    out
}

#[derive(Debug, Clone, Serialize)]
pub struct CycleQualityResult {
    /// Mean Q at each sequencing cycle (1 entry per cycle, in cycle order).
    /// Always present; on single-end input this carries the entire signal.
    pub read1_mean_quality: Vec<f64>,
    /// Mean Q at each sequencing cycle for read 2. Absent on single-end input.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read2_mean_quality: Option<Vec<f64>>,
    /// Cap on cycles tracked. Reads longer than this contribute only their
    /// first `max_cycles_tracked` bases (in sequencing order).
    pub max_cycles_tracked: u32,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::bam_flags::{BAM_FDUP, BAM_FSECONDARY, BAM_FUNMAP};
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    fn build(flags: u16, quals: Vec<u8>) -> Record {
        let mut record = Record::new();
        let len = quals.len() as u32;
        let cigar = CigarString(vec![Cigar::Match(len)]);
        let bases = vec![b'A'; len as usize];
        record.set(b"r", Some(&cigar), &bases, &quals);
        record.set_flags(flags);
        record
    }

    #[test]
    fn forward_strand_walked_in_order() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(0, vec![10, 20, 30, 40]));
        assert_eq!(a.r1.sum[0], 10);
        assert_eq!(a.r1.sum[3], 40);
        assert_eq!(a.r1.count[0], 1);
        assert_eq!(a.r1.count[3], 1);
        assert_eq!(a.r1.count[4], 0);
    }

    #[test]
    fn reverse_strand_walked_back_to_front() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(BAM_FREVERSE, vec![10, 20, 30, 40]));
        // Cycle 0 = last base off the sequencer in BAM = qual[3] = 40.
        assert_eq!(a.r1.sum[0], 40);
        assert_eq!(a.r1.sum[3], 10);
    }

    #[test]
    fn read2_lands_in_r2_slot() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(BAM_FREAD2, vec![15, 25]));
        assert_eq!(a.r2.sum[0], 15);
        assert_eq!(a.r2.sum[1], 25);
        assert_eq!(a.r1.count[0], 0);
    }

    #[test]
    fn long_reads_are_capped_at_max_cycles() {
        let mut a = CycleQualityAccum::new();
        let len = MAX_CYCLES + 50;
        a.process_read(&build(0, vec![20u8; len]));
        // Last tracked cycle should have exactly one observation.
        assert_eq!(a.r1.count[MAX_CYCLES - 1], 1);
        // No buffer overflow; arrays are still MAX_CYCLES long.
        assert_eq!(a.r1.sum.len(), MAX_CYCLES);
    }

    #[test]
    fn duplicates_secondary_and_unmapped_are_skipped() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(BAM_FUNMAP, vec![20]));
        a.process_read(&build(BAM_FSECONDARY, vec![20]));
        a.process_read(&build(BAM_FDUP, vec![20]));
        for c in &a.r1.count {
            assert_eq!(*c, 0);
        }
    }

    #[test]
    fn merge_is_additive() {
        let mut a = CycleQualityAccum::new();
        let mut b = CycleQualityAccum::new();
        a.process_read(&build(0, vec![10, 20]));
        b.process_read(&build(0, vec![30, 40]));
        a.merge(b);
        assert_eq!(a.r1.sum[0], 40);
        assert_eq!(a.r1.sum[1], 60);
        assert_eq!(a.r1.count[0], 2);
    }

    #[test]
    fn into_result_truncates_to_observed_prefix() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(0, vec![20, 22, 24]));
        let r = a.into_result(false);
        assert_eq!(r.read1_mean_quality.len(), 3);
        assert!((r.read1_mean_quality[0] - 20.0).abs() < 1e-12);
        assert!((r.read1_mean_quality[2] - 24.0).abs() < 1e-12);
        assert!(r.read2_mean_quality.is_none());
    }

    #[test]
    fn into_result_paired_end_emits_r2_array() {
        let mut a = CycleQualityAccum::new();
        a.process_read(&build(0, vec![20]));
        a.process_read(&build(BAM_FREAD2, vec![25]));
        let r = a.into_result(true);
        assert!(r.read2_mean_quality.is_some());
        assert_eq!(r.read2_mean_quality.unwrap()[0], 25.0);
    }
}
