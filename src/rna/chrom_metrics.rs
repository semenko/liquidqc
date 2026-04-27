//! Per-chromosome read counts + proper-pair insert-size mean.
//!
//! Surfaces signals implicit in `idxstats` per-contig: chrM (necrosis proxy),
//! chrY (sex / transplant chimerism), CNV bias via per-chromosome insert-size
//! mean. Tid → contig name resolved at finalize from the BAM header.

use rust_htslib::bam::record::Record;
use serde::Serialize;
use std::collections::BTreeMap;

use crate::rna::bam_flags::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FSECONDARY,
    BAM_FSUPPLEMENTARY, BAM_FUNMAP,
};

#[derive(Debug, Default, Clone)]
struct TidCounts {
    total_records: u64,
    mapped: u64,
    duplicates: u64,
    /// Proper pairs counted on the leftmost mate only (`insert_size > 0`).
    proper_pair_count: u64,
    /// Sum of `|insert_size|` over the same population.
    proper_pair_insert_size_sum: u128,
}

impl TidCounts {
    fn merge(&mut self, other: TidCounts) {
        self.total_records = self.total_records.saturating_add(other.total_records);
        self.mapped = self.mapped.saturating_add(other.mapped);
        self.duplicates = self.duplicates.saturating_add(other.duplicates);
        self.proper_pair_count = self
            .proper_pair_count
            .saturating_add(other.proper_pair_count);
        self.proper_pair_insert_size_sum = self
            .proper_pair_insert_size_sum
            .saturating_add(other.proper_pair_insert_size_sum);
    }
}

/// Per-worker accumulator. Maintained as a `Vec<TidCounts>` indexed by BAM tid;
/// the BAM header reference count is passed at construction so the inner vec
/// can be allocated zeroed-out and indexed without bounds-checking on every
/// record.
#[derive(Debug, Clone)]
pub struct ChromMetricsAccum {
    counts: Vec<TidCounts>,
}

impl ChromMetricsAccum {
    pub fn new(num_tids: usize) -> Self {
        Self {
            counts: vec![TidCounts::default(); num_tids],
        }
    }

    /// Per-record hot path. No-op for unmapped (`tid() < 0`) records or
    /// records past the header reference range (defensive).
    pub fn process_read(&mut self, record: &Record) {
        let tid = record.tid();
        if tid < 0 {
            return;
        }
        let idx = tid as usize;
        if idx >= self.counts.len() {
            return;
        }
        let c = &mut self.counts[idx];
        c.total_records = c.total_records.saturating_add(1);

        let flags = record.flags();
        let is_unmapped = flags & BAM_FUNMAP != 0;
        let is_secondary = flags & BAM_FSECONDARY != 0;
        let is_suppl = flags & BAM_FSUPPLEMENTARY != 0;
        let is_qc_fail = flags & BAM_FQCFAIL != 0;
        let is_dup = flags & BAM_FDUP != 0;

        if !is_unmapped {
            c.mapped = c.mapped.saturating_add(1);
        }
        if is_dup {
            c.duplicates = c.duplicates.saturating_add(1);
        }

        // Proper-pair insert-size: leftmost mate of a properly paired,
        // primary, on-same-chrom alignment. Mirrors `is_leftmost_proper_pair_mate`
        // in fragmentomics::common but with no MAPQ filter (we want a chrom-
        // level signal that includes the lower-MAPQ tail).
        if is_unmapped
            || is_secondary
            || is_suppl
            || is_qc_fail
            || flags & BAM_FMUNMAP != 0
            || flags & BAM_FPAIRED == 0
            || flags & BAM_FPROPER_PAIR == 0
        {
            return;
        }
        if record.tid() != record.mtid() {
            return;
        }
        let isize = record.insert_size();
        if isize <= 0 {
            return;
        }
        c.proper_pair_count = c.proper_pair_count.saturating_add(1);
        c.proper_pair_insert_size_sum = c
            .proper_pair_insert_size_sum
            .saturating_add(isize.unsigned_abs() as u128);
    }

    pub fn merge(&mut self, other: ChromMetricsAccum) {
        if self.counts.len() < other.counts.len() {
            self.counts.resize(other.counts.len(), TidCounts::default());
        }
        for (i, oc) in other.counts.into_iter().enumerate() {
            self.counts[i].merge(oc);
        }
    }

    /// Resolve tid → contig name and contig length, drop entries with no
    /// records, and emit a `BTreeMap<contig, …>` keyed by header order.
    pub fn into_result(self, bam_header_refs: &[(String, u64)]) -> ChromMetricsResult {
        let mut total_mapped: u64 = 0;
        for c in &self.counts {
            total_mapped = total_mapped.saturating_add(c.mapped);
        }

        let mut entries = BTreeMap::new();
        for (i, c) in self.counts.into_iter().enumerate() {
            if c.total_records == 0 {
                continue;
            }
            let (name, length) = match bam_header_refs.get(i) {
                Some((n, l)) => (n.clone(), *l),
                None => continue,
            };
            let proper_pair_insert_size_mean = if c.proper_pair_count > 0 {
                (c.proper_pair_insert_size_sum as f64) / (c.proper_pair_count as f64)
            } else {
                0.0
            };
            let read_fraction = if total_mapped > 0 {
                c.mapped as f64 / total_mapped as f64
            } else {
                0.0
            };
            entries.insert(
                name,
                ChromMetricsEntry {
                    contig_length: length,
                    total_records: c.total_records,
                    mapped: c.mapped,
                    duplicates: c.duplicates,
                    read_fraction,
                    proper_pair_count: c.proper_pair_count,
                    proper_pair_insert_size_mean,
                },
            );
        }
        ChromMetricsResult {
            total_mapped,
            contigs: entries,
        }
    }
}

/// Per-contig aggregate, schema-stable.
#[derive(Debug, Clone, Serialize)]
pub struct ChromMetricsEntry {
    pub contig_length: u64,
    pub total_records: u64,
    pub mapped: u64,
    pub duplicates: u64,
    /// Fraction of mapped reads landing on this contig.
    pub read_fraction: f64,
    pub proper_pair_count: u64,
    /// Mean `|insert_size|` over the proper-pair leftmost-mate population
    /// on this contig. `0.0` when `proper_pair_count == 0`.
    pub proper_pair_insert_size_mean: f64,
}

/// Sample-level result.
#[derive(Debug, Clone, Serialize)]
pub struct ChromMetricsResult {
    pub total_mapped: u64,
    pub contigs: BTreeMap<String, ChromMetricsEntry>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{CigarString, Record};

    fn build_record(tid: i32, flags: u16, isize: i64, mtid: i32) -> Record {
        let mut record = Record::new();
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(75)]);
        let bases = vec![b'A'; 75];
        let quals = vec![20u8; 75];
        record.set(b"r", Some(&cigar), &bases, &quals);
        record.set_flags(flags);
        record.set_tid(tid);
        record.set_mtid(mtid);
        record.set_pos(0);
        record.set_mpos(0);
        record.set_insert_size(isize);
        record
    }

    #[test]
    fn unmapped_records_are_skipped() {
        let mut a = ChromMetricsAccum::new(3);
        let r = build_record(-1, BAM_FUNMAP, 0, -1);
        a.process_read(&r);
        for c in &a.counts {
            assert_eq!(c.total_records, 0);
        }
    }

    #[test]
    fn mapped_record_counts_against_its_contig() {
        let mut a = ChromMetricsAccum::new(3);
        let r = build_record(1, 0, 0, 1);
        a.process_read(&r);
        assert_eq!(a.counts[1].total_records, 1);
        assert_eq!(a.counts[1].mapped, 1);
        assert_eq!(a.counts[0].total_records, 0);
    }

    #[test]
    fn duplicate_flag_is_counted_separately() {
        let mut a = ChromMetricsAccum::new(2);
        let r = build_record(0, BAM_FDUP, 0, 0);
        a.process_read(&r);
        assert_eq!(a.counts[0].mapped, 1);
        assert_eq!(a.counts[0].duplicates, 1);
    }

    #[test]
    fn proper_pair_isize_only_counted_on_leftmost_mate() {
        let mut a = ChromMetricsAccum::new(2);
        // Leftmost mate: insert_size > 0 → counted.
        let r1 = build_record(0, BAM_FPAIRED | BAM_FPROPER_PAIR, 200, 0);
        a.process_read(&r1);
        // Rightmost mate: insert_size < 0 → skipped (avoids double counting).
        let r2 = build_record(0, BAM_FPAIRED | BAM_FPROPER_PAIR, -200, 0);
        a.process_read(&r2);
        assert_eq!(a.counts[0].proper_pair_count, 1);
        assert_eq!(a.counts[0].proper_pair_insert_size_sum, 200);
    }

    #[test]
    fn merge_is_additive() {
        let mut a = ChromMetricsAccum::new(2);
        let mut b = ChromMetricsAccum::new(2);
        a.counts[0].mapped = 5;
        a.counts[0].total_records = 5;
        a.counts[0].proper_pair_count = 2;
        a.counts[0].proper_pair_insert_size_sum = 400;
        b.counts[0].mapped = 3;
        b.counts[0].total_records = 3;
        b.counts[0].proper_pair_count = 1;
        b.counts[0].proper_pair_insert_size_sum = 100;
        a.merge(b);
        assert_eq!(a.counts[0].mapped, 8);
        assert_eq!(a.counts[0].proper_pair_count, 3);
        assert_eq!(a.counts[0].proper_pair_insert_size_sum, 500);
    }

    #[test]
    fn into_result_drops_zero_record_contigs_and_resolves_names() {
        let mut a = ChromMetricsAccum::new(3);
        a.counts[0].total_records = 10;
        a.counts[0].mapped = 9;
        a.counts[0].proper_pair_count = 5;
        a.counts[0].proper_pair_insert_size_sum = 1000;
        a.counts[2].total_records = 1;
        a.counts[2].mapped = 1;
        let header = vec![
            ("chr1".to_string(), 1000u64),
            ("chr2".to_string(), 2000u64),
            ("chrM".to_string(), 16569u64),
        ];
        let r = a.into_result(&header);
        assert_eq!(r.total_mapped, 10);
        assert!(r.contigs.contains_key("chr1"));
        assert!(!r.contigs.contains_key("chr2"));
        assert!(r.contigs.contains_key("chrM"));
        let c1 = &r.contigs["chr1"];
        assert_eq!(c1.contig_length, 1000);
        assert_eq!(c1.proper_pair_insert_size_mean, 200.0);
        assert!((c1.read_fraction - 0.9).abs() < 1e-12);
    }
}
