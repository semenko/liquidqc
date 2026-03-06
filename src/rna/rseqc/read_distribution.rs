//! Read distribution across genomic features.
//!
//! Reimplementation of RSeQC's `read_distribution.py`: classifies BAM read tags
//! into CDS exons, 5'/3' UTRs, introns, and intergenic regions using a BED12
//! gene model. Tags (CIGAR M-blocks) are classified by midpoint with priority
//! CDS > UTR > Intron > Intergenic.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;

use crate::gtf::Gene;

// ===================================================================
// Region types
// ===================================================================

/// A genomic interval: [start, end) on a chromosome (0-based, half-open).
#[derive(Debug, Clone, Copy)]
struct Interval {
    start: u64,
    end: u64,
}

/// Merged intervals for a single chromosome, sorted and non-overlapping.
#[derive(Debug, Clone, Default)]
pub struct ChromIntervals {
    intervals: Vec<Interval>,
}

impl ChromIntervals {
    /// Add a raw interval (may overlap existing ones).
    fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push(Interval { start, end });
        }
    }

    /// Sort and merge overlapping intervals.
    fn merge(&mut self) {
        if self.intervals.is_empty() {
            return;
        }
        self.intervals.sort_by_key(|i| (i.start, i.end));
        let mut merged = Vec::with_capacity(self.intervals.len());
        let mut current = self.intervals[0];
        for iv in &self.intervals[1..] {
            if iv.start <= current.end {
                current.end = current.end.max(iv.end);
            } else {
                merged.push(current);
                current = *iv;
            }
        }
        merged.push(current);
        self.intervals = merged;
    }

    /// Subtract another set of intervals from this one.
    /// Both must be pre-merged and sorted.
    fn subtract(&mut self, other: &ChromIntervals) {
        if self.intervals.is_empty() || other.intervals.is_empty() {
            return;
        }
        let mut result = Vec::new();
        let mut j = 0;
        for iv in &self.intervals {
            let mut start = iv.start;
            let end = iv.end;
            while j < other.intervals.len() && other.intervals[j].end <= start {
                j += 1;
            }
            let mut k = j;
            while k < other.intervals.len() && other.intervals[k].start < end {
                let sub = &other.intervals[k];
                if sub.start > start {
                    result.push(Interval {
                        start,
                        end: sub.start.min(end),
                    });
                }
                start = sub.end;
                k += 1;
            }
            if start < end {
                result.push(Interval { start, end });
            }
        }
        self.intervals = result;
    }

    /// Total bases covered by these intervals.
    pub fn total_bases(&self) -> u64 {
        self.intervals.iter().map(|i| i.end - i.start).sum()
    }

    /// Check if a point falls strictly inside any interval (binary search).
    /// Uses strict containment: `start < point < end`, matching RSeQC's
    /// bx-python `Intersecter.find(mid, mid)` which checks
    /// `(self.start < end) and (self.end > start)` with start==end==mid,
    /// i.e. `self.start < mid and self.end > mid`.
    pub fn contains(&self, point: u64) -> bool {
        self.intervals
            .binary_search_by(|iv| {
                if point <= iv.start {
                    std::cmp::Ordering::Greater
                } else if point >= iv.end {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .is_ok()
    }
}

/// Per-chromosome region sets for all feature types.
///
/// This is the core data structure for read distribution analysis. It can be
/// built from either a BED12 file or from GTF gene annotations.
#[derive(Debug, Default)]
pub struct RegionSets {
    pub cds_exon: HashMap<String, ChromIntervals>,
    pub utr_5: HashMap<String, ChromIntervals>,
    pub utr_3: HashMap<String, ChromIntervals>,
    pub intron: HashMap<String, ChromIntervals>,
    pub tss_up_1kb: HashMap<String, ChromIntervals>,
    pub tss_up_5kb: HashMap<String, ChromIntervals>,
    pub tss_up_10kb: HashMap<String, ChromIntervals>,
    pub tes_down_1kb: HashMap<String, ChromIntervals>,
    pub tes_down_5kb: HashMap<String, ChromIntervals>,
    pub tes_down_10kb: HashMap<String, ChromIntervals>,
}

// ===================================================================
// BED12 parsing and region extraction
// ===================================================================

/// A parsed BED12 transcript record.
#[derive(Debug)]
struct Bed12Record {
    chrom: String,
    tx_start: u64,
    tx_end: u64,
    strand: char,
    cds_start: u64,
    cds_end: u64,
    exon_starts: Vec<u64>,
    exon_ends: Vec<u64>,
}

impl Bed12Record {
    /// Parse a BED12 line. Returns None if malformed.
    fn parse(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }
        let chrom = fields[0].to_uppercase();
        let tx_start: u64 = fields[1].parse().ok()?;
        let tx_end: u64 = fields[2].parse().ok()?;
        let strand = fields[5].chars().next().unwrap_or('+');
        let cds_start: u64 = fields[6].parse().ok()?;
        let cds_end: u64 = fields[7].parse().ok()?;
        let block_count: usize = fields[9].parse().ok()?;
        let block_sizes: Vec<u64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();
        let block_starts: Vec<u64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        if block_sizes.len() < block_count || block_starts.len() < block_count {
            return None;
        }

        let exon_starts: Vec<u64> = block_starts
            .iter()
            .take(block_count)
            .map(|s| tx_start + s)
            .collect();
        let exon_ends: Vec<u64> = exon_starts
            .iter()
            .zip(block_sizes.iter().take(block_count))
            .map(|(s, sz)| s + sz)
            .collect();

        Some(Bed12Record {
            chrom,
            tx_start,
            tx_end,
            strand,
            cds_start,
            cds_end,
            exon_starts,
            exon_ends,
        })
    }

    /// Get CDS exon intervals (intersection of exon blocks with [cds_start, cds_end)).
    fn get_cds_exons(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            if *end <= self.cds_start || *start >= self.cds_end {
                continue;
            }
            let s = (*start).max(self.cds_start);
            let e = (*end).min(self.cds_end);
            if s < e {
                result.push((s, e));
            }
        }
        result
    }

    /// Get 5' UTR intervals.
    fn get_utr_5(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            match self.strand {
                '+' => {
                    if *start < self.cds_start {
                        result.push((*start, (*end).min(self.cds_start)));
                    }
                }
                '-' => {
                    if *end > self.cds_end {
                        result.push(((*start).max(self.cds_end), *end));
                    }
                }
                _ => {}
            }
        }
        result
    }

    /// Get 3' UTR intervals.
    fn get_utr_3(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            match self.strand {
                '+' => {
                    if *end > self.cds_end {
                        result.push(((*start).max(self.cds_end), *end));
                    }
                }
                '-' => {
                    if *start < self.cds_start {
                        result.push((*start, (*end).min(self.cds_start)));
                    }
                }
                _ => {}
            }
        }
        result
    }

    /// Get intron intervals (gaps between exon blocks).
    fn get_introns(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for i in 0..self.exon_starts.len().saturating_sub(1) {
            let start = self.exon_ends[i];
            let end = self.exon_starts[i + 1];
            if start < end {
                result.push((start, end));
            }
        }
        result
    }

    /// Get TSS upstream region.
    fn get_tss_upstream(&self, size: u64) -> (u64, u64) {
        match self.strand {
            '-' => (self.tx_end, self.tx_end + size),
            _ => (self.tx_start.saturating_sub(size), self.tx_start),
        }
    }

    /// Get TES downstream region.
    fn get_tes_downstream(&self, size: u64) -> (u64, u64) {
        match self.strand {
            '-' => (self.tx_start.saturating_sub(size), self.tx_start),
            _ => (self.tx_end, self.tx_end + size),
        }
    }
}

// ===================================================================
// Build region sets from BED12 file
// ===================================================================

/// Build all region sets from a BED12 gene model file.
pub fn build_regions_from_bed(bed_path: &str) -> Result<RegionSets> {
    let content = crate::io::read_to_string(bed_path)
        .with_context(|| format!("Failed to read BED file: {}", bed_path))?;

    let mut regions = RegionSets::default();

    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }
        let rec = match Bed12Record::parse(line) {
            Some(r) => r,
            None => continue,
        };

        // CDS exons
        for (s, e) in rec.get_cds_exons() {
            regions
                .cds_exon
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // 5' UTR
        for (s, e) in rec.get_utr_5() {
            regions
                .utr_5
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // 3' UTR
        for (s, e) in rec.get_utr_3() {
            regions
                .utr_3
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // Introns
        for (s, e) in rec.get_introns() {
            regions
                .intron
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // TSS upstream regions
        for (size, map) in [
            (1000u64, &mut regions.tss_up_1kb),
            (5000, &mut regions.tss_up_5kb),
            (10000, &mut regions.tss_up_10kb),
        ] {
            let (s, e) = rec.get_tss_upstream(size);
            map.entry(rec.chrom.clone()).or_default().add(s, e);
        }

        // TES downstream regions
        for (size, map) in [
            (1000u64, &mut regions.tes_down_1kb),
            (5000, &mut regions.tes_down_5kb),
            (10000, &mut regions.tes_down_10kb),
        ] {
            let (s, e) = rec.get_tes_downstream(size);
            map.entry(rec.chrom.clone()).or_default().add(s, e);
        }
    }

    // Merge all regions
    for map in [
        &mut regions.cds_exon,
        &mut regions.utr_5,
        &mut regions.utr_3,
        &mut regions.intron,
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        for intervals in map.values_mut() {
            intervals.merge();
        }
    }

    // Priority-based subtraction to make regions mutually exclusive
    // utr_5 -= cds_exon
    subtract_regions(&mut regions.utr_5, &regions.cds_exon);
    // utr_3 -= cds_exon
    subtract_regions(&mut regions.utr_3, &regions.cds_exon);
    // intron -= cds_exon, utr_5, utr_3
    subtract_regions(&mut regions.intron, &regions.cds_exon);
    subtract_regions(&mut regions.intron, &regions.utr_5);
    subtract_regions(&mut regions.intron, &regions.utr_3);
    // intergenic -= cds_exon, utr_5, utr_3, intron
    for intergenic in [
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        subtract_regions(intergenic, &regions.cds_exon);
        subtract_regions(intergenic, &regions.utr_5);
        subtract_regions(intergenic, &regions.utr_3);
        subtract_regions(intergenic, &regions.intron);
    }

    Ok(regions)
}

/// Build all region sets from GTF gene annotations.
///
/// Derives the same region types as [`build_regions_from_bed`] by using
/// transcript-level exon and CDS information from the GTF:
///
/// - **CDS exons**: intersection of transcript exon blocks with CDS range
/// - **5'/3' UTR**: exon portions outside CDS range (strand-aware)
/// - **Introns**: gaps between consecutive exon blocks per transcript
/// - **TSS/TES flanking**: upstream/downstream of transcript start/end (strand-aware)
///
/// Non-coding transcripts (no CDS) have a virtual CDS spanning the entire
/// transcript (`tx_start` to `tx_end`), matching the Perl gtf2bed convention
/// where `thickStart == txStart` and `thickEnd == txEnd`. This means all
/// exon bases become CDS exon bases rather than UTR.
pub fn build_regions_from_genes(genes: &IndexMap<String, Gene>) -> RegionSets {
    let mut regions = RegionSets::default();

    for gene in genes.values() {
        for tx in &gene.transcripts {
            let chrom = tx.chrom.to_uppercase();
            let strand = tx.strand;

            // Convert GTF 1-based inclusive exons to 0-based half-open
            let exon_starts: Vec<u64> = tx.exons.iter().map(|&(s, _)| s - 1).collect();
            let exon_ends: Vec<u64> = tx.exons.iter().map(|&(_, e)| e).collect(); // GTF end inclusive -> BED end exclusive = same value

            // CDS range (0-based half-open).
            // Non-coding transcripts: synthesize cds_start = tx_start and
            // cds_end = tx_end, so the entire transcript span is treated as
            // CDS-like. This means all exon bases become CDS exon bases
            // (instead of UTR), matching the Perl gtf2bed convention where
            // non-coding transcripts have thickStart=txStart, thickEnd=txEnd.
            let tx_start = tx.start - 1; // GTF 1-based -> 0-based
            let tx_end = tx.end; // GTF inclusive end -> exclusive end

            let (cds_start, cds_end) =
                if let (Some(cds_s), Some(cds_e)) = (tx.cds_start, tx.cds_end) {
                    (cds_s - 1, cds_e) // GTF 1-based inclusive -> 0-based half-open
                } else {
                    // Non-coding: treat entire transcript as CDS-like
                    // (matches Perl gtf2bed thickStart=txStart, thickEnd=txEnd)
                    (tx_start, tx_end)
                };

            // CDS exons: intersection of exon blocks with CDS range
            for (&es, &ee) in exon_starts.iter().zip(exon_ends.iter()) {
                if ee <= cds_start || es >= cds_end {
                    continue;
                }
                let s = es.max(cds_start);
                let e = ee.min(cds_end);
                if s < e {
                    regions.cds_exon.entry(chrom.clone()).or_default().add(s, e);
                }
            }

            // 5' UTR: exon portions before CDS (strand-aware)
            for (&es, &ee) in exon_starts.iter().zip(exon_ends.iter()) {
                match strand {
                    '+' => {
                        if es < cds_start {
                            let e = ee.min(cds_start);
                            regions.utr_5.entry(chrom.clone()).or_default().add(es, e);
                        }
                    }
                    '-' => {
                        if ee > cds_end {
                            let s = es.max(cds_end);
                            regions.utr_5.entry(chrom.clone()).or_default().add(s, ee);
                        }
                    }
                    _ => {}
                }
            }

            // 3' UTR: exon portions after CDS (strand-aware)
            for (&es, &ee) in exon_starts.iter().zip(exon_ends.iter()) {
                match strand {
                    '+' => {
                        if ee > cds_end {
                            let s = es.max(cds_end);
                            regions.utr_3.entry(chrom.clone()).or_default().add(s, ee);
                        }
                    }
                    '-' => {
                        if es < cds_start {
                            let e = ee.min(cds_start);
                            regions.utr_3.entry(chrom.clone()).or_default().add(es, e);
                        }
                    }
                    _ => {}
                }
            }

            // Introns: gaps between consecutive exon blocks
            for i in 0..exon_starts.len().saturating_sub(1) {
                let start = exon_ends[i];
                let end = exon_starts[i + 1];
                if start < end {
                    regions
                        .intron
                        .entry(chrom.clone())
                        .or_default()
                        .add(start, end);
                }
            }

            // TSS upstream regions (strand-aware)
            for (size, map) in [
                (1000u64, &mut regions.tss_up_1kb),
                (5000, &mut regions.tss_up_5kb),
                (10000, &mut regions.tss_up_10kb),
            ] {
                let (s, e) = match strand {
                    '-' => (tx_end, tx_end + size),
                    _ => (tx_start.saturating_sub(size), tx_start),
                };
                map.entry(chrom.clone()).or_default().add(s, e);
            }

            // TES downstream regions (strand-aware)
            for (size, map) in [
                (1000u64, &mut regions.tes_down_1kb),
                (5000, &mut regions.tes_down_5kb),
                (10000, &mut regions.tes_down_10kb),
            ] {
                let (s, e) = match strand {
                    '-' => (tx_start.saturating_sub(size), tx_start),
                    _ => (tx_end, tx_end + size),
                };
                map.entry(chrom.clone()).or_default().add(s, e);
            }
        }
    }

    // Merge all regions
    for map in [
        &mut regions.cds_exon,
        &mut regions.utr_5,
        &mut regions.utr_3,
        &mut regions.intron,
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        for intervals in map.values_mut() {
            intervals.merge();
        }
    }

    // Priority-based subtraction to make regions mutually exclusive
    subtract_regions(&mut regions.utr_5, &regions.cds_exon);
    subtract_regions(&mut regions.utr_3, &regions.cds_exon);
    subtract_regions(&mut regions.intron, &regions.cds_exon);
    subtract_regions(&mut regions.intron, &regions.utr_5);
    subtract_regions(&mut regions.intron, &regions.utr_3);
    for intergenic in [
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        subtract_regions(intergenic, &regions.cds_exon);
        subtract_regions(intergenic, &regions.utr_5);
        subtract_regions(intergenic, &regions.utr_3);
        subtract_regions(intergenic, &regions.intron);
    }

    regions
}

/// Subtract `sub` regions from `target` regions for all chromosomes.
fn subtract_regions(
    target: &mut HashMap<String, ChromIntervals>,
    sub: &HashMap<String, ChromIntervals>,
) {
    for (chrom, intervals) in target.iter_mut() {
        if let Some(sub_intervals) = sub.get(chrom) {
            intervals.subtract(sub_intervals);
        }
    }
}

/// Result of read distribution analysis.
#[derive(Debug)]
pub struct ReadDistributionResult {
    /// Total reads processed (after filtering).
    pub total_reads: u64,
    /// Total tags (CIGAR M-blocks) across all reads.
    pub total_tags: u64,
    /// Region statistics: (name, total_bases, tag_count).
    pub regions: Vec<(String, u64, u64)>,
    /// Number of tags not assigned to any region.
    pub unassigned_tags: u64,
}

/// Write read distribution results to a file in RSeQC-compatible format.
pub fn write_read_distribution(result: &ReadDistributionResult, output_path: &Path) -> Result<()> {
    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?,
    );

    let assigned = result.total_tags - result.unassigned_tags;

    // Header stats
    writeln!(writer, "{:<30}{}", "Total Reads", result.total_reads)?;
    writeln!(writer, "{:<30}{}", "Total Tags", result.total_tags)?;
    writeln!(writer, "{:<30}{}", "Total Assigned Tags", assigned)?;

    // Separator
    writeln!(writer, "{}", "=".repeat(69))?;

    // Table header
    writeln!(
        writer,
        "{:<20}{:<20}{:<20}{:<20}",
        "Group", "Total_bases", "Tag_count", "Tags/Kb"
    )?;

    // Region rows
    for (name, bases, tags) in &result.regions {
        let tags_per_kb = *tags as f64 * 1000.0 / (*bases as f64 + 1.0);
        writeln!(
            writer,
            "{:<20}{:<20}{:<20}{:<18.2}",
            name, bases, tags, tags_per_kb
        )?;
    }

    // Footer separator
    writeln!(writer, "{}", "=".repeat(69))?;

    writer.flush()?;
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chrom_intervals_merge() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(15, 25);
        ci.add(30, 40);
        ci.merge();
        assert_eq!(ci.intervals.len(), 2);
        assert_eq!(ci.intervals[0].start, 10);
        assert_eq!(ci.intervals[0].end, 25);
        assert_eq!(ci.intervals[1].start, 30);
        assert_eq!(ci.intervals[1].end, 40);
    }

    #[test]
    fn test_chrom_intervals_subtract() {
        let mut a = ChromIntervals::default();
        a.add(10, 50);
        a.merge();

        let mut b = ChromIntervals::default();
        b.add(20, 30);
        b.merge();

        a.subtract(&b);
        assert_eq!(a.intervals.len(), 2);
        assert_eq!(a.intervals[0].start, 10);
        assert_eq!(a.intervals[0].end, 20);
        assert_eq!(a.intervals[1].start, 30);
        assert_eq!(a.intervals[1].end, 50);
    }

    #[test]
    fn test_chrom_intervals_contains() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(30, 40);
        ci.merge();
        // RSeQC uses strict containment: point > start && point < end
        // (matching bx-python IntervalTree find(mid, mid) semantics)
        assert!(!ci.contains(10)); // boundary: not contained
        assert!(ci.contains(11)); // inside [10, 20)
        assert!(ci.contains(19)); // inside [10, 20)
        assert!(!ci.contains(20)); // boundary: not contained
        assert!(ci.contains(35)); // inside [30, 40)
        assert!(!ci.contains(5)); // outside
        assert!(!ci.contains(25)); // outside (gap between intervals)
        assert!(!ci.contains(30)); // boundary: not contained
        assert!(ci.contains(31)); // inside [30, 40)
        assert!(!ci.contains(40)); // boundary: not contained
    }

    #[test]
    fn test_bed12_cds_exons() {
        // Simple: one exon, fully CDS
        let line = "chr1\t100\t200\tgene1\t0\t+\t100\t200\t0\t1\t100,\t0,";
        let rec = Bed12Record::parse(line).unwrap();
        let cds = rec.get_cds_exons();
        assert_eq!(cds, vec![(100, 200)]);
    }

    #[test]
    fn test_bed12_utr() {
        // Two exons: first is 5'UTR, second is CDS on + strand
        let line = "chr1\t100\t300\tgene1\t0\t+\t200\t300\t0\t2\t50,50,\t0,150,";
        let rec = Bed12Record::parse(line).unwrap();
        let utr5 = rec.get_utr_5();
        assert_eq!(utr5, vec![(100, 150)]); // first exon: 100-150, all before CDS at 200
        let utr3 = rec.get_utr_3();
        assert!(utr3.is_empty()); // no exon portions after CDS end (300)
    }

    #[test]
    fn test_total_bases() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(30, 50);
        ci.merge();
        assert_eq!(ci.total_bases(), 30); // 10 + 20
    }

    #[test]
    fn test_noncoding_transcript_exons_as_cds() {
        use crate::gtf::{Gene, Transcript};

        // Non-coding transcript on + strand: exons at [101,200] and [301,400] (1-based inclusive)
        // tx_start=101, tx_end=400 (1-based). No CDS.
        // With gtf2bed convention: thickStart=txStart, thickEnd=txEnd →
        // all exon bases become CDS exon bases.
        let tx_plus = Transcript {
            transcript_id: "TX_NC_PLUS".to_string(),
            chrom: "chr1".to_string(),
            start: 101,
            end: 400,
            strand: '+',
            exons: vec![(101, 200), (301, 400)],
            cds_start: None,
            cds_end: None,
        };

        // Non-coding transcript on - strand: exons at [501,600] and [701,800] (1-based inclusive)
        // All exon bases should also become CDS exon bases.
        let tx_minus = Transcript {
            transcript_id: "TX_NC_MINUS".to_string(),
            chrom: "chr1".to_string(),
            start: 501,
            end: 800,
            strand: '-',
            exons: vec![(501, 600), (701, 800)],
            cds_start: None,
            cds_end: None,
        };

        let mut genes = IndexMap::new();
        genes.insert(
            "GENE_NC".to_string(),
            Gene {
                gene_id: "GENE_NC".to_string(),
                chrom: "chr1".to_string(),
                start: 101,
                end: 800,
                strand: '+',
                exons: vec![],
                effective_length: 0,
                attributes: std::collections::HashMap::new(),
                transcripts: vec![tx_plus, tx_minus],
            },
        );

        let regions = build_regions_from_genes(&genes);

        // Non-coding transcripts: all exon bases should be CDS exons (400 total)
        let cds_bases = regions.cds_exon.get("CHR1").map_or(0, |c| c.total_bases());
        assert_eq!(
            cds_bases, 400,
            "Non-coding transcript exons should all be CDS exons"
        );

        // No UTR for non-coding transcripts (CDS spans entire transcript)
        let utr5_bases = regions.utr_5.get("CHR1").map_or(0, |c| c.total_bases());
        assert_eq!(utr5_bases, 0, "Non-coding transcripts should have no 5'UTR");

        let utr3_bases = regions.utr_3.get("CHR1").map_or(0, |c| c.total_bases());
        assert_eq!(utr3_bases, 0, "Non-coding transcripts should have no 3'UTR");
    }
}
