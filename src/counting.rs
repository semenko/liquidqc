//! Read counting engine for BAM files.
//!
//! Assigns reads from a BAM file to genes based on GTF annotation,
//! producing four count vectors matching the dupRadar approach:
//! 1. All reads including multimappers (with duplicates)
//! 2. All reads including multimappers (without duplicates)
//! 3. Uniquely mapped reads only (with duplicates)
//! 4. Uniquely mapped reads only (without duplicates)
//!
//! This implements a simplified featureCounts-compatible counting strategy.

use crate::gtf::Gene;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

/// Flag indicating the read is a PCR or optical duplicate (0x400).
const BAM_FDUP: u16 = 0x400;
/// Flag indicating the read is unmapped (0x4).
const BAM_FUNMAP: u16 = 0x4;
/// Flag indicating the read failed quality checks (0x200).
const BAM_FQCFAIL: u16 = 0x200;
/// Flag indicating a secondary alignment (0x100).
const BAM_FSECONDARY: u16 = 0x100;
/// Flag indicating a supplementary alignment (0x800).
const BAM_FSUPPLEMENTARY: u16 = 0x800;
/// Flag indicating the read is paired (0x1).
const BAM_FPAIRED: u16 = 0x1;
/// Flag indicating it is the first read in a pair (0x40).
const BAM_FREAD1: u16 = 0x40;
/// Flag indicating the read is mapped in a proper pair (0x2).
const BAM_FPROPER_PAIR: u16 = 0x2;
/// Flag indicating the mate is unmapped (0x8).
const BAM_FMUNMAP: u16 = 0x8;
/// Flag indicating the read is reverse-complemented (0x10).
const BAM_FREVERSE: u16 = 0x10;

/// Counts for a single gene across the four counting modes.
#[derive(Debug, Clone, Default)]
pub struct GeneCounts {
    /// Count with multimappers, with duplicates
    pub all_multi: u64,
    /// Count with multimappers, without duplicates (dups excluded)
    pub nodup_multi: u64,
    /// Count without multimappers, with duplicates
    pub all_unique: u64,
    /// Count without multimappers, without duplicates (dups excluded)
    pub nodup_unique: u64,
}

/// Result of the counting step, including per-gene counts and total mapped reads.
#[derive(Debug)]
pub struct CountResult {
    /// Per-gene counts indexed by gene_id
    pub gene_counts: IndexMap<String, GeneCounts>,
    /// Total mapped reads (including multimappers, including duplicates)
    pub total_reads_multi_dup: u64,
    /// Total mapped reads (including multimappers, excluding duplicates)
    #[allow(dead_code)]
    pub total_reads_multi_nodup: u64,
    /// Total mapped reads (unique only, including duplicates)
    pub total_reads_unique_dup: u64,
    /// Total mapped reads (unique only, excluding duplicates)
    #[allow(dead_code)]
    pub total_reads_unique_nodup: u64,
}

/// An interval tree node for efficient overlap queries.
/// Uses a simple sorted-interval approach with binary search.
#[derive(Debug, Clone)]
struct GeneInterval {
    /// Start position (0-based, half-open for internal use)
    start: u64,
    /// End position (0-based, half-open)
    end: u64,
    /// Gene ID index in the genes map
    gene_id: String,
    /// Strand
    strand: char,
}

/// A simple interval index for a single chromosome.
/// Intervals are sorted by start position for binary search.
#[derive(Debug)]
struct ChromIndex {
    intervals: Vec<GeneInterval>,
}

impl ChromIndex {
    fn new(mut intervals: Vec<GeneInterval>) -> Self {
        intervals.sort_by_key(|iv| (iv.start, iv.end));
        ChromIndex { intervals }
    }

    /// Find all gene intervals overlapping the query range [start, end).
    /// Returns gene_ids of overlapping genes (may contain duplicates if
    /// multiple exons of the same gene overlap).
    fn query(&self, start: u64, end: u64) -> Vec<&GeneInterval> {
        let mut results = Vec::new();

        // Binary search for the first interval that could overlap
        // An interval overlaps [start, end) if interval.start < end AND interval.end > start
        let search_start = self
            .intervals
            .partition_point(|iv| iv.end <= start);

        for iv in &self.intervals[search_start..] {
            if iv.start >= end {
                break;
            }
            // iv.start < end (guaranteed since we didn't break)
            // iv.end > start (guaranteed since partition_point excludes iv.end <= start)
            results.push(iv);
        }

        results
    }
}

/// Build a spatial index from gene annotations for fast overlap queries.
fn build_index(genes: &IndexMap<String, Gene>) -> HashMap<String, ChromIndex> {
    let mut chrom_intervals: HashMap<String, Vec<GeneInterval>> = HashMap::new();

    for gene in genes.values() {
        // Add each exon as an interval (featureCounts assigns reads at exon level)
        for exon in &gene.exons {
            let interval = GeneInterval {
                // Convert from 1-based inclusive GTF to 0-based half-open
                start: exon.start - 1,
                end: exon.end,
                gene_id: gene.gene_id.clone(),
                strand: exon.strand,
            };
            chrom_intervals
                .entry(exon.chrom.clone())
                .or_default()
                .push(interval);
        }
    }

    chrom_intervals
        .into_iter()
        .map(|(chrom, intervals)| (chrom, ChromIndex::new(intervals)))
        .collect()
}

/// Determine if a read's strand matches the expected gene strand,
/// given the library strandedness.
///
/// # Arguments
/// * `read_reverse` - Whether the read is mapped to the reverse strand
/// * `is_read1` - Whether this is the first read in a pair (for paired-end)
/// * `paired` - Whether the library is paired-end
/// * `gene_strand` - The strand of the gene ('+' or '-')
/// * `stranded` - Library strandedness (0=unstranded, 1=forward, 2=reverse)
fn strand_matches(
    read_reverse: bool,
    is_read1: bool,
    paired: bool,
    gene_strand: char,
    stranded: u8,
) -> bool {
    if stranded == 0 || gene_strand == '.' {
        return true; // Unstranded: everything matches
    }

    // Determine the "effective" strand of the read
    // For paired-end: read1 maps to the transcript strand, read2 maps to the opposite
    // For single-end: the read maps directly
    let read_on_plus = !read_reverse;

    let effective_plus = if paired {
        if is_read1 {
            read_on_plus
        } else {
            !read_on_plus // Read2 is the complement
        }
    } else {
        read_on_plus
    };

    match stranded {
        1 => {
            // Forward stranded: read (or read1) aligns to the same strand as the gene
            (effective_plus && gene_strand == '+') || (!effective_plus && gene_strand == '-')
        }
        2 => {
            // Reverse stranded: read (or read1) aligns to the opposite strand of the gene
            (!effective_plus && gene_strand == '+') || (effective_plus && gene_strand == '-')
        }
        _ => true,
    }
}

/// Count reads from a BAM file and assign them to genes.
///
/// Performs four simultaneous counting modes matching dupRadar's approach:
/// - With/without multimappers
/// - With/without PCR duplicates
///
/// # Arguments
/// * `bam_path` - Path to the duplicate-marked BAM file
/// * `genes` - Gene annotation map from GTF parsing
/// * `stranded` - Library strandedness (0, 1, or 2)
/// * `paired` - Whether the library is paired-end
/// * `threads` - Number of threads for BAM reading
pub fn count_reads(
    bam_path: &str,
    genes: &IndexMap<String, Gene>,
    stranded: u8,
    paired: bool,
    threads: usize,
) -> Result<CountResult> {
    // Build spatial index
    let index = build_index(genes);

    // Initialize gene counts
    let mut gene_counts: IndexMap<String, GeneCounts> = IndexMap::new();
    for gene_id in genes.keys() {
        gene_counts.insert(gene_id.clone(), GeneCounts::default());
    }

    // Open BAM file
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;

    if threads > 1 {
        bam.set_threads(threads - 1)
            .context("Failed to set BAM reader threads")?;
    }

    // Get chromosome names from header
    let header = bam.header().clone();
    let tid_to_name: Vec<String> = (0..header.target_count())
        .map(|tid| {
            String::from_utf8_lossy(header.tid2name(tid)).to_string()
        })
        .collect();

    // Track statistics
    let mut total_reads: u64 = 0;
    let mut total_mapped: u64 = 0;
    let mut total_dup: u64 = 0;
    let mut total_multi: u64 = 0;

    // Totals for N calculation (mapped reads per mode)
    let mut n_multi_dup: u64 = 0;
    let mut n_multi_nodup: u64 = 0;
    let mut n_unique_dup: u64 = 0;
    let mut n_unique_nodup: u64 = 0;

    // Iterate over all records
    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result.context("Error reading BAM record")?;
        total_reads += 1;

        if total_reads % 5_000_000 == 0 {
            debug!("Processed {} reads...", total_reads);
        }

        let flags = record.flags();

        // Skip unmapped reads
        if flags & BAM_FUNMAP != 0 {
            continue;
        }

        // Skip secondary and supplementary alignments
        if flags & BAM_FSECONDARY != 0 || flags & BAM_FSUPPLEMENTARY != 0 {
            continue;
        }

        // Skip QC-failed reads
        if flags & BAM_FQCFAIL != 0 {
            continue;
        }

        // For paired-end data, apply additional filters
        if paired {
            // Must be paired
            if flags & BAM_FPAIRED == 0 {
                continue;
            }
            // Must be properly paired
            if flags & BAM_FPROPER_PAIR == 0 {
                continue;
            }
            // Mate must be mapped
            if flags & BAM_FMUNMAP != 0 {
                continue;
            }
            // Only count the first in pair to avoid double-counting
            if flags & BAM_FREAD1 == 0 {
                continue;
            }
        }

        total_mapped += 1;

        let is_dup = flags & BAM_FDUP != 0;
        if is_dup {
            total_dup += 1;
        }

        // Determine if the read is a multimapper
        // NH tag stores the number of reported alignments
        let is_multi = match record.aux(b"NH") {
            Ok(rust_htslib::bam::record::Aux::U8(nh)) => nh > 1,
            Ok(rust_htslib::bam::record::Aux::U16(nh)) => nh > 1,
            Ok(rust_htslib::bam::record::Aux::U32(nh)) => nh > 1,
            Ok(rust_htslib::bam::record::Aux::I8(nh)) => nh > 1,
            Ok(rust_htslib::bam::record::Aux::I16(nh)) => nh > 1,
            Ok(rust_htslib::bam::record::Aux::I32(nh)) => nh > 1,
            _ => false, // No NH tag = assume unique
        };
        if is_multi {
            total_multi += 1;
        }

        let is_reverse = flags & BAM_FREVERSE != 0;
        let is_read1 = flags & BAM_FREAD1 != 0;

        // Update total mapped reads per mode.
        //
        // IMPORTANT: N (total mapped reads for RPKM) must match featureCounts behavior.
        // In featureCounts, N = sum(all stats) - Unassigned_Unmapped.
        // This means N counts ALL mapped reads regardless of whether they're multimappers
        // or assigned to a feature. The multimapper/duplicate dimensions only affect
        // which reads are *counted toward genes*, not the N denominator.
        //
        // Therefore:
        // - n_multi_dup = n_unique_dup = all mapped reads (incl dups)
        // - n_multi_nodup = n_unique_nodup = all mapped non-dup reads
        //
        // The "multi" vs "unique" distinction only applies to gene assignment below.
        n_multi_dup += 1;
        n_unique_dup += 1;
        if !is_dup {
            n_multi_nodup += 1;
            n_unique_nodup += 1;
        }

        // Get the chromosome name
        let tid = record.tid();
        if tid < 0 || tid as usize >= tid_to_name.len() {
            continue;
        }
        let chrom = &tid_to_name[tid as usize];

        // Get read alignment position (0-based)
        let read_start = record.pos() as u64;
        let read_end = record.cigar().end_pos() as u64;

        // Look up overlapping genes
        let chrom_idx = match index.get(chrom) {
            Some(idx) => idx,
            None => continue, // Chromosome not in annotation
        };

        let overlaps = chrom_idx.query(read_start, read_end);

        if overlaps.is_empty() {
            continue;
        }

        // Deduplicate gene hits (a read may overlap multiple exons of the same gene)
        let mut assigned_genes: Vec<&str> = overlaps
            .iter()
            .filter(|iv| strand_matches(is_reverse, is_read1, paired, iv.strand, stranded))
            .map(|iv| iv.gene_id.as_str())
            .collect();
        assigned_genes.sort_unstable();
        assigned_genes.dedup();

        // featureCounts default: if a read overlaps multiple genes, it is ambiguous
        // and NOT assigned (unless countMultiMappingReads is set differently).
        // However, for the multimapper dimension, featureCounts still counts
        // ambiguous reads to each gene when countMultiMappingReads=TRUE.
        //
        // For simplicity and matching featureCounts defaults:
        // - If the read maps to exactly one gene: assign it
        // - If the read maps to multiple genes: skip (ambiguous)
        if assigned_genes.len() != 1 {
            continue;
        }

        let gene_id = assigned_genes[0];

        if let Some(counts) = gene_counts.get_mut(gene_id) {
            // Mode 1: Include multimappers (all reads assigned to gene)
            counts.all_multi += 1;
            if !is_dup {
                counts.nodup_multi += 1;
            }

            // Mode 2: Exclude multimappers (unique reads only)
            if !is_multi {
                counts.all_unique += 1;
                if !is_dup {
                    counts.nodup_unique += 1;
                }
            }
        }
    }

    info!(
        "Read {} total reads, {} mapped, {} duplicates, {} multimappers",
        total_reads, total_mapped, total_dup, total_multi
    );

    Ok(CountResult {
        gene_counts,
        total_reads_multi_dup: n_multi_dup,
        total_reads_multi_nodup: n_multi_nodup,
        total_reads_unique_dup: n_unique_dup,
        total_reads_unique_nodup: n_unique_nodup,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_matching_unstranded() {
        // Unstranded: everything matches
        assert!(strand_matches(false, true, false, '+', 0));
        assert!(strand_matches(true, true, false, '+', 0));
        assert!(strand_matches(false, true, false, '-', 0));
        assert!(strand_matches(true, true, false, '-', 0));
    }

    #[test]
    fn test_strand_matching_forward_single() {
        // Forward stranded, single-end
        // Read on + strand should match + gene
        assert!(strand_matches(false, true, false, '+', 1));
        // Read on - strand should match - gene
        assert!(strand_matches(true, true, false, '-', 1));
        // Read on + strand should NOT match - gene
        assert!(!strand_matches(false, true, false, '-', 1));
        // Read on - strand should NOT match + gene
        assert!(!strand_matches(true, true, false, '+', 1));
    }

    #[test]
    fn test_strand_matching_reverse_single() {
        // Reverse stranded, single-end
        // Read on + strand should match - gene (reversed)
        assert!(strand_matches(false, true, false, '-', 2));
        // Read on - strand should match + gene (reversed)
        assert!(strand_matches(true, true, false, '+', 2));
    }

    #[test]
    fn test_strand_matching_forward_paired() {
        // Forward stranded, paired-end
        // Read1 on + strand: effective + -> matches + gene
        assert!(strand_matches(false, true, true, '+', 1));
        // Read2 on + strand: effective - -> matches - gene
        assert!(strand_matches(false, false, true, '-', 1));
        // Read1 on - strand: effective - -> matches - gene
        assert!(strand_matches(true, true, true, '-', 1));
        // Read2 on - strand: effective + -> matches + gene
        assert!(strand_matches(true, false, true, '+', 1));
    }
}
