//! Per-read accumulator structs for single-pass RSeQC integration.
//!
//! Each RSeQC tool has an accumulator that collects per-read data during the
//! main BAM counting loop. Accumulators are created per chromosome worker and
//! merged after parallel processing, just like `ChromResult`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{Hash, Hasher};

use rust_htslib::bam;

use super::bam_stat::BamStatResult;
use super::common::{self, KnownJunctionSet, ReferenceJunctions};
use super::infer_experiment::{GeneModel, InferExperimentResult};
use super::inner_distance::{
    build_histogram, ExonBitset, InnerDistanceResult, PairRecord, TranscriptTree,
};
use super::junction_annotation::{Junction, JunctionClass, JunctionResults};
use super::junction_saturation::SaturationResult;
use super::read_distribution::{ChromIntervals, ReadDistributionResult, RegionSets};
use super::read_duplication::ReadDuplicationResult;

// ===================================================================
// BAM flag constants
// ===================================================================

const BAM_FPAIRED: u16 = 0x1;
const BAM_FPROPER_PAIR: u16 = 0x2;
const BAM_FUNMAP: u16 = 0x4;
const BAM_FREVERSE: u16 = 0x10;
const BAM_FREAD1: u16 = 0x40;
const BAM_FREAD2: u16 = 0x80;
const BAM_FSECONDARY: u16 = 0x100;
const BAM_FQCFAIL: u16 = 0x200;
const BAM_FDUP: u16 = 0x400;
const BAM_FSUPPLEMENTARY: u16 = 0x800;

// ===================================================================
// Shared references to annotation data
// ===================================================================

/// Read-only annotation data shared across all chromosome workers.
///
/// Each field is `Option` — `None` when the corresponding tool is disabled
/// or the annotation source doesn't support it (e.g., BED mode has no GTF).
pub struct RseqcAnnotations<'a> {
    /// Gene model for infer_experiment.
    pub gene_model: Option<&'a GeneModel>,
    /// Reference junctions for junction_annotation.
    pub ref_junctions: Option<&'a ReferenceJunctions>,
    /// Known junction set for junction_saturation (used at result-conversion time).
    #[allow(dead_code)]
    pub known_junctions: Option<&'a KnownJunctionSet>,
    /// Genomic region sets for read_distribution.
    pub rd_regions: Option<&'a RegionSets>,
    /// Exon bitset for inner_distance.
    pub exon_bitset: Option<&'a ExonBitset>,
    /// Transcript tree for inner_distance.
    pub transcript_tree: Option<&'a TranscriptTree>,
    /// Chromosomes present in the known junction set (precomputed for fast lookup).
    pub ref_chroms: Option<&'a HashSet<String>>,
}

/// Per-tool configuration parameters.
///
/// Some fields are only used at accumulator construction time or during
/// result conversion, not during per-read processing.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RseqcConfig {
    /// MAPQ cutoff for read quality filtering.
    pub mapq_cut: u8,
    /// Maximum reads to sample for infer_experiment.
    pub infer_experiment_sample_size: u64,
    /// Minimum intron size for junction filtering.
    pub min_intron: u64,
    /// Minimum coverage for junction saturation.
    pub junction_saturation_min_coverage: u32,
    /// Sampling start percentage for junction saturation.
    pub junction_saturation_sample_start: u32,
    /// Sampling end percentage for junction saturation.
    pub junction_saturation_sample_end: u32,
    /// Sampling step percentage for junction saturation.
    pub junction_saturation_sample_step: u32,
    /// Maximum read pairs to sample for inner distance.
    pub inner_distance_sample_size: u64,
    /// Lower bound of inner distance histogram.
    pub inner_distance_lower_bound: i64,
    /// Upper bound of inner distance histogram.
    pub inner_distance_upper_bound: i64,
    /// Bin width for inner distance histogram.
    pub inner_distance_step: i64,
    /// Which tools are enabled.
    pub bam_stat_enabled: bool,
    pub infer_experiment_enabled: bool,
    pub read_duplication_enabled: bool,
    pub read_distribution_enabled: bool,
    pub junction_annotation_enabled: bool,
    pub junction_saturation_enabled: bool,
    pub inner_distance_enabled: bool,
}

// ===================================================================
// Per-tool accumulators
// ===================================================================

/// bam_stat accumulator — simple flag/MAPQ counting.
#[derive(Debug, Default)]
pub struct BamStatAccum {
    pub total_records: u64,
    pub qc_failed: u64,
    pub duplicates: u64,
    pub non_primary: u64,
    pub unmapped: u64,
    pub non_unique: u64,
    pub unique: u64,
    pub read_1: u64,
    pub read_2: u64,
    pub forward: u64,
    pub reverse: u64,
    pub splice: u64,
    pub non_splice: u64,
    pub proper_pairs: u64,
    pub proper_pair_diff_chrom: u64,
    pub mapq_distribution: BTreeMap<u8, u64>,
}

impl BamStatAccum {
    /// Process a single BAM record. Called for EVERY record (before counting filters).
    pub fn process_read(&mut self, record: &bam::Record, mapq_cut: u8) {
        let flags = record.flags();
        self.total_records += 1;

        // 1. QC-failed
        if flags & BAM_FQCFAIL != 0 {
            self.qc_failed += 1;
            return;
        }

        // 2. Duplicate
        if flags & BAM_FDUP != 0 {
            self.duplicates += 1;
            return;
        }

        // 3. Secondary (non-primary) — NOT supplementary
        if flags & BAM_FSECONDARY != 0 {
            self.non_primary += 1;
            return;
        }

        // 4. Unmapped
        if flags & BAM_FUNMAP != 0 {
            self.unmapped += 1;
            return;
        }

        // Mapped primary non-QC-fail non-dup: record MAPQ distribution
        let mapq = record.mapq();
        *self.mapq_distribution.entry(mapq).or_insert(0) += 1;

        // 5. MAPQ classification
        if mapq < mapq_cut {
            self.non_unique += 1;
            return;
        }

        // Uniquely mapped
        self.unique += 1;

        if flags & BAM_FREAD1 != 0 {
            self.read_1 += 1;
        }
        if flags & BAM_FREAD2 != 0 {
            self.read_2 += 1;
        }
        if flags & BAM_FREVERSE != 0 {
            self.reverse += 1;
        } else {
            self.forward += 1;
        }

        // Splice detection: CIGAR N operation
        let has_splice = record
            .cigar()
            .iter()
            .any(|op| matches!(op, rust_htslib::bam::record::Cigar::RefSkip(_)));
        if has_splice {
            self.splice += 1;
        } else {
            self.non_splice += 1;
        }

        // Proper pair analysis
        if flags & BAM_FPAIRED != 0 && flags & BAM_FPROPER_PAIR != 0 {
            self.proper_pairs += 1;
            if record.tid() != record.mtid() {
                self.proper_pair_diff_chrom += 1;
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: BamStatAccum) {
        self.total_records += other.total_records;
        self.qc_failed += other.qc_failed;
        self.duplicates += other.duplicates;
        self.non_primary += other.non_primary;
        self.unmapped += other.unmapped;
        self.non_unique += other.non_unique;
        self.unique += other.unique;
        self.read_1 += other.read_1;
        self.read_2 += other.read_2;
        self.forward += other.forward;
        self.reverse += other.reverse;
        self.splice += other.splice;
        self.non_splice += other.non_splice;
        self.proper_pairs += other.proper_pairs;
        self.proper_pair_diff_chrom += other.proper_pair_diff_chrom;
        for (mapq, count) in other.mapq_distribution {
            *self.mapq_distribution.entry(mapq).or_insert(0) += count;
        }
    }
}

// -------------------------------------------------------------------
// infer_experiment accumulator
// -------------------------------------------------------------------

/// infer_experiment accumulator — strand protocol inference via sampling.
#[derive(Debug, Default)]
pub struct InferExpAccum {
    /// Paired-end strand class counts (keys: "1++", "1--", "2+-", "2-+", etc.)
    pub p_strandness: HashMap<String, u64>,
    /// Single-end strand class counts (keys: "++", "--", "+-", "-+", etc.)
    pub s_strandness: HashMap<String, u64>,
    /// Number of usable reads sampled so far.
    pub count: u64,
    /// Maximum reads to sample.
    pub sample_size: u64,
}

impl InferExpAccum {
    /// Create a new accumulator with the given sample size.
    pub fn new(sample_size: u64) -> Self {
        InferExpAccum {
            sample_size,
            ..Default::default()
        }
    }

    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        model: &GeneModel,
        mapq_cut: u8,
    ) {
        // Already reached sample size
        if self.count >= self.sample_size {
            return;
        }

        let flags = record.flags();

        // Skip QC-fail, dup, secondary, unmapped, supplementary
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
            || flags & BAM_FSUPPLEMENTARY != 0
        {
            return;
        }

        if record.mapq() < mapq_cut {
            return;
        }

        let map_strand = if record.is_reverse() { '-' } else { '+' };

        // Compute query alignment length (M+I+=+X) to match RSeQC's qlen
        let read_start = record.pos() as u64;
        let qalen: u64 = record
            .cigar()
            .iter()
            .filter_map(|op| {
                use rust_htslib::bam::record::Cigar::*;
                match op {
                    Match(len) | Ins(len) | Equal(len) | Diff(len) => Some(*len as u64),
                    _ => None,
                }
            })
            .sum();
        let read_end = read_start + qalen;

        let strands = model.find_strands(chrom, read_start, read_end);
        if strands.is_empty() {
            return;
        }

        let strand_str: String = strands
            .iter()
            .map(|&s| (s as char).to_string())
            .collect::<Vec<String>>()
            .join(":");

        if record.is_paired() {
            let read_id = if record.is_first_in_template() {
                "1"
            } else {
                "2"
            };
            let key = format!("{}{}{}", read_id, map_strand, strand_str);
            *self.p_strandness.entry(key).or_insert(0) += 1;
        } else {
            let key = format!("{}{}", map_strand, strand_str);
            *self.s_strandness.entry(key).or_insert(0) += 1;
        }

        self.count += 1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: InferExpAccum) {
        for (key, count) in other.p_strandness {
            *self.p_strandness.entry(key).or_insert(0) += count;
        }
        for (key, count) in other.s_strandness {
            *self.s_strandness.entry(key).or_insert(0) += count;
        }
        self.count += other.count;
    }
}

// -------------------------------------------------------------------
// read_duplication accumulator
// -------------------------------------------------------------------

/// read_duplication accumulator — sequence-based and position-based dedup.
#[derive(Debug, Default)]
pub struct ReadDupAccum {
    /// Sequence hash → occurrence count (hash-based dedup to save memory).
    pub seq_dup: HashMap<u128, u64>,
    /// Position key → occurrence count.
    pub pos_dup: HashMap<String, u64>,
}

impl ReadDupAccum {
    /// Process a single BAM record.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, mapq_cut: u8) {
        let flags = record.flags();

        // Filter: unmapped, QC-fail. Does NOT skip dup or secondary (intentional).
        if flags & BAM_FUNMAP != 0 || flags & BAM_FQCFAIL != 0 {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        // Sequence-based: hash the uppercased sequence
        let seq = record.seq().as_bytes();
        let seq_hash = hash_sequence(&seq);
        *self.seq_dup.entry(seq_hash).or_insert(0) += 1;

        // Position-based: build key from CIGAR
        let pos = record.pos();
        let cigar = record.cigar();
        let key = build_position_key(chrom, pos, &cigar);
        *self.pos_dup.entry(key).or_insert(0) += 1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: ReadDupAccum) {
        for (hash, count) in other.seq_dup {
            *self.seq_dup.entry(hash).or_insert(0) += count;
        }
        for (key, count) in other.pos_dup {
            *self.pos_dup.entry(key).or_insert(0) += count;
        }
    }
}

/// Hash a read sequence using SipHash-128 for deduplication.
/// Collision probability for 100M sequences is negligible (~1 in 10^19).
fn hash_sequence(seq: &[u8]) -> u128 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for &b in seq {
        b.to_ascii_uppercase().hash(&mut hasher);
    }
    // DefaultHasher produces u64; extend to u128 by double-hashing with length
    let h1 = hasher.finish();
    let mut hasher2 = std::collections::hash_map::DefaultHasher::new();
    seq.len().hash(&mut hasher2);
    h1.hash(&mut hasher2);
    let h2 = hasher2.finish();
    (h1 as u128) << 64 | (h2 as u128)
}

/// Build position key matching RSeQC's `fetch_exon` + position key logic.
fn build_position_key(chrom: &str, pos: i64, cigar: &bam::record::CigarStringView) -> String {
    use rust_htslib::bam::record::Cigar;

    let mut key = format!("{}:{}:", chrom, pos);
    let mut ref_pos = pos;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end = ref_pos + *len as i64;
                key.push_str(&format!("{}-{}:", ref_pos, end));
                ref_pos = end;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                // RSeQC bug: S advances reference position
                ref_pos += *len as i64;
            }
            Cigar::Ins(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    key
}

// -------------------------------------------------------------------
// read_distribution accumulator
// -------------------------------------------------------------------

/// read_distribution accumulator — region classification counters.
#[derive(Debug, Default)]
pub struct ReadDistAccum {
    pub total_reads: u64,
    pub total_tags: u64,
    pub cds_tags: u64,
    pub utr5_tags: u64,
    pub utr3_tags: u64,
    pub intron_tags: u64,
    pub tss_1k_tags: u64,
    pub tss_5k_tags: u64,
    pub tss_10k_tags: u64,
    pub tes_1k_tags: u64,
    pub tes_5k_tags: u64,
    pub tes_10k_tags: u64,
    pub unassigned: u64,
}

impl ReadDistAccum {
    /// Process a single BAM record.
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, regions: &RegionSets) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped. No MAPQ filter.
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }

        self.total_reads += 1;

        // Extract exon blocks (RSeQC-compatible: M-only, S advances)
        let exon_blocks = fetch_exon_blocks_rseqc(record);

        for (block_start, block_end) in exon_blocks {
            self.total_tags += 1;
            let midpoint = block_start + (block_end - block_start) / 2;

            // Priority cascade: CDS > UTR > Intron > TSS/TES > unassigned
            if point_in(&regions.cds_exon, chrom_upper, midpoint) {
                self.cds_tags += 1;
            } else if point_in(&regions.utr_5, chrom_upper, midpoint) {
                self.utr5_tags += 1;
            } else if point_in(&regions.utr_3, chrom_upper, midpoint) {
                self.utr3_tags += 1;
            } else if point_in(&regions.intron, chrom_upper, midpoint) {
                self.intron_tags += 1;
            } else if point_in(&regions.tss_up_1kb, chrom_upper, midpoint) {
                self.tss_1k_tags += 1;
            } else if point_in(&regions.tes_down_1kb, chrom_upper, midpoint) {
                self.tes_1k_tags += 1;
            } else if point_in(&regions.tss_up_5kb, chrom_upper, midpoint) {
                self.tss_5k_tags += 1;
            } else if point_in(&regions.tes_down_5kb, chrom_upper, midpoint) {
                self.tes_5k_tags += 1;
            } else if point_in(&regions.tss_up_10kb, chrom_upper, midpoint) {
                self.tss_10k_tags += 1;
            } else if point_in(&regions.tes_down_10kb, chrom_upper, midpoint) {
                self.tes_10k_tags += 1;
            } else {
                self.unassigned += 1;
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: ReadDistAccum) {
        self.total_reads += other.total_reads;
        self.total_tags += other.total_tags;
        self.cds_tags += other.cds_tags;
        self.utr5_tags += other.utr5_tags;
        self.utr3_tags += other.utr3_tags;
        self.intron_tags += other.intron_tags;
        self.tss_1k_tags += other.tss_1k_tags;
        self.tss_5k_tags += other.tss_5k_tags;
        self.tss_10k_tags += other.tss_10k_tags;
        self.tes_1k_tags += other.tes_1k_tags;
        self.tes_5k_tags += other.tes_5k_tags;
        self.tes_10k_tags += other.tes_10k_tags;
        self.unassigned += other.unassigned;
    }
}

// -------------------------------------------------------------------
// junction_annotation accumulator
// -------------------------------------------------------------------

/// junction_annotation accumulator — junction classification.
#[derive(Debug, Default)]
pub struct JuncAnnotAccum {
    /// Per-junction read counts and classification.
    pub junction_counts: HashMap<Junction, (u64, JunctionClass)>,
    pub total_events: u64,
    pub known_events: u64,
    pub partial_novel_events: u64,
    pub complete_novel_events: u64,
    pub filtered_events: u64,
}

impl JuncAnnotAccum {
    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        ref_junctions: &ReferenceJunctions,
        min_intron: u64,
        mapq_cut: u8,
    ) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        let start_pos = record.pos() as u64;
        let cigar = record.cigar();
        let introns = common::fetch_introns(start_pos, cigar.as_ref());

        for (intron_start, intron_end) in introns {
            self.total_events += 1;

            let intron_size = intron_end.saturating_sub(intron_start);
            if intron_size < min_intron {
                self.filtered_events += 1;
                continue;
            }

            let class = classify_junction(chrom_upper, intron_start, intron_end, ref_junctions);

            match class {
                JunctionClass::Annotated => self.known_events += 1,
                JunctionClass::PartialNovel => self.partial_novel_events += 1,
                JunctionClass::CompleteNovel => self.complete_novel_events += 1,
            }

            let junction = Junction {
                chrom: chrom_upper.to_string(),
                intron_start,
                intron_end,
            };
            self.junction_counts
                .entry(junction)
                .and_modify(|(c, _)| *c += 1)
                .or_insert((1, class));
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: JuncAnnotAccum) {
        for (junction, (count, class)) in other.junction_counts {
            self.junction_counts
                .entry(junction)
                .and_modify(|(c, _)| *c += count)
                .or_insert((count, class));
        }
        self.total_events += other.total_events;
        self.known_events += other.known_events;
        self.partial_novel_events += other.partial_novel_events;
        self.complete_novel_events += other.complete_novel_events;
        self.filtered_events += other.filtered_events;
    }
}

/// Classify a junction (duplicated from junction_annotation to avoid circular deps).
fn classify_junction(
    chrom: &str,
    intron_start: u64,
    intron_end: u64,
    reference: &ReferenceJunctions,
) -> JunctionClass {
    let start_known = reference
        .intron_starts
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_start));
    let end_known = reference
        .intron_ends
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_end));

    match (start_known, end_known) {
        (true, true) => JunctionClass::Annotated,
        (false, false) => JunctionClass::CompleteNovel,
        _ => JunctionClass::PartialNovel,
    }
}

// -------------------------------------------------------------------
// junction_saturation accumulator
// -------------------------------------------------------------------

/// junction_saturation accumulator — collects all junction observations.
///
/// Uses packed `(chrom_upper, start, end)` tuples instead of heap-allocated
/// Strings. The shuffle + subsampling phase runs after all workers merge.
#[derive(Debug, Default)]
pub struct JuncSatAccum {
    /// All junction observation keys: `"CHROM:start-end"`.
    pub observations: Vec<String>,
}

impl JuncSatAccum {
    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        ref_chroms: &HashSet<String>,
        min_intron: u64,
        mapq_cut: u8,
    ) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        // Skip chromosomes not in reference
        if !ref_chroms.contains(chrom_upper) {
            return;
        }

        let start = record.pos() as u64;
        let cigar = record.cigar();
        let introns = common::fetch_introns(start, cigar.as_ref());

        for (istart, iend) in introns {
            if iend - istart < min_intron {
                continue;
            }
            let key = format!("{}:{}-{}", chrom_upper, istart, iend);
            self.observations.push(key);
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: JuncSatAccum) {
        self.observations.extend(other.observations);
    }
}

// -------------------------------------------------------------------
// inner_distance accumulator
// -------------------------------------------------------------------

/// A single read pair's inner distance record (same as inner_distance::PairRecord).
#[derive(Debug)]
pub struct InnerDistPair {
    /// Read name.
    pub name: String,
    /// Inner distance (None if different chromosomes).
    pub distance: Option<i64>,
    /// Classification string.
    pub classification: String,
}

/// inner_distance accumulator — paired-end inner distance sampling.
#[derive(Debug, Default)]
pub struct InnerDistAccum {
    /// Per-pair detail records.
    pub pairs: Vec<InnerDistPair>,
    /// Distances for histogram building.
    pub distances: Vec<i64>,
    /// Number of pairs processed so far.
    pub pair_num: u64,
    /// Maximum pairs to sample.
    pub sample_size: u64,
}

impl InnerDistAccum {
    /// Create a new accumulator with the given sample size.
    pub fn new(sample_size: u64) -> Self {
        InnerDistAccum {
            sample_size,
            ..Default::default()
        }
    }

    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        exon_bitset: &ExonBitset,
        transcript_tree: &TranscriptTree,
        mapq_cut: u8,
    ) {
        // Already reached sample size
        if self.pair_num >= self.sample_size {
            return;
        }

        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, unpaired, mate unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if flags & BAM_FPAIRED == 0 || record.is_mate_unmapped() {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        let read1_start = record.pos() as u64;
        let read2_start = record.mpos() as u64;

        // Mate dedup: skip if mate has lower position (already processed)
        if read2_start < read1_start {
            return;
        }
        // Same position: skip if this is read1
        if read2_start == read1_start && record.is_first_in_template() {
            return;
        }

        self.pair_num += 1;

        let read_name = String::from_utf8_lossy(record.qname()).to_string();

        // Different chromosomes
        if record.tid() != record.mtid() {
            self.pairs.push(InnerDistPair {
                name: read_name,
                distance: None,
                classification: "sameChrom=No".to_string(),
            });
            return;
        }

        // Compute read1_end from CIGAR
        let read1_end = record.cigar().end_pos() as u64;

        // Compute inner distance
        let inner_dist: i64 = if read2_start >= read1_end {
            (read2_start - read1_end) as i64
        } else {
            // Overlap: count exonic positions of read1 in overlap region
            let exon_blocks = fetch_exon_blocks_rseqc(record);
            let mut overlap_count: i64 = 0;
            for (ex_start, ex_end) in &exon_blocks {
                let ov_start = (*ex_start).max(read2_start);
                let ov_end = (*ex_end).min(read1_end);
                if ov_start < ov_end {
                    overlap_count += (ov_end - ov_start) as i64;
                }
            }
            -overlap_count
        };

        // Check transcript membership
        let read1_genes =
            transcript_tree.find_overlapping(chrom_upper, read1_end.saturating_sub(1));
        let read2_genes = transcript_tree.find_overlapping(chrom_upper, read2_start);
        let common_genes: HashSet<_> = read1_genes.intersection(&read2_genes).collect();

        let classification: String;

        if common_genes.is_empty() {
            classification = "sameTranscript=No,dist=genomic".to_string();
        } else if inner_dist > 0 {
            if !exon_bitset.has_chrom(chrom_upper) {
                classification = "unknownChromosome,dist=genomic".to_string();
            } else {
                let exonic_bases =
                    exon_bitset.count_exonic_bases(chrom_upper, read1_end, read2_start);

                if exonic_bases as i64 == inner_dist {
                    classification = "sameTranscript=Yes,sameExon=Yes,dist=mRNA".to_string();
                } else if exonic_bases > 0 {
                    classification = "sameTranscript=Yes,sameExon=No,dist=mRNA".to_string();
                    let mrna_dist = exonic_bases as i64;
                    self.pairs.push(InnerDistPair {
                        name: read_name,
                        distance: Some(mrna_dist),
                        classification,
                    });
                    self.distances.push(mrna_dist);
                    return;
                } else {
                    classification = "sameTranscript=Yes,nonExonic=Yes,dist=genomic".to_string();
                }
            }
        } else {
            classification = "readPairOverlap".to_string();
        }

        self.pairs.push(InnerDistPair {
            name: read_name,
            distance: Some(inner_dist),
            classification,
        });
        self.distances.push(inner_dist);
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: InnerDistAccum) {
        self.pairs.extend(other.pairs);
        self.distances.extend(other.distances);
        self.pair_num += other.pair_num;
    }
}

// ===================================================================
// Shared CIGAR helper
// ===================================================================

/// Extract exon blocks from CIGAR matching RSeQC's `bam_cigar.fetch_exon()`.
///
/// Only M (Match) creates blocks. D/N/S advance reference position.
/// =/X/I/H/P are ignored. The S-advances behavior is an RSeQC bug we replicate.
fn fetch_exon_blocks_rseqc(record: &bam::Record) -> Vec<(u64, u64)> {
    let mut exons = Vec::new();
    let mut chrom_st = record.pos() as u64;

    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) => {
                let start = chrom_st;
                chrom_st += *len as u64;
                exons.push((start, chrom_st));
            }
            Del(len) | RefSkip(len) => chrom_st += *len as u64,
            SoftClip(len) => chrom_st += *len as u64, // RSeQC bug
            _ => {}
        }
    }

    exons
}

// ===================================================================
// Top-level accumulator bundle
// ===================================================================

/// Bundle of all RSeQC accumulators. Each field is `Option` — `None` when
/// that tool is disabled.
#[derive(Debug, Default)]
pub struct RseqcAccumulators {
    pub bam_stat: Option<BamStatAccum>,
    pub infer_exp: Option<InferExpAccum>,
    pub read_dup: Option<ReadDupAccum>,
    pub read_dist: Option<ReadDistAccum>,
    pub junc_annot: Option<JuncAnnotAccum>,
    pub junc_sat: Option<JuncSatAccum>,
    pub inner_dist: Option<InnerDistAccum>,
}

impl RseqcAccumulators {
    /// Create an empty set with no accumulators enabled.
    pub fn empty() -> Self {
        Self {
            bam_stat: None,
            infer_exp: None,
            read_dup: None,
            read_dist: None,
            junc_annot: None,
            junc_sat: None,
            inner_dist: None,
        }
    }

    /// Create accumulators for all enabled tools.
    pub fn new(config: &RseqcConfig) -> Self {
        RseqcAccumulators {
            bam_stat: if config.bam_stat_enabled {
                Some(BamStatAccum::default())
            } else {
                None
            },
            infer_exp: if config.infer_experiment_enabled {
                Some(InferExpAccum::new(config.infer_experiment_sample_size))
            } else {
                None
            },
            read_dup: if config.read_duplication_enabled {
                Some(ReadDupAccum::default())
            } else {
                None
            },
            read_dist: if config.read_distribution_enabled {
                Some(ReadDistAccum::default())
            } else {
                None
            },
            junc_annot: if config.junction_annotation_enabled {
                Some(JuncAnnotAccum::default())
            } else {
                None
            },
            junc_sat: if config.junction_saturation_enabled {
                Some(JuncSatAccum::default())
            } else {
                None
            },
            inner_dist: if config.inner_distance_enabled {
                Some(InnerDistAccum::new(config.inner_distance_sample_size))
            } else {
                None
            },
        }
    }

    /// Dispatch a BAM record to all enabled tool accumulators.
    ///
    /// This is called for EVERY record before the counting.rs filter cascade,
    /// so each tool applies its own filters internally.
    #[allow(clippy::too_many_arguments)]
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        chrom_upper: &str,
        annotations: &RseqcAnnotations,
        config: &RseqcConfig,
    ) {
        // bam_stat: sees all records, applies its own filters
        if let Some(ref mut accum) = self.bam_stat {
            accum.process_read(record, config.mapq_cut);
        }

        // read_duplication: needs chrom for position key, applies its own filters
        if let Some(ref mut accum) = self.read_dup {
            accum.process_read(record, chrom, config.mapq_cut);
        }

        // infer_experiment: needs gene model overlap
        if let (Some(ref mut accum), Some(model)) = (&mut self.infer_exp, annotations.gene_model) {
            accum.process_read(record, chrom, model, config.mapq_cut);
        }

        // read_distribution: needs region sets, uses uppercased chrom
        if let (Some(ref mut accum), Some(regions)) = (&mut self.read_dist, annotations.rd_regions)
        {
            accum.process_read(record, chrom_upper, regions);
        }

        // junction_annotation: needs reference junctions, uses uppercased chrom
        if let (Some(ref mut accum), Some(ref_junctions)) =
            (&mut self.junc_annot, annotations.ref_junctions)
        {
            accum.process_read(
                record,
                chrom_upper,
                ref_junctions,
                config.min_intron,
                config.mapq_cut,
            );
        }

        // junction_saturation: needs known junction chroms, uses uppercased chrom
        if let (Some(ref mut accum), Some(ref_chroms)) =
            (&mut self.junc_sat, &annotations.ref_chroms)
        {
            accum.process_read(
                record,
                chrom_upper,
                ref_chroms,
                config.min_intron,
                config.mapq_cut,
            );
        }

        // inner_distance: needs exon bitset + transcript tree, uses uppercased chrom
        if let (Some(ref mut accum), Some(exon_bitset), Some(transcript_tree)) = (
            &mut self.inner_dist,
            annotations.exon_bitset,
            annotations.transcript_tree,
        ) {
            accum.process_read(
                record,
                chrom_upper,
                exon_bitset,
                transcript_tree,
                config.mapq_cut,
            );
        }
    }

    /// Merge another set of accumulators into this one.
    pub fn merge(&mut self, other: RseqcAccumulators) {
        if let (Some(ref mut a), Some(b)) = (&mut self.bam_stat, other.bam_stat) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.infer_exp, other.infer_exp) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.read_dup, other.read_dup) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.read_dist, other.read_dist) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.junc_annot, other.junc_annot) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.junc_sat, other.junc_sat) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.inner_dist, other.inner_dist) {
            a.merge(b);
        }
    }
}

// ===================================================================
// RegionSets point-query helper
// ===================================================================

/// Check if a point falls within any interval for the given chromosome
/// in a region map (HashMap<String, ChromIntervals>).
fn point_in(region_map: &HashMap<String, ChromIntervals>, chrom: &str, point: u64) -> bool {
    region_map.get(chrom).is_some_and(|ci| ci.contains(point))
}

// ===================================================================
// Converter methods: accumulator → result types for output functions
// ===================================================================

impl BamStatAccum {
    /// Convert accumulated counters into a `BamStatResult` for output.
    pub fn into_result(self) -> BamStatResult {
        BamStatResult {
            total_records: self.total_records,
            qc_failed: self.qc_failed,
            duplicates: self.duplicates,
            non_primary: self.non_primary,
            unmapped: self.unmapped,
            non_unique: self.non_unique,
            unique: self.unique,
            read_1: self.read_1,
            read_2: self.read_2,
            forward: self.forward,
            reverse: self.reverse,
            splice: self.splice,
            non_splice: self.non_splice,
            proper_pairs: self.proper_pairs,
            proper_pair_diff_chrom: self.proper_pair_diff_chrom,
            mapq_distribution: self.mapq_distribution,
        }
    }
}

impl InferExpAccum {
    /// Convert accumulated strand counts into an `InferExperimentResult`.
    pub fn into_result(self) -> InferExperimentResult {
        let p_total: u64 = self.p_strandness.values().sum();
        let s_total: u64 = self.s_strandness.values().sum();
        let total = p_total + s_total;

        if total == 0 {
            return InferExperimentResult {
                total_sampled: 0,
                library_type: String::from("Undetermined"),
                frac_failed: 0.0,
                frac_protocol1: 0.0,
                frac_protocol2: 0.0,
            };
        }

        // PE keys for spec1: "1++", "1--", "2+-", "2-+"
        // PE keys for spec2: "1+-", "1-+", "2++", "2--"
        // SE keys for spec1: "++", "--"
        // SE keys for spec2: "+-", "-+"
        let (library_type, spec1, spec2) = if p_total > 0 && s_total > 0 {
            // Mixed PE and SE
            let pe_spec1 = *self.p_strandness.get("1++,1--,2+-,2-+").unwrap_or(&0);
            let pe_spec2 = *self.p_strandness.get("1+-,1-+,2++,2--").unwrap_or(&0);
            let se_spec1 = *self.s_strandness.get("++,--").unwrap_or(&0);
            let se_spec2 = *self.s_strandness.get("+-,-+").unwrap_or(&0);
            (
                "Mixture".to_string(),
                pe_spec1 + se_spec1,
                pe_spec2 + se_spec2,
            )
        } else if p_total > 0 {
            let spec1 = *self.p_strandness.get("1++,1--,2+-,2-+").unwrap_or(&0);
            let spec2 = *self.p_strandness.get("1+-,1-+,2++,2--").unwrap_or(&0);
            ("PairEnd".to_string(), spec1, spec2)
        } else {
            let spec1 = *self.s_strandness.get("++,--").unwrap_or(&0);
            let spec2 = *self.s_strandness.get("+-,-+").unwrap_or(&0);
            ("SingleEnd".to_string(), spec1, spec2)
        };

        let determined = spec1 + spec2;
        let failed = total - determined;
        let total_f = total as f64;

        InferExperimentResult {
            total_sampled: total,
            library_type,
            frac_failed: failed as f64 / total_f,
            frac_protocol1: spec1 as f64 / total_f,
            frac_protocol2: spec2 as f64 / total_f,
        }
    }
}

impl ReadDupAccum {
    /// Convert accumulated hash maps into a `ReadDuplicationResult`.
    ///
    /// The accumulators use u128 hash keys for sequence dedup (memory-efficient),
    /// so this just builds the duplication-level histograms from the raw counts.
    pub fn into_result(self) -> ReadDuplicationResult {
        let pos_histogram = build_dup_histogram(&self.pos_dup);
        let seq_histogram = build_dup_histogram_from_hash(&self.seq_dup);
        ReadDuplicationResult {
            pos_histogram,
            seq_histogram,
        }
    }
}

/// Build duplication-level histogram from a count map.
/// Key = duplication level, Value = number of positions/sequences at that level.
fn build_dup_histogram(counts: &HashMap<String, u64>) -> BTreeMap<u64, u64> {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

/// Build duplication-level histogram from the hash-based sequence map.
fn build_dup_histogram_from_hash(counts: &HashMap<u128, u64>) -> BTreeMap<u64, u64> {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

impl ReadDistAccum {
    /// Convert accumulated tag counters into a `ReadDistributionResult`.
    pub fn into_result(self, regions: &RegionSets) -> ReadDistributionResult {
        fn sum_bases(map: &HashMap<String, ChromIntervals>) -> u64 {
            map.values().map(|ci| ci.total_bases()).sum()
        }

        let region_data = vec![
            (
                "CDS_Exons".to_string(),
                sum_bases(&regions.cds_exon),
                self.cds_tags,
            ),
            (
                "5'UTR_Exons".to_string(),
                sum_bases(&regions.utr_5),
                self.utr5_tags,
            ),
            (
                "3'UTR_Exons".to_string(),
                sum_bases(&regions.utr_3),
                self.utr3_tags,
            ),
            (
                "Introns".to_string(),
                sum_bases(&regions.intron),
                self.intron_tags,
            ),
            (
                "TSS_up_1kb".to_string(),
                sum_bases(&regions.tss_up_1kb),
                self.tss_1k_tags,
            ),
            (
                "TSS_up_5kb".to_string(),
                sum_bases(&regions.tss_up_5kb),
                self.tss_5k_tags,
            ),
            (
                "TSS_up_10kb".to_string(),
                sum_bases(&regions.tss_up_10kb),
                self.tss_10k_tags,
            ),
            (
                "TES_down_1kb".to_string(),
                sum_bases(&regions.tes_down_1kb),
                self.tes_1k_tags,
            ),
            (
                "TES_down_5kb".to_string(),
                sum_bases(&regions.tes_down_5kb),
                self.tes_5k_tags,
            ),
            (
                "TES_down_10kb".to_string(),
                sum_bases(&regions.tes_down_10kb),
                self.tes_10k_tags,
            ),
        ];

        ReadDistributionResult {
            total_reads: self.total_reads,
            total_tags: self.total_tags,
            regions: region_data,
            unassigned_tags: self.unassigned,
        }
    }
}

impl JuncAnnotAccum {
    /// Convert accumulated junction data into `JunctionResults`.
    pub fn into_result(self) -> JunctionResults {
        JunctionResults {
            junctions: self.junction_counts,
            total_events: self.total_events,
            known_events: self.known_events,
            partial_novel_events: self.partial_novel_events,
            complete_novel_events: self.complete_novel_events,
            filtered_events: self.filtered_events,
        }
    }
}

impl JuncSatAccum {
    /// Post-process accumulated observations into a `SaturationResult`.
    ///
    /// Performs the shuffle (Phase 2) and incremental subsampling (Phase 3)
    /// that were previously done in the standalone `junction_saturation()` function.
    pub fn into_result(
        mut self,
        known_junctions: &KnownJunctionSet,
        sample_start: u32,
        sample_end: u32,
        sample_step: u32,
        min_coverage: u32,
    ) -> SaturationResult {
        use rand::seq::SliceRandom;
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;

        // Phase 2: deterministic shuffle
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        self.observations.shuffle(&mut rng);

        // Build percentage series
        let mut percentages: Vec<u32> = (sample_start..=sample_end)
            .step_by(sample_step as usize)
            .collect();
        if *percentages.last().unwrap_or(&0) != 100 {
            percentages.push(100);
        }

        // Phase 3: incremental sampling
        let total = self.observations.len();
        let mut unique_junctions: HashMap<String, u32> = HashMap::new();
        let mut prev_end = 0;
        let mut known_counts = Vec::with_capacity(percentages.len());
        let mut novel_counts = Vec::with_capacity(percentages.len());
        let mut all_counts = Vec::with_capacity(percentages.len());

        for &pct in &percentages {
            let index_end = total * pct as usize / 100;
            for obs in &self.observations[prev_end..index_end] {
                *unique_junctions.entry(obs.clone()).or_insert(0) += 1;
            }
            prev_end = index_end;

            let mut known = 0usize;
            let mut novel = 0usize;
            for (key, &count) in &unique_junctions {
                if known_junctions.junctions.contains(key) {
                    if count >= min_coverage {
                        known += 1;
                    }
                } else {
                    novel += 1;
                }
            }
            known_counts.push(known);
            novel_counts.push(novel);
            all_counts.push(unique_junctions.len());
        }

        SaturationResult {
            percentages,
            known_counts,
            novel_counts,
            all_counts,
        }
    }
}

impl InnerDistAccum {
    /// Convert accumulated pair data into an `InnerDistanceResult`.
    pub fn into_result(self, lower_bound: i64, upper_bound: i64, step: i64) -> InnerDistanceResult {
        let pairs: Vec<PairRecord> = self
            .pairs
            .into_iter()
            .map(|p| PairRecord {
                name: p.name,
                distance: p.distance,
                classification: p.classification,
            })
            .collect();
        let histogram = build_histogram(&self.distances, lower_bound, upper_bound, step);
        InnerDistanceResult {
            pairs,
            histogram,
            total_pairs: self.pair_num,
        }
    }
}
