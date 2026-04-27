//! Per-sample v1 JSON envelope writer.
//!
//! Defines [`SampleEnvelope`], the canonical liquidqc v1 output schema, plus
//! one *Block view struct per inherited accumulator. Each block has a
//! `from_*` constructor that translates the inherited Debug-only result
//! types into stable, JSON-friendly fields.
//!
//! The envelope is written to `<outdir>/<sample_id>.liquidqc.json`. Schema
//! lives at `schema/v1/liquidqc.schema.json`; bump [`SCHEMA_VERSION`] in
//! lockstep with the schema file.
//!
//! Field naming and ordering mirror the schema's `required` array so a
//! reviewer can scan the struct against the schema directly.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::Path;

use crate::cli::Strandedness;
use crate::rna::dupradar::counting::CountResult;
use crate::rna::dupradar::fitting::FitResult;
use crate::rna::preseq::PreseqResult;
use crate::rna::qualimap::QualimapResult;
use crate::rna::rseqc::bam_stat::BamStatResult;
use crate::rna::rseqc::infer_experiment::InferExperimentResult;
use crate::rna::rseqc::inner_distance::InnerDistanceResult;
use crate::rna::rseqc::junction_annotation::JunctionResults;
use crate::rna::rseqc::junction_saturation::SaturationResult;
use crate::rna::rseqc::read_distribution::ReadDistributionResult;
use crate::rna::rseqc::read_duplication::ReadDuplicationResult;
use crate::rna::rseqc::tin::TinResults;

/// Schema version this struct serializes to.
///
/// Must equal the `schema_version` example string in
/// `schema/v1/liquidqc.schema.json`. Bump together when the wire format changes.
pub const SCHEMA_VERSION: &str = "0.1.0-stub";

// ============================================================================
// Top-level envelope
// ============================================================================

/// Per-sample v1 envelope. One file per BAM input.
///
/// Schema-required fields appear first, in the order listed in the schema's
/// `required` array. Optional metric blocks (covered by `additionalProperties:
/// true`) follow, each skipped from output when `None`.
#[derive(Debug, Serialize)]
pub struct SampleEnvelope {
    pub extractor_version: String,
    pub schema_version: String,
    pub git_commit: String,
    pub sample_id: String,
    pub bam_path: String,
    pub bam_md5: String,
    pub bam_size_bytes: u64,
    pub gtf_path: String,
    pub gtf_md5: String,
    pub reference_fasta_md5: Option<String>,
    pub library_prep: String,
    pub paired_end: bool,
    pub strandedness: String,
    pub read_length_max: u64,
    pub read_length_mean: f64,
    pub read_count_total: u64,
    pub read_count_after_filters: u64,
    pub filters: Filters,
    pub qc_flags: Vec<String>,
    pub runtime_seconds: f64,
    pub peak_rss_mb: f64,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub bam_stat: Option<BamStatBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub infer_experiment: Option<InferExperimentBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_distribution: Option<ReadDistributionBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_duplication: Option<ReadDuplicationBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub junction_annotation: Option<JunctionAnnotationBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub junction_saturation: Option<JunctionSaturationBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inner_distance: Option<InnerDistanceBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tin: Option<TinBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub preseq: Option<PreseqBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qualimap: Option<QualimapBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dupradar: Option<DupradarBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub featurecounts: Option<FeatureCountsBlock>,
}

/// Filters applied to the BAM read stream before metric computation.
///
/// These reflect static defaults baked into the inherited single-pass code in
/// Phase 1. Future phases can switch them on/off via CLI without changing the
/// envelope shape.
#[derive(Debug, Serialize)]
pub struct Filters {
    pub min_mapq: u8,
    pub require_proper_pair: bool,
    pub drop_secondary: bool,
    pub drop_supplementary: bool,
    pub drop_duplicates: bool,
}

impl Filters {
    /// Phase 1 filter defaults derived from the inherited counting code.
    pub fn from_mapq(mapq_cut: u8) -> Self {
        Self {
            min_mapq: mapq_cut,
            require_proper_pair: false,
            drop_secondary: true,
            drop_supplementary: true,
            drop_duplicates: false,
        }
    }
}

// ============================================================================
// Metric blocks — view structs
// ============================================================================

/// bam_stat / flagstat / samtools-stats summary (subset of [`BamStatResult`]).
///
/// Excludes the large per-cycle arrays and per-chromosome map that already
/// live in the bam_stat / idxstats / stats text outputs.
#[derive(Debug, Serialize)]
pub struct BamStatBlock {
    pub total_records: u64,
    pub primary: u64,
    pub primary_mapped: u64,
    pub primary_duplicates: u64,
    pub mapped: u64,
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_failed: u64,
    pub duplicates: u64,
    pub unique: u64,
    pub multimappers: u64,
    pub proper_pairs: u64,
    pub forward: u64,
    pub reverse: u64,
    pub splice: u64,
    pub non_splice: u64,
    pub paired_flagstat: u64,
    pub read1_flagstat: u64,
    pub read2_flagstat: u64,
    pub singletons: u64,
    pub mate_diff_chr: u64,
    pub max_read_length: u64,
    pub mean_read_length: f64,
    pub bases_mapped: u64,
}

impl BamStatBlock {
    pub fn from_result(r: &BamStatResult) -> Self {
        let mean_read_length = if r.primary_count > 0 {
            r.total_len as f64 / r.primary_count as f64
        } else {
            0.0
        };
        Self {
            total_records: r.total_records,
            primary: r.primary_count,
            primary_mapped: r.primary_mapped,
            primary_duplicates: r.primary_duplicates,
            mapped: r.mapped,
            unmapped: r.unmapped,
            secondary: r.secondary,
            supplementary: r.supplementary,
            qc_failed: r.qc_failed,
            duplicates: r.duplicates,
            unique: r.unique,
            multimappers: r.non_unique,
            proper_pairs: r.proper_pairs,
            forward: r.forward,
            reverse: r.reverse,
            splice: r.splice,
            non_splice: r.non_splice,
            paired_flagstat: r.paired_flagstat,
            read1_flagstat: r.read1_flagstat,
            read2_flagstat: r.read2_flagstat,
            singletons: r.singletons,
            mate_diff_chr: r.mate_diff_chr,
            max_read_length: r.max_len,
            mean_read_length,
            bases_mapped: r.bases_mapped,
        }
    }
}

/// Strandedness-inference summary (full [`InferExperimentResult`]).
#[derive(Debug, Serialize)]
pub struct InferExperimentBlock {
    pub total_sampled: u64,
    pub library_type: String,
    pub frac_failed: f64,
    pub frac_protocol1: f64,
    pub frac_protocol2: f64,
}

impl InferExperimentBlock {
    pub fn from_result(r: &InferExperimentResult) -> Self {
        Self {
            total_sampled: r.total_sampled,
            library_type: r.library_type.clone(),
            frac_failed: r.frac_failed,
            frac_protocol1: r.frac_protocol1,
            frac_protocol2: r.frac_protocol2,
        }
    }
}

/// Read-distribution table (all of [`ReadDistributionResult`]).
#[derive(Debug, Serialize)]
pub struct ReadDistributionRegion {
    pub name: String,
    pub total_bases: u64,
    pub tag_count: u64,
}

#[derive(Debug, Serialize)]
pub struct ReadDistributionBlock {
    pub total_reads: u64,
    pub total_tags: u64,
    pub assigned_tags: u64,
    pub unassigned_tags: u64,
    pub regions: Vec<ReadDistributionRegion>,
}

impl ReadDistributionBlock {
    pub fn from_result(r: &ReadDistributionResult) -> Self {
        Self {
            total_reads: r.total_reads,
            total_tags: r.total_tags,
            assigned_tags: r.total_tags.saturating_sub(r.unassigned_tags),
            unassigned_tags: r.unassigned_tags,
            regions: r
                .regions
                .iter()
                .map(|(name, bases, tags)| ReadDistributionRegion {
                    name: name.clone(),
                    total_bases: *bases,
                    tag_count: *tags,
                })
                .collect(),
        }
    }
}

/// Read-duplication histograms.
///
/// `BTreeMap<u64, u64>` serializes as a JSON object whose keys are
/// stringified u64s — `serde_json` handles the integer-keyed-map case
/// directly, so no intermediate map is needed.
#[derive(Debug, Serialize)]
pub struct ReadDuplicationBlock {
    pub position_histogram: BTreeMap<u64, u64>,
    pub sequence_histogram: BTreeMap<u64, u64>,
}

impl ReadDuplicationBlock {
    pub fn from_result(r: &ReadDuplicationResult) -> Self {
        Self {
            position_histogram: r.pos_histogram.clone(),
            sequence_histogram: r.seq_histogram.clone(),
        }
    }
}

/// Junction-annotation summary. The full per-junction map is not embedded —
/// that lives in the `junctions.xls` / `junctions.bed` outputs.
#[derive(Debug, Serialize)]
pub struct JunctionAnnotationBlock {
    pub total_events: u64,
    pub known_events: u64,
    pub partial_novel_events: u64,
    pub complete_novel_events: u64,
    pub filtered_events: u64,
    pub distinct_junctions_total: u64,
    pub distinct_junctions_known: u64,
    pub distinct_junctions_partial_novel: u64,
    pub distinct_junctions_novel: u64,
}

impl JunctionAnnotationBlock {
    pub fn from_result(r: &JunctionResults) -> Self {
        let counts = r.junction_counts();
        Self {
            total_events: r.total_events,
            known_events: r.known_events,
            partial_novel_events: r.partial_novel_events,
            complete_novel_events: r.complete_novel_events,
            filtered_events: r.filtered_events,
            distinct_junctions_total: counts.total,
            distinct_junctions_known: counts.known,
            distinct_junctions_partial_novel: counts.partial_novel,
            distinct_junctions_novel: counts.novel,
        }
    }
}

/// Junction-saturation curve (all four parallel vectors from [`SaturationResult`]).
#[derive(Debug, Serialize)]
pub struct JunctionSaturationBlock {
    pub percentages: Vec<u32>,
    pub known_counts: Vec<u64>,
    pub novel_counts: Vec<u64>,
    pub all_counts: Vec<u64>,
}

impl JunctionSaturationBlock {
    pub fn from_result(r: &SaturationResult) -> Self {
        Self {
            percentages: r.percentages.clone(),
            known_counts: r.known_counts.iter().map(|&v| v as u64).collect(),
            novel_counts: r.novel_counts.iter().map(|&v| v as u64).collect(),
            all_counts: r.all_counts.iter().map(|&v| v as u64).collect(),
        }
    }
}

/// Inner-distance histogram + summary stats. Per-pair detail (`pairs`) is not
/// embedded; it lives in the `inner_distance.txt` output.
#[derive(Debug, Serialize)]
pub struct InnerDistanceHistogramBin {
    pub lower: i64,
    pub upper: i64,
    pub count: u64,
}

#[derive(Debug, Serialize)]
pub struct InnerDistanceBlock {
    pub total_pairs: u64,
    pub histogram: Vec<InnerDistanceHistogramBin>,
    pub mean: Option<f64>,
    pub median: Option<f64>,
}

impl InnerDistanceBlock {
    pub fn from_result(r: &InnerDistanceResult) -> Self {
        let histogram = r
            .histogram
            .iter()
            .map(|&(lo, hi, count)| InnerDistanceHistogramBin {
                lower: lo,
                upper: hi,
                count,
            })
            .collect();
        let distances: Vec<f64> = r
            .pairs
            .iter()
            .filter_map(|p| p.distance.map(|d| d as f64))
            .collect();
        let mean = if distances.is_empty() {
            None
        } else {
            Some(distances.iter().sum::<f64>() / distances.len() as f64)
        };
        let median = if distances.is_empty() {
            None
        } else {
            Some(crate::io::median(&distances))
        };
        Self {
            total_pairs: r.total_pairs,
            histogram,
            mean,
            median,
        }
    }
}

/// TIN summary statistics over transcripts that passed the coverage threshold.
/// The per-transcript TIN table lives in `tin.txt`.
#[derive(Debug, Serialize)]
pub struct TinBlock {
    pub transcripts_total: u64,
    pub transcripts_passed: u64,
    pub mean: Option<f64>,
    pub median: Option<f64>,
    pub stdev: Option<f64>,
}

impl TinBlock {
    pub fn from_result(r: &TinResults) -> Self {
        let transcripts_passed = r.transcripts.iter().filter(|t| t.passed_threshold).count();
        if transcripts_passed == 0 {
            return Self {
                transcripts_total: r.transcripts.len() as u64,
                transcripts_passed: 0,
                mean: None,
                median: None,
                stdev: None,
            };
        }
        let (mean, median, stdev) = crate::rna::rseqc::tin::summary_stats(r);
        Self {
            transcripts_total: r.transcripts.len() as u64,
            transcripts_passed: transcripts_passed as u64,
            mean: Some(mean),
            median: Some(median),
            stdev: Some(stdev),
        }
    }
}

/// Preseq complexity-extrapolation curve.
#[derive(Debug, Serialize)]
pub struct PreseqCurvePoint {
    pub total_reads: f64,
    pub expected_distinct: f64,
    pub lower_ci: f64,
    pub upper_ci: f64,
}

#[derive(Debug, Serialize)]
pub struct PreseqBlock {
    pub curve: Vec<PreseqCurvePoint>,
}

impl PreseqBlock {
    pub fn from_result(r: &PreseqResult) -> Self {
        Self {
            curve: r
                .curve
                .iter()
                .map(|&(t, d, lo, hi)| PreseqCurvePoint {
                    total_reads: t,
                    expected_distinct: d,
                    lower_ci: lo,
                    upper_ci: hi,
                })
                .collect(),
        }
    }
}

/// Qualimap-equivalent counters from the merged result.
#[derive(Debug, Serialize)]
pub struct QualimapBlock {
    pub primary_alignments: u64,
    pub secondary_alignments: u64,
    pub not_aligned: u64,
    pub alignment_not_unique: u64,
    pub exonic_reads: u64,
    pub intronic_reads: u64,
    pub intergenic_reads: u64,
    pub overlapping_exon_reads: u64,
    pub ambiguous_reads: u64,
    pub no_feature: u64,
    pub read_count: u64,
    pub reads_at_junctions: u64,
    pub known_junction_events: u64,
    pub partly_known_junction_events: u64,
    pub left_proper_in_pair: u64,
    pub right_proper_in_pair: u64,
    pub both_proper_in_pair: u64,
}

impl QualimapBlock {
    pub fn from_result(r: &QualimapResult) -> Self {
        Self {
            primary_alignments: r.primary_alignments,
            secondary_alignments: r.secondary_alignments,
            not_aligned: r.not_aligned,
            alignment_not_unique: r.alignment_not_unique,
            exonic_reads: r.exonic_reads,
            intronic_reads: r.intronic_reads,
            intergenic_reads: r.intergenic_reads,
            overlapping_exon_reads: r.overlapping_exon_reads,
            ambiguous_reads: r.ambiguous_reads,
            no_feature: r.no_feature,
            read_count: r.read_count,
            reads_at_junctions: r.reads_at_junctions,
            known_junction_events: r.known_junction_events,
            partly_known_junction_events: r.partly_known_junction_events,
            left_proper_in_pair: r.left_proper_in_pair,
            right_proper_in_pair: r.right_proper_in_pair,
            both_proper_in_pair: r.both_proper_in_pair,
        }
    }
}

/// dupRadar fit summary.
#[derive(Debug, Serialize)]
pub struct DupradarBlock {
    pub total_genes: u64,
    pub genes_with_reads: u64,
    pub genes_with_duplication: u64,
    pub fit: Option<DupradarFit>,
}

#[derive(Debug, Serialize)]
pub struct DupradarFit {
    pub beta0: f64,
    pub beta1: f64,
    pub intercept: f64,
    pub slope: f64,
}

impl DupradarBlock {
    pub fn from_parts(
        total_genes: u64,
        genes_with_reads: u64,
        genes_with_duplication: u64,
        fit: Option<&FitResult>,
    ) -> Self {
        Self {
            total_genes,
            genes_with_reads,
            genes_with_duplication,
            fit: fit.map(|f| DupradarFit {
                beta0: f.beta0,
                beta1: f.beta1,
                intercept: f.intercept,
                slope: f.slope,
            }),
        }
    }
}

/// featureCounts assignment statistics.
#[derive(Debug, Serialize)]
pub struct FeatureCountsBlock {
    pub fragments_total: u64,
    pub fragments_assigned: u64,
    pub fragments_ambiguous: u64,
    pub fragments_no_features: u64,
    pub reads_assigned: u64,
    pub reads_ambiguous: u64,
    pub reads_no_features: u64,
    pub reads_multimapping: u64,
    pub reads_unmapped: u64,
    pub reads_singleton: u64,
    pub reads_chimera: u64,
    pub biotype_assigned: u64,
    pub biotype_ambiguous: u64,
    pub biotype_no_features: u64,
    pub biotype_counts: IndexMap<String, u64>,
}

impl FeatureCountsBlock {
    pub fn from_result(c: &CountResult) -> Self {
        let mut biotype_counts: IndexMap<String, u64> = IndexMap::new();
        for (i, name) in c.biotype_names.iter().enumerate() {
            if let Some(count) = c.biotype_reads.get(i) {
                biotype_counts.insert(name.clone(), *count);
            }
        }
        Self {
            fragments_total: c.stat_total_fragments,
            fragments_assigned: c.stat_assigned,
            fragments_ambiguous: c.stat_ambiguous,
            fragments_no_features: c.stat_no_features,
            reads_assigned: c.fc_assigned,
            reads_ambiguous: c.fc_ambiguous,
            reads_no_features: c.fc_no_features,
            reads_multimapping: c.fc_multimapping,
            reads_unmapped: c.fc_unmapped,
            reads_singleton: c.fc_singleton,
            reads_chimera: c.fc_chimera,
            biotype_assigned: c.fc_biotype_assigned,
            biotype_ambiguous: c.fc_biotype_ambiguous,
            biotype_no_features: c.fc_biotype_no_features,
            biotype_counts,
        }
    }
}

// ============================================================================
// Builder + writer
// ============================================================================

/// Inputs assembled in `process_single_bam` and handed to [`build`].
///
/// Borrowed references for everything that's already allocated; owned values
/// only for the small computed fields. Building the envelope is pure data
/// shuffling — no I/O.
pub struct BuildInputs<'a> {
    pub extractor_version: &'a str,
    pub git_commit: &'a str,
    pub sample_id: &'a str,
    pub bam_path: &'a str,
    pub bam_md5: &'a str,
    pub bam_size_bytes: u64,
    pub gtf_path: &'a str,
    pub gtf_md5: &'a str,
    pub reference_fasta_md5: Option<&'a str>,
    pub library_prep: &'a str,
    pub paired_end: bool,
    pub strandedness: Strandedness,
    pub filters: Filters,
    pub runtime_seconds: f64,
    pub peak_rss_mb: f64,
    pub qc_flags: Vec<String>,

    pub bam_stat: Option<&'a BamStatResult>,
    pub infer_experiment: Option<&'a InferExperimentResult>,
    pub read_distribution: Option<&'a ReadDistributionResult>,
    pub read_duplication: Option<&'a ReadDuplicationResult>,
    pub junction_annotation: Option<&'a JunctionResults>,
    pub junction_saturation: Option<&'a SaturationResult>,
    pub inner_distance: Option<&'a InnerDistanceResult>,
    pub tin: Option<&'a TinResults>,
    pub preseq: Option<&'a PreseqResult>,
    pub qualimap: Option<&'a QualimapResult>,
    pub count_result: Option<&'a CountResult>,
    pub dupradar_fit: Option<&'a FitResult>,
    pub dupradar_genes_with_reads: u64,
    pub dupradar_genes_with_duplication: u64,
}

/// Build a [`SampleEnvelope`] from collected accumulator results.
pub fn build(inputs: BuildInputs<'_>) -> SampleEnvelope {
    let bam_stat_block = inputs.bam_stat.map(BamStatBlock::from_result);
    let read_length_max = bam_stat_block
        .as_ref()
        .map(|b| b.max_read_length)
        .unwrap_or(0);
    let read_length_mean = bam_stat_block
        .as_ref()
        .map(|b| b.mean_read_length)
        .unwrap_or(0.0);
    let read_count_total = inputs
        .bam_stat
        .map(|b| b.total_records)
        .or_else(|| inputs.count_result.map(|c| c.stat_total_reads))
        .unwrap_or(0);
    // "after filters" in Phase 1 = primary mapped reads with MAPQ >= cutoff
    // (the fraction the inherited unique counter reports). This is the
    // population every per-tool metric uses.
    let read_count_after_filters = inputs.bam_stat.map(|b| b.unique).unwrap_or(0);

    SampleEnvelope {
        extractor_version: inputs.extractor_version.to_string(),
        schema_version: SCHEMA_VERSION.to_string(),
        git_commit: inputs.git_commit.to_string(),
        sample_id: inputs.sample_id.to_string(),
        bam_path: inputs.bam_path.to_string(),
        bam_md5: inputs.bam_md5.to_string(),
        bam_size_bytes: inputs.bam_size_bytes,
        gtf_path: inputs.gtf_path.to_string(),
        gtf_md5: inputs.gtf_md5.to_string(),
        reference_fasta_md5: inputs.reference_fasta_md5.map(|s| s.to_string()),
        library_prep: inputs.library_prep.to_string(),
        paired_end: inputs.paired_end,
        strandedness: inputs.strandedness.to_string(),
        read_length_max,
        read_length_mean,
        read_count_total,
        read_count_after_filters,
        filters: inputs.filters,
        qc_flags: inputs.qc_flags,
        runtime_seconds: inputs.runtime_seconds,
        peak_rss_mb: inputs.peak_rss_mb,
        bam_stat: bam_stat_block,
        infer_experiment: inputs
            .infer_experiment
            .map(InferExperimentBlock::from_result),
        read_distribution: inputs
            .read_distribution
            .map(ReadDistributionBlock::from_result),
        read_duplication: inputs
            .read_duplication
            .map(ReadDuplicationBlock::from_result),
        junction_annotation: inputs
            .junction_annotation
            .map(JunctionAnnotationBlock::from_result),
        junction_saturation: inputs
            .junction_saturation
            .map(JunctionSaturationBlock::from_result),
        inner_distance: inputs.inner_distance.map(InnerDistanceBlock::from_result),
        tin: inputs.tin.map(TinBlock::from_result),
        preseq: inputs.preseq.map(PreseqBlock::from_result),
        qualimap: inputs.qualimap.map(QualimapBlock::from_result),
        dupradar: inputs.count_result.map(|c| {
            DupradarBlock::from_parts(
                c.gene_counts.len() as u64,
                inputs.dupradar_genes_with_reads,
                inputs.dupradar_genes_with_duplication,
                inputs.dupradar_fit,
            )
        }),
        featurecounts: inputs.count_result.map(FeatureCountsBlock::from_result),
    }
}

/// Serialize an envelope to a pretty-printed JSON file.
pub fn write(envelope: &SampleEnvelope, path: &Path) -> Result<()> {
    let f = std::fs::File::create(path)
        .with_context(|| format!("creating envelope file: {}", path.display()))?;
    let mut writer = std::io::BufWriter::new(f);
    serde_json::to_writer_pretty(&mut writer, envelope)
        .with_context(|| format!("serializing envelope: {}", path.display()))?;
    use std::io::Write;
    writer.write_all(b"\n")?;
    writer.flush()?;
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Catch silent SCHEMA_VERSION drift: the constant must appear verbatim
    /// somewhere in the schema JSON. Both are bumped together at every wire
    /// format change.
    #[test]
    fn schema_version_appears_in_schema_file() {
        let schema_text = include_str!("../schema/v1/liquidqc.schema.json");
        assert!(
            schema_text.contains(SCHEMA_VERSION),
            "SCHEMA_VERSION constant `{SCHEMA_VERSION}` is not present in schema file"
        );
    }

    #[test]
    fn empty_envelope_has_all_required_top_level_keys() {
        // Build a minimal envelope with no metric blocks; all schema-required
        // fields must still be present in the JSON output.
        let env = SampleEnvelope {
            extractor_version: "0.0.0".into(),
            schema_version: SCHEMA_VERSION.into(),
            git_commit: "deadbee".into(),
            sample_id: "test".into(),
            bam_path: "/tmp/test.bam".into(),
            bam_md5: "00000000000000000000000000000000".into(),
            bam_size_bytes: 0,
            gtf_path: "/tmp/test.gtf".into(),
            gtf_md5: "00000000000000000000000000000000".into(),
            reference_fasta_md5: None,
            library_prep: "unknown".into(),
            paired_end: false,
            strandedness: "unknown".into(),
            read_length_max: 0,
            read_length_mean: 0.0,
            read_count_total: 0,
            read_count_after_filters: 0,
            filters: Filters::from_mapq(30),
            qc_flags: vec![],
            runtime_seconds: 0.0,
            peak_rss_mb: 0.0,
            bam_stat: None,
            infer_experiment: None,
            read_distribution: None,
            read_duplication: None,
            junction_annotation: None,
            junction_saturation: None,
            inner_distance: None,
            tin: None,
            preseq: None,
            qualimap: None,
            dupradar: None,
            featurecounts: None,
        };
        let v: serde_json::Value = serde_json::to_value(&env).unwrap();
        let schema_text = include_str!("../schema/v1/liquidqc.schema.json");
        let schema: serde_json::Value = serde_json::from_str(schema_text).unwrap();
        let required = schema["required"].as_array().expect("schema required[]");
        for r in required {
            let key = r.as_str().unwrap();
            assert!(
                v.get(key).is_some(),
                "envelope missing schema-required field `{key}`"
            );
        }
    }

    #[test]
    fn filters_from_mapq_defaults_match_phase1_contract() {
        let f = Filters::from_mapq(30);
        assert_eq!(f.min_mapq, 30);
        assert!(!f.require_proper_pair);
        assert!(f.drop_secondary);
        assert!(f.drop_supplementary);
        assert!(!f.drop_duplicates);
    }
}
