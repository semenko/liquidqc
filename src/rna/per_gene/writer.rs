//! Per-gene Parquet writer (Apache Arrow + Parquet).
//!
//! Owns the `PerGeneRow` row schema, the `PerGeneFileMeta` KV metadata,
//! conversion of per-worker merged `PerGeneAccumulator` state into rows
//! (`build_rows`), and the Arrow schema + streaming `ArrowWriter`
//! invocation (`write_per_gene_parquet`).
//!
//! Sample-id is stored ONLY in the Parquet file metadata, not as a
//! column, per the implementation brief — so cross-sample concatenation
//! is a pure stack with sample_id pulled from file metadata.

use anyhow::Result;
use indexmap::IndexMap;
use std::path::Path;

use crate::gtf::Gene;
use crate::rna::dupradar::counting::GeneCounts;

use super::state::{PerGeneState, KMER4_CARD};

/// One row of the per-gene Parquet table. Populated by
/// `PerGeneAccumulator::finalize`.
#[derive(Debug, Default, Clone)]
pub struct PerGeneRow {
    pub gene_id: String,
    pub gene_symbol: Option<String>,
    pub chrom: String,
    pub gene_start: u64,
    pub gene_end: u64,
    pub strand: String,
    pub effective_length: u64,
    pub biotype: Option<String>,

    pub fragment_count_assigned: u32,
    pub primary_reads: u32,
    pub primary_reads_unique: u32,
    pub primary_reads_unique_nodup: u32,
    pub primary_reads_multi: u32,
    pub fc_reads: u32,

    pub tlen_bin_lt_50: u32,
    pub tlen_bin_50_80: u32,
    pub tlen_bin_80_120: u32,
    pub tlen_bin_120_160: u32,
    pub tlen_bin_160_200: u32,
    pub tlen_bin_200_300: u32,
    pub tlen_bin_300_500: u32,
    pub tlen_bin_gt_500: u32,
    pub tlen_overflow: u32,
    pub tlen_count: u32,
    pub tlen_mean: Option<f64>,
    pub tlen_median_approx: Option<f64>,
    pub tlen_min: Option<u32>,
    pub tlen_max: Option<u32>,

    pub em_pair_count_used: u32,
    pub em_pair_count_skipped_non_acgt: u32,
    pub em_pair_count_skipped_oob: u32,
    /// Sparse 4-mer histogram. Only nonzero entries.
    pub em_5p: Vec<(String, u32)>,
    pub em_3p: Vec<(String, u32)>,
    pub em_shannon_5p: Option<f64>,
    pub em_shannon_3p: Option<f64>,
    pub em_jsd_5p_3p: Option<f64>,

    pub decile_coverage: [u32; super::state::DECILE_COUNT],

    pub soft_clip_rate_5p: f64,
    pub soft_clip_rate_3p: f64,
    pub primary_reads_with_5p_clip: u32,
    pub primary_reads_with_3p_clip: u32,

    pub exon_aligned_bp: u64,
    pub intron_aligned_bp: u64,
    /// Phase 4: aligned bp inside any merged CDS bbox.
    pub cds_aligned_bp: u64,
    /// Phase 4: `exon_aligned_bp - cds_aligned_bp` (exonic-but-not-CDS).
    pub utr_aligned_bp: u64,
    /// `cds_aligned_bp / exon_aligned_bp`. None when no exonic bp.
    pub cds_utr_ratio: Option<f64>,
    /// `intron_aligned_bp / (intron_aligned_bp + exon_aligned_bp)`.
    pub intron_retention_rate: Option<f64>,
}

/// Key/value metadata to attach to the Parquet file.
#[derive(Debug, Clone, Default)]
pub struct PerGeneFileMeta {
    pub sample_id: String,
    pub extractor_version: String,
    pub git_commit: String,
    pub per_gene_schema_version: String,
    pub bam_md5: String,
    pub gtf_md5: String,
    pub library_prep: String,
    pub min_gene_reads_threshold: u32,
    pub n_genes_total: u64,
    pub n_genes_emitted: u64,
    pub paired_end: bool,
    pub strandedness: String,
}

/// Stable list of column names in declaration order. Used by both the
/// writer and the envelope `per_gene.columns` field so consumers can
/// validate without opening the Parquet file.
pub fn column_names() -> Vec<String> {
    vec![
        "gene_id",
        "gene_symbol",
        "chrom",
        "gene_start",
        "gene_end",
        "strand",
        "effective_length",
        "biotype",
        "fragment_count_assigned",
        "primary_reads",
        "primary_reads_unique",
        "primary_reads_unique_nodup",
        "primary_reads_multi",
        "fc_reads",
        "tlen_bin_lt_50",
        "tlen_bin_50_80",
        "tlen_bin_80_120",
        "tlen_bin_120_160",
        "tlen_bin_160_200",
        "tlen_bin_200_300",
        "tlen_bin_300_500",
        "tlen_bin_gt_500",
        "tlen_overflow",
        "tlen_count",
        "tlen_mean",
        "tlen_median_approx",
        "tlen_min",
        "tlen_max",
        "em_pair_count_used",
        "em_pair_count_skipped_non_acgt",
        "em_pair_count_skipped_oob",
        "em_5p",
        "em_3p",
        "em_shannon_5p",
        "em_shannon_3p",
        "em_jsd_5p_3p",
        "decile_coverage",
        "soft_clip_rate_5p",
        "soft_clip_rate_3p",
        "primary_reads_with_5p_clip",
        "primary_reads_with_3p_clip",
        "exon_aligned_bp",
        "intron_aligned_bp",
        "cds_aligned_bp",
        "utr_aligned_bp",
        "cds_utr_ratio",
        "intron_retention_rate",
    ]
    .into_iter()
    .map(String::from)
    .collect()
}

/// Output of `PerGeneAccumulator::finalize` — Parquet rows after
/// `--min-gene-reads` filtering, plus headline counts for the envelope.
#[derive(Debug)]
pub struct PerGeneOutput {
    pub rows: Vec<PerGeneRow>,
    pub n_genes_total: u64,
    pub n_genes_emitted: u64,
    pub min_gene_reads_threshold: u32,
}

/// Build the per-row table from per-worker per-gene state plus the
/// per-sample featureCounts gene_counts (for `primary_reads_unique` and
/// friends, which are already accumulated by the inherited counting
/// pipeline). Genes with `state.primary_reads < min_gene_reads` are
/// dropped from the table.
///
/// `genes` and `gene_counts` are indexed by `gene_id` and iterate in the
/// `GeneIdInterner` insertion order — i.e., the same order as
/// `states`. This invariant is enforced by debug_assert.
pub fn build_rows(
    states: &[PerGeneState],
    min_gene_reads: u32,
    genes: &IndexMap<String, Gene>,
    gene_counts: &IndexMap<String, GeneCounts>,
) -> PerGeneOutput {
    debug_assert_eq!(
        states.len(),
        genes.len(),
        "PerGeneAccumulator state length ({}) must match gene count ({})",
        states.len(),
        genes.len()
    );

    let n_genes_total = genes.len() as u64;
    let mut rows = Vec::new();

    for ((gene_id, gene), state) in genes.iter().zip(states.iter()) {
        if state.primary_reads < min_gene_reads {
            continue;
        }
        let counts = gene_counts.get(gene_id).cloned().unwrap_or_default();
        rows.push(build_row(gene_id, gene, &counts, state));
    }

    let n_genes_emitted = rows.len() as u64;
    PerGeneOutput {
        rows,
        n_genes_total,
        n_genes_emitted,
        min_gene_reads_threshold: min_gene_reads,
    }
}

fn build_row(gene_id: &str, gene: &Gene, counts: &GeneCounts, state: &PerGeneState) -> PerGeneRow {
    use crate::rna::fragmentomics::common::{
        jsd_log2_from_probs, normalize_kmer_counts, shannon_entropy_log2_from_probs,
    };
    // Normalize the per-gene 4-mer histograms once and reuse for entropy + JSD.
    let em_5p_view = state.em_5p_view();
    let em_3p_view = state.em_3p_view();
    let probs_5p = normalize_kmer_counts(em_5p_view);
    let probs_3p = normalize_kmer_counts(em_3p_view);
    let em_shannon_5p = probs_5p.and_then(shannon_entropy_log2_from_probs);
    let em_shannon_3p = probs_3p.and_then(shannon_entropy_log2_from_probs);
    let em_jsd_5p_3p = match (probs_5p, probs_3p) {
        (Some(pa), Some(pb)) => Some(jsd_log2_from_probs(&pa, &pb)),
        _ => None,
    };
    PerGeneRow {
        gene_id: gene_id.to_string(),
        gene_symbol: gene.attributes.get("gene_name").cloned(),
        chrom: gene.chrom.clone(),
        gene_start: gene.start,
        gene_end: gene.end,
        strand: gene.strand.to_string(),
        effective_length: gene.effective_length,
        biotype: gene
            .attributes
            .get("gene_biotype")
            .or_else(|| gene.attributes.get("gene_type"))
            .cloned(),

        fragment_count_assigned: state.fragment_count_assigned,
        primary_reads: state.primary_reads,
        primary_reads_unique: counts.all_unique.min(u32::MAX as u64) as u32,
        primary_reads_unique_nodup: counts.nodup_unique.min(u32::MAX as u64) as u32,
        primary_reads_multi: (counts.all_multi.saturating_sub(counts.all_unique))
            .min(u32::MAX as u64) as u32,
        fc_reads: counts.fc_reads.min(u32::MAX as u64) as u32,

        tlen_bin_lt_50: state.tlen_bins[0],
        tlen_bin_50_80: state.tlen_bins[1],
        tlen_bin_80_120: state.tlen_bins[2],
        tlen_bin_120_160: state.tlen_bins[3],
        tlen_bin_160_200: state.tlen_bins[4],
        tlen_bin_200_300: state.tlen_bins[5],
        tlen_bin_300_500: state.tlen_bins[6],
        tlen_bin_gt_500: state.tlen_bins[7],
        tlen_overflow: state.tlen_overflow,
        tlen_count: state.tlen_count,
        tlen_mean: state.tlen_mean(),
        tlen_median_approx: state.tlen_median_approx(),
        tlen_min: if state.tlen_count == 0 {
            None
        } else {
            Some(state.tlen_min)
        },
        tlen_max: if state.tlen_count == 0 {
            None
        } else {
            Some(state.tlen_max)
        },

        em_pair_count_used: state.em_pair_count_used,
        em_pair_count_skipped_non_acgt: state.em_pair_count_skipped_non_acgt,
        em_pair_count_skipped_oob: state.em_pair_count_skipped_oob,
        em_5p: kmer_array_to_pairs(em_5p_view),
        em_3p: kmer_array_to_pairs(em_3p_view),
        em_shannon_5p,
        em_shannon_3p,
        em_jsd_5p_3p,

        decile_coverage: state.decile_coverage,

        soft_clip_rate_5p: state.soft_clip_rate_5p().unwrap_or(0.0),
        soft_clip_rate_3p: state.soft_clip_rate_3p().unwrap_or(0.0),
        primary_reads_with_5p_clip: state.primary_reads_with_5p_clip,
        primary_reads_with_3p_clip: state.primary_reads_with_3p_clip,

        exon_aligned_bp: state.exon_aligned_bp,
        intron_aligned_bp: state.intron_aligned_bp,
        cds_aligned_bp: state.cds_aligned_bp,
        utr_aligned_bp: state.utr_aligned_bp(),
        cds_utr_ratio: state.cds_utr_ratio(),
        intron_retention_rate: state.intron_retention_rate(),
    }
}

/// Convert a 256-entry u32 4-mer count array into a sparse list of
/// `(kmer_string, count)` pairs (only nonzero entries). Order is
/// 4-mer-index ascending (AAAA → TTTT).
fn kmer_array_to_pairs(arr: &[u32; KMER4_CARD]) -> Vec<(String, u32)> {
    crate::rna::fragmentomics::common::kmer_array_iter_nonzero(arr).collect()
}

/// Build the Arrow schema for the per-gene Parquet table. Field order
/// matches `column_names()` exactly, and is part of the per-gene
/// `1.0.0-stub` schema contract.
pub fn build_arrow_schema() -> arrow::datatypes::Schema {
    use arrow::datatypes::{DataType, Field};
    use std::sync::Arc;

    // Map<Utf8, UInt32> for em_5p / em_3p sparse 4-mer histograms.
    // The Map struct field is a non-nullable `Struct<keys: Utf8, values: UInt32>`.
    // Note: arrow's `MapBuilder` defaults to nullable values; matching that
    // here keeps the schema and the builder output consistent.
    let map_struct_field = |inner_name: &str| -> Field {
        Field::new(
            inner_name,
            DataType::Struct(
                vec![
                    Field::new("keys", DataType::Utf8, false),
                    Field::new("values", DataType::UInt32, true),
                ]
                .into(),
            ),
            false,
        )
    };
    let map_field = |outer_name: &str| -> Field {
        Field::new(
            outer_name,
            DataType::Map(Arc::new(map_struct_field("entries")), false),
            true,
        )
    };

    // List<UInt32> for decile_coverage (length 10, but emitted as variable-length
    // for ergonomics; consumers verify length 10 from schema metadata).
    // arrow's `ListBuilder` defaults to nullable inner items; matching that.
    let decile_field = Field::new(
        "decile_coverage",
        DataType::List(Arc::new(Field::new("item", DataType::UInt32, true))),
        false,
    );

    arrow::datatypes::Schema::new(vec![
        Field::new("gene_id", DataType::Utf8, false),
        Field::new("gene_symbol", DataType::Utf8, true),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("gene_start", DataType::UInt64, false),
        Field::new("gene_end", DataType::UInt64, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("effective_length", DataType::UInt64, false),
        Field::new("biotype", DataType::Utf8, true),
        Field::new("fragment_count_assigned", DataType::UInt32, false),
        Field::new("primary_reads", DataType::UInt32, false),
        Field::new("primary_reads_unique", DataType::UInt32, false),
        Field::new("primary_reads_unique_nodup", DataType::UInt32, false),
        Field::new("primary_reads_multi", DataType::UInt32, false),
        Field::new("fc_reads", DataType::UInt32, false),
        Field::new("tlen_bin_lt_50", DataType::UInt32, false),
        Field::new("tlen_bin_50_80", DataType::UInt32, false),
        Field::new("tlen_bin_80_120", DataType::UInt32, false),
        Field::new("tlen_bin_120_160", DataType::UInt32, false),
        Field::new("tlen_bin_160_200", DataType::UInt32, false),
        Field::new("tlen_bin_200_300", DataType::UInt32, false),
        Field::new("tlen_bin_300_500", DataType::UInt32, false),
        Field::new("tlen_bin_gt_500", DataType::UInt32, false),
        Field::new("tlen_overflow", DataType::UInt32, false),
        Field::new("tlen_count", DataType::UInt32, false),
        Field::new("tlen_mean", DataType::Float64, true),
        Field::new("tlen_median_approx", DataType::Float64, true),
        Field::new("tlen_min", DataType::UInt32, true),
        Field::new("tlen_max", DataType::UInt32, true),
        Field::new("em_pair_count_used", DataType::UInt32, false),
        Field::new("em_pair_count_skipped_non_acgt", DataType::UInt32, false),
        Field::new("em_pair_count_skipped_oob", DataType::UInt32, false),
        map_field("em_5p"),
        map_field("em_3p"),
        Field::new("em_shannon_5p", DataType::Float64, true),
        Field::new("em_shannon_3p", DataType::Float64, true),
        Field::new("em_jsd_5p_3p", DataType::Float64, true),
        decile_field,
        Field::new("soft_clip_rate_5p", DataType::Float64, false),
        Field::new("soft_clip_rate_3p", DataType::Float64, false),
        Field::new("primary_reads_with_5p_clip", DataType::UInt32, false),
        Field::new("primary_reads_with_3p_clip", DataType::UInt32, false),
        Field::new("exon_aligned_bp", DataType::UInt64, false),
        Field::new("intron_aligned_bp", DataType::UInt64, false),
        Field::new("cds_aligned_bp", DataType::UInt64, false),
        Field::new("utr_aligned_bp", DataType::UInt64, false),
        Field::new("cds_utr_ratio", DataType::Float64, true),
        Field::new("intron_retention_rate", DataType::Float64, true),
    ])
}

/// Convert a slice of `PerGeneRow` into a single Arrow `RecordBatch`
/// matching `build_arrow_schema()`.
fn rows_to_record_batch(
    rows: &[PerGeneRow],
    schema: std::sync::Arc<arrow::datatypes::Schema>,
) -> Result<arrow::record_batch::RecordBatch> {
    use arrow::array::{
        ArrayRef, Float64Array, ListBuilder, MapBuilder, StringArray, StringBuilder, UInt32Array,
        UInt32Builder, UInt64Array,
    };
    use std::sync::Arc;

    let n = rows.len();

    let gene_id: ArrayRef = Arc::new(StringArray::from_iter_values(
        rows.iter().map(|r| r.gene_id.as_str()),
    ));
    let gene_symbol: ArrayRef = Arc::new(StringArray::from_iter(
        rows.iter().map(|r| r.gene_symbol.as_deref()),
    ));
    let chrom: ArrayRef = Arc::new(StringArray::from_iter_values(
        rows.iter().map(|r| r.chrom.as_str()),
    ));
    let gene_start: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.gene_start),
    ));
    let gene_end: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.gene_end),
    ));
    let strand: ArrayRef = Arc::new(StringArray::from_iter_values(
        rows.iter().map(|r| r.strand.as_str()),
    ));
    let effective_length: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.effective_length),
    ));
    let biotype: ArrayRef = Arc::new(StringArray::from_iter(
        rows.iter().map(|r| r.biotype.as_deref()),
    ));

    let fragment_count_assigned: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.fragment_count_assigned),
    ));
    let primary_reads: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.primary_reads),
    ));
    let primary_reads_unique: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.primary_reads_unique),
    ));
    let primary_reads_unique_nodup: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.primary_reads_unique_nodup),
    ));
    let primary_reads_multi: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.primary_reads_multi),
    ));
    let fc_reads: ArrayRef = Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.fc_reads),
    ));

    let mk_u32 = |f: fn(&PerGeneRow) -> u32| -> ArrayRef {
        Arc::new(UInt32Array::from_iter_values(rows.iter().map(f)))
    };
    let tlen_bin_lt_50 = mk_u32(|r| r.tlen_bin_lt_50);
    let tlen_bin_50_80 = mk_u32(|r| r.tlen_bin_50_80);
    let tlen_bin_80_120 = mk_u32(|r| r.tlen_bin_80_120);
    let tlen_bin_120_160 = mk_u32(|r| r.tlen_bin_120_160);
    let tlen_bin_160_200 = mk_u32(|r| r.tlen_bin_160_200);
    let tlen_bin_200_300 = mk_u32(|r| r.tlen_bin_200_300);
    let tlen_bin_300_500 = mk_u32(|r| r.tlen_bin_300_500);
    let tlen_bin_gt_500 = mk_u32(|r| r.tlen_bin_gt_500);
    let tlen_overflow = mk_u32(|r| r.tlen_overflow);
    let tlen_count = mk_u32(|r| r.tlen_count);

    let tlen_mean: ArrayRef = Arc::new(Float64Array::from_iter(rows.iter().map(|r| r.tlen_mean)));
    let tlen_median_approx: ArrayRef = Arc::new(Float64Array::from_iter(
        rows.iter().map(|r| r.tlen_median_approx),
    ));
    let tlen_min: ArrayRef = Arc::new(UInt32Array::from_iter(rows.iter().map(|r| r.tlen_min)));
    let tlen_max: ArrayRef = Arc::new(UInt32Array::from_iter(rows.iter().map(|r| r.tlen_max)));

    let em_pair_count_used = mk_u32(|r| r.em_pair_count_used);
    let em_pair_count_skipped_non_acgt = mk_u32(|r| r.em_pair_count_skipped_non_acgt);
    let em_pair_count_skipped_oob = mk_u32(|r| r.em_pair_count_skipped_oob);

    // em_5p / em_3p as Map<Utf8, UInt32>.
    let build_map = |get: fn(&PerGeneRow) -> &Vec<(String, u32)>| -> Result<ArrayRef> {
        let mut b: MapBuilder<StringBuilder, UInt32Builder> =
            MapBuilder::new(None, StringBuilder::new(), UInt32Builder::new());
        for r in rows.iter() {
            let entries = get(r);
            for (k, v) in entries {
                b.keys().append_value(k);
                b.values().append_value(*v);
            }
            b.append(true)?; // entry present (sparse, may be empty)
        }
        Ok(Arc::new(b.finish()) as ArrayRef)
    };
    let em_5p_arr = build_map(|r| &r.em_5p)?;
    let em_3p_arr = build_map(|r| &r.em_3p)?;

    let em_shannon_5p: ArrayRef = Arc::new(Float64Array::from_iter(
        rows.iter().map(|r| r.em_shannon_5p),
    ));
    let em_shannon_3p: ArrayRef = Arc::new(Float64Array::from_iter(
        rows.iter().map(|r| r.em_shannon_3p),
    ));
    let em_jsd_5p_3p: ArrayRef =
        Arc::new(Float64Array::from_iter(rows.iter().map(|r| r.em_jsd_5p_3p)));

    // decile_coverage as List<UInt32>, exactly 10 entries per row.
    let decile_coverage_arr: ArrayRef = {
        let mut b: ListBuilder<UInt32Builder> = ListBuilder::new(UInt32Builder::new());
        for r in rows.iter() {
            for v in r.decile_coverage.iter() {
                b.values().append_value(*v);
            }
            b.append(true);
        }
        Arc::new(b.finish()) as ArrayRef
    };

    let soft_clip_rate_5p: ArrayRef = Arc::new(Float64Array::from_iter_values(
        rows.iter().map(|r| r.soft_clip_rate_5p),
    ));
    let soft_clip_rate_3p: ArrayRef = Arc::new(Float64Array::from_iter_values(
        rows.iter().map(|r| r.soft_clip_rate_3p),
    ));
    let primary_reads_with_5p_clip = mk_u32(|r| r.primary_reads_with_5p_clip);
    let primary_reads_with_3p_clip = mk_u32(|r| r.primary_reads_with_3p_clip);

    let exon_aligned_bp: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.exon_aligned_bp),
    ));
    let intron_aligned_bp: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.intron_aligned_bp),
    ));
    let cds_aligned_bp: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.cds_aligned_bp),
    ));
    let utr_aligned_bp: ArrayRef = Arc::new(UInt64Array::from_iter_values(
        rows.iter().map(|r| r.utr_aligned_bp),
    ));
    let cds_utr_ratio: ArrayRef = Arc::new(Float64Array::from_iter(
        rows.iter().map(|r| r.cds_utr_ratio),
    ));
    let intron_retention_rate: ArrayRef = Arc::new(Float64Array::from_iter(
        rows.iter().map(|r| r.intron_retention_rate),
    ));

    let columns = vec![
        gene_id,
        gene_symbol,
        chrom,
        gene_start,
        gene_end,
        strand,
        effective_length,
        biotype,
        fragment_count_assigned,
        primary_reads,
        primary_reads_unique,
        primary_reads_unique_nodup,
        primary_reads_multi,
        fc_reads,
        tlen_bin_lt_50,
        tlen_bin_50_80,
        tlen_bin_80_120,
        tlen_bin_120_160,
        tlen_bin_160_200,
        tlen_bin_200_300,
        tlen_bin_300_500,
        tlen_bin_gt_500,
        tlen_overflow,
        tlen_count,
        tlen_mean,
        tlen_median_approx,
        tlen_min,
        tlen_max,
        em_pair_count_used,
        em_pair_count_skipped_non_acgt,
        em_pair_count_skipped_oob,
        em_5p_arr,
        em_3p_arr,
        em_shannon_5p,
        em_shannon_3p,
        em_jsd_5p_3p,
        decile_coverage_arr,
        soft_clip_rate_5p,
        soft_clip_rate_3p,
        primary_reads_with_5p_clip,
        primary_reads_with_3p_clip,
        exon_aligned_bp,
        intron_aligned_bp,
        cds_aligned_bp,
        utr_aligned_bp,
        cds_utr_ratio,
        intron_retention_rate,
    ];

    debug_assert_eq!(columns.len(), schema.fields().len());
    debug_assert_eq!(n, columns.first().map_or(0, |c| c.len()));
    Ok(arrow::record_batch::RecordBatch::try_new(schema, columns)?)
}

/// Per-gene table schema version embedded in Parquet KV metadata.
/// Bumped to 1.1.0-stub in Phase 4 with the addition of `cds_aligned_bp`,
/// `utr_aligned_bp`, `cds_utr_ratio`, and the rename of `intron_exon_ratio`
/// to `intron_retention_rate`.
pub const PER_GENE_SCHEMA_VERSION: &str = "1.1.0-stub";

/// Write the per-gene Parquet file with `SNAPPY` compression and
/// row-group size 10_000. File-level KV metadata embeds sample-id and
/// provenance fields. Returns the md5 of the written file.
pub fn write_per_gene_parquet(
    rows: &[PerGeneRow],
    meta: &PerGeneFileMeta,
    path: &Path,
) -> Result<String> {
    use anyhow::Context;
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::metadata::KeyValue;
    use parquet::file::properties::WriterProperties;
    use std::fs::File;
    use std::sync::Arc;

    let schema = Arc::new(build_arrow_schema());
    let batch = rows_to_record_batch(rows, schema.clone())?;

    let kv = vec![
        KeyValue::new(
            "liquidqc.sample_id".to_string(),
            Some(meta.sample_id.clone()),
        ),
        KeyValue::new(
            "liquidqc.extractor_version".to_string(),
            Some(meta.extractor_version.clone()),
        ),
        KeyValue::new(
            "liquidqc.git_commit".to_string(),
            Some(meta.git_commit.clone()),
        ),
        KeyValue::new(
            "liquidqc.per_gene_schema_version".to_string(),
            Some(meta.per_gene_schema_version.clone()),
        ),
        KeyValue::new("liquidqc.bam_md5".to_string(), Some(meta.bam_md5.clone())),
        KeyValue::new("liquidqc.gtf_md5".to_string(), Some(meta.gtf_md5.clone())),
        KeyValue::new(
            "liquidqc.library_prep".to_string(),
            Some(meta.library_prep.clone()),
        ),
        KeyValue::new(
            "liquidqc.min_gene_reads_threshold".to_string(),
            Some(meta.min_gene_reads_threshold.to_string()),
        ),
        KeyValue::new(
            "liquidqc.n_genes_total".to_string(),
            Some(meta.n_genes_total.to_string()),
        ),
        KeyValue::new(
            "liquidqc.n_genes_emitted".to_string(),
            Some(meta.n_genes_emitted.to_string()),
        ),
        KeyValue::new(
            "liquidqc.paired_end".to_string(),
            Some(meta.paired_end.to_string()),
        ),
        KeyValue::new(
            "liquidqc.strandedness".to_string(),
            Some(meta.strandedness.clone()),
        ),
    ];

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_row_count(Some(10_000))
        .set_key_value_metadata(Some(kv))
        .build();

    let file = File::create(path)
        .with_context(|| format!("Failed to create per-gene Parquet file: {}", path.display()))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .context("Failed to construct Parquet ArrowWriter")?;
    if !rows.is_empty() {
        writer
            .write(&batch)
            .context("Failed to write per-gene Parquet batch")?;
    }
    writer
        .close()
        .context("Failed to close per-gene Parquet writer")?;

    crate::hash::md5(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_pairs_drops_zero_keys() {
        let mut arr = [0u32; KMER4_CARD];
        arr[0] = 7; // AAAA
        arr[255] = 3; // TTTT
        let pairs = kmer_array_to_pairs(&arr);
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], ("AAAA".to_string(), 7));
        assert_eq!(pairs[1], ("TTTT".to_string(), 3));
    }

    #[test]
    fn arrow_schema_columns_match_column_names() {
        let schema = build_arrow_schema();
        let names = column_names();
        assert_eq!(schema.fields().len(), names.len());
        for (field, name) in schema.fields().iter().zip(names.iter()) {
            assert_eq!(field.name(), name);
        }
    }

    fn make_row(gene_id: &str, primary_reads: u32) -> PerGeneRow {
        let mut row = PerGeneRow {
            gene_id: gene_id.to_string(),
            chrom: "chr1".to_string(),
            strand: "+".to_string(),
            primary_reads,
            ..Default::default()
        };
        // Populate a couple of em_5p sparse entries to exercise the Map column.
        row.em_5p = vec![("AAAA".to_string(), 5), ("CCCC".to_string(), 2)];
        row.em_3p = vec![("GGGG".to_string(), 3)];
        row.decile_coverage = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        row.tlen_min = Some(50);
        row.tlen_max = Some(180);
        row.tlen_mean = Some(120.0);
        row
    }

    #[test]
    fn parquet_round_trip_three_rows() {
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
        use std::fs::File;
        use tempfile::tempdir;

        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("per_gene.parquet");

        let rows = vec![
            make_row("ENSG00000001", 10),
            make_row("ENSG00000002", 25),
            make_row("ENSG00000003", 100),
        ];
        let meta = PerGeneFileMeta {
            sample_id: "test".to_string(),
            extractor_version: "0.1.0".to_string(),
            git_commit: "abc1234".to_string(),
            per_gene_schema_version: PER_GENE_SCHEMA_VERSION.to_string(),
            bam_md5: "0".repeat(32),
            gtf_md5: "1".repeat(32),
            library_prep: "unknown".to_string(),
            min_gene_reads_threshold: 1,
            n_genes_total: 3,
            n_genes_emitted: 3,
            paired_end: false,
            strandedness: "unknown".to_string(),
        };
        let md5 = write_per_gene_parquet(&rows, &meta, &path).expect("write");
        assert_eq!(md5.len(), 32);

        // Read back and verify schema + row count + KV metadata.
        let f = File::open(&path).expect("open");
        let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("builder");

        let parquet_meta = builder.metadata();
        let kv = parquet_meta
            .file_metadata()
            .key_value_metadata()
            .expect("kv");
        let kv_map: std::collections::HashMap<&str, Option<&str>> = kv
            .iter()
            .map(|k| (k.key.as_str(), k.value.as_deref()))
            .collect();
        assert_eq!(
            kv_map.get("liquidqc.sample_id").copied().flatten(),
            Some("test")
        );
        assert_eq!(
            kv_map
                .get("liquidqc.per_gene_schema_version")
                .copied()
                .flatten(),
            Some(PER_GENE_SCHEMA_VERSION),
        );
        assert_eq!(
            kv_map.get("liquidqc.n_genes_emitted").copied().flatten(),
            Some("3")
        );

        let arrow_schema = builder.schema().clone();
        assert_eq!(arrow_schema.fields().len(), column_names().len());

        let reader = builder.build().expect("reader");
        let mut total_rows = 0;
        for batch_result in reader {
            let batch = batch_result.expect("batch");
            total_rows += batch.num_rows();
        }
        assert_eq!(total_rows, 3);
    }

    #[test]
    fn parquet_zero_rows_writes_valid_file() {
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
        use std::fs::File;
        use tempfile::tempdir;

        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("empty.parquet");
        let meta = PerGeneFileMeta {
            sample_id: "empty".to_string(),
            extractor_version: "0.1.0".to_string(),
            git_commit: "abc".to_string(),
            per_gene_schema_version: PER_GENE_SCHEMA_VERSION.to_string(),
            bam_md5: "0".repeat(32),
            gtf_md5: "0".repeat(32),
            library_prep: "unknown".to_string(),
            min_gene_reads_threshold: 1000000,
            n_genes_total: 100,
            n_genes_emitted: 0,
            paired_end: true,
            strandedness: "forward".to_string(),
        };
        write_per_gene_parquet(&[], &meta, &path).expect("write empty");
        let f = File::open(&path).expect("open");
        let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("builder");
        assert_eq!(builder.metadata().file_metadata().num_rows(), 0);
        assert_eq!(builder.schema().fields().len(), column_names().len());
    }
}
