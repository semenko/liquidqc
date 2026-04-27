//! Phase 3 Tier 2 per-gene Parquet integration tests.
//!
//! Exercises end-to-end emission of `<sample_id>.per_gene.parquet`
//! alongside the JSON envelope's `per_gene` block. Reads the Parquet
//! back via the official parquet-rs reader to validate schema, KV
//! metadata, row count, and field semantics.

use std::path::Path;
use std::process::Command;

use arrow::array::Array;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn binary() -> String {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_liquidqc") {
        return path;
    }
    let status = Command::new("cargo")
        .args(["build"])
        .status()
        .expect("Failed to run cargo build");
    assert!(status.success(), "cargo build failed");
    let path = Path::new("target/debug/liquidqc");
    assert!(path.exists(), "Binary not found at {:?}", path);
    path.to_str().unwrap().to_string()
}

fn run_with(args: &[&str], label: &str) -> (serde_json::Value, std::path::PathBuf) {
    let outdir =
        std::env::temp_dir().join(format!("liquidqc-phase3-{}-{}", label, std::process::id()));
    let _ = std::fs::remove_dir_all(&outdir);
    std::fs::create_dir_all(&outdir).unwrap();

    let mut full_args: Vec<&str> =
        vec!["rna", "tests/data/test.bam", "--gtf", "tests/data/test.gtf"];
    full_args.extend_from_slice(args);
    full_args.extend_from_slice(&["--outdir", outdir.to_str().unwrap(), "--skip-dup-check"]);

    let out = Command::new(binary())
        .args(&full_args)
        .output()
        .expect("Failed to execute liquidqc rna");
    assert!(
        out.status.success(),
        "liquidqc rna failed: args={:?} status={:?} stderr={}",
        full_args,
        out.status,
        String::from_utf8_lossy(&out.stderr),
    );

    let envelope_path = outdir.join("test.liquidqc.json");
    let text = std::fs::read_to_string(&envelope_path)
        .unwrap_or_else(|e| panic!("envelope missing at {envelope_path:?}: {e}"));
    (
        serde_json::from_str(&text).expect("envelope is not valid JSON"),
        outdir,
    )
}

#[test]
fn per_gene_parquet_sibling_exists_with_envelope_block() {
    let (env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "exists",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    assert!(
        parquet_path.exists(),
        "per-gene Parquet should exist at {parquet_path:?}"
    );
    let metadata = std::fs::metadata(&parquet_path).expect("stat parquet");
    assert!(metadata.len() > 0, "Parquet file should be non-empty");

    let block = env
        .get("per_gene")
        .expect("envelope must include per_gene block");
    assert_eq!(
        block["parquet_path"].as_str(),
        Some("test.per_gene.parquet"),
        "envelope per_gene.parquet_path should be the sibling filename"
    );
    assert_eq!(
        block["min_gene_reads_threshold"].as_u64(),
        Some(1),
        "min_gene_reads_threshold should reflect --min-gene-reads"
    );
    assert_eq!(
        block["schema_version_per_gene"].as_str(),
        Some("1.1.0-stub"),
        "schema_version_per_gene should be 1.1.0-stub"
    );
    let md5 = block["parquet_md5"].as_str().expect("parquet_md5 missing");
    assert_eq!(md5.len(), 32, "parquet_md5 should be 32-hex");
    assert!(
        md5.chars().all(|c| c.is_ascii_hexdigit()),
        "parquet_md5 should be lowercase hex; got {md5:?}"
    );
}

#[test]
fn per_gene_parquet_schema_columns_match_envelope() {
    let (env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "schema",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    let f = std::fs::File::open(&parquet_path).expect("open parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    let arrow_schema = builder.schema().clone();

    let envelope_columns: Vec<&str> = env["per_gene"]["columns"]
        .as_array()
        .expect("envelope per_gene.columns missing")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();

    assert_eq!(
        arrow_schema.fields().len(),
        envelope_columns.len(),
        "Arrow schema field count must match envelope columns"
    );
    for (field, expected) in arrow_schema.fields().iter().zip(envelope_columns.iter()) {
        assert_eq!(
            field.name(),
            expected,
            "column order in Parquet schema must match envelope per_gene.columns"
        );
    }

    // Spot-check a handful of column types — the full listing is documented
    // in the implementation plan; here we verify the load-bearing ones.
    use arrow::datatypes::DataType;
    let by_name: std::collections::HashMap<&str, &DataType> = arrow_schema
        .fields()
        .iter()
        .map(|f| (f.name().as_str(), f.data_type()))
        .collect();
    assert_eq!(by_name.get("gene_id"), Some(&&DataType::Utf8));
    assert_eq!(by_name.get("primary_reads"), Some(&&DataType::UInt32));
    assert_eq!(by_name.get("tlen_mean"), Some(&&DataType::Float64));
    matches!(
        by_name.get("decile_coverage"),
        Some(DataType::List(_) | DataType::FixedSizeList(_, _))
    );
    matches!(by_name.get("em_5p"), Some(DataType::Map(_, _)));
}

#[test]
fn per_gene_parquet_kv_metadata_is_consistent_with_envelope() {
    let (env, outdir) = run_with(
        &[
            "--library-prep",
            "neb_next_ultra_ii_directional",
            "--min-gene-reads",
            "1",
        ],
        "kv",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    let f = std::fs::File::open(&parquet_path).expect("open parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    let kv = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .expect("KV metadata missing");
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
        Some("1.1.0-stub")
    );
    assert_eq!(
        kv_map.get("liquidqc.library_prep").copied().flatten(),
        Some("neb_next_ultra_ii_directional")
    );

    let n_emitted_kv = kv_map
        .get("liquidqc.n_genes_emitted")
        .copied()
        .flatten()
        .expect("n_genes_emitted in KV");
    let n_emitted_env = env["per_gene"]["n_genes_emitted"].as_u64().unwrap();
    assert_eq!(
        n_emitted_kv,
        n_emitted_env.to_string(),
        "n_genes_emitted must agree between Parquet KV and envelope"
    );
}

#[test]
fn per_gene_parquet_row_count_matches_envelope() {
    let (env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "rowcount",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    let f = std::fs::File::open(&parquet_path).expect("open parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    let parquet_rows = builder.metadata().file_metadata().num_rows();
    let envelope_emitted = env["per_gene"]["n_genes_emitted"].as_u64().unwrap();
    assert_eq!(
        parquet_rows as u64, envelope_emitted,
        "Parquet num_rows ({parquet_rows}) must equal envelope n_genes_emitted ({envelope_emitted})"
    );
    assert!(
        parquet_rows > 0,
        "fixture should produce at least one gene row"
    );
}

#[test]
fn high_threshold_filters_to_zero_rows() {
    let (env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1000000"],
        "highthreshold",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    assert!(
        parquet_path.exists(),
        "Parquet must exist even when no rows survive"
    );
    let f = std::fs::File::open(&parquet_path).expect("open parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    assert_eq!(
        builder.metadata().file_metadata().num_rows(),
        0,
        "Parquet should have 0 rows above the threshold"
    );
    assert_eq!(env["per_gene"]["n_genes_emitted"].as_u64(), Some(0));
}

#[test]
fn single_end_fixture_produces_null_tlen_mean_per_gene() {
    let (_env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "single-end-tlen",
    );
    let parquet_path = outdir.join("test.per_gene.parquet");
    let f = std::fs::File::open(&parquet_path).expect("open parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    let mut reader = builder.build().expect("reader");
    let batch = reader
        .next()
        .expect("at least one batch")
        .expect("batch ok");
    let tlen_count_idx = batch
        .schema()
        .fields()
        .iter()
        .position(|f| f.name() == "tlen_count")
        .expect("tlen_count column");
    let tlen_mean_idx = batch
        .schema()
        .fields()
        .iter()
        .position(|f| f.name() == "tlen_mean")
        .expect("tlen_mean column");
    let tlen_count = batch
        .column(tlen_count_idx)
        .as_any()
        .downcast_ref::<arrow::array::UInt32Array>()
        .expect("UInt32 column");
    let tlen_mean = batch
        .column(tlen_mean_idx)
        .as_any()
        .downcast_ref::<arrow::array::Float64Array>()
        .expect("Float64 column");
    for i in 0..batch.num_rows() {
        // Single-end fixture: every row has tlen_count == 0 and tlen_mean == null.
        assert_eq!(tlen_count.value(i), 0, "row {i} tlen_count should be 0");
        assert!(tlen_mean.is_null(i), "row {i} tlen_mean should be null");
    }
}
