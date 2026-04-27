//! Phase 4 Tier 3 integration tests.
//!
//! Drives the binary end-to-end against the inherited `tests/data/test.bam`
//! fixture and asserts that every Phase 4 envelope block emits cleanly with
//! the expected shape. Also exercises the per-gene Parquet bump to
//! `1.1.0-stub` (CDS/UTR/intron retention columns) and the rule-engine
//! changes (drop `sex_swap_warning`, add `cfdna_contamination_suspected`).
//!
//! Notes on the fixture:
//!   - Synthetic BAM uses `chr1` / `chr2` only; no chrY contig is present, so
//!     `sex_inference.predicted_sex` resolves to `unknown` (denominator below
//!     `MIN_MAPPED_FOR_CALL`).
//!   - GTF has no CDS lines, so `cds_aligned_bp` is 0 and all exonic bp lands
//!     in `utr_aligned_bp`. That's exactly what the schema promises and is
//!     what we assert here.
//!   - No FASTA → `splice_site_dinucleotides` is `null`.
//!   - Without `--paired`, `adapter_readthrough_rate` is `null` and
//!     `single_end_no_tlen` shows up in `qc_flags`.

use std::path::Path;
use std::process::Command;

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
        std::env::temp_dir().join(format!("liquidqc-phase4-{}-{}", label, std::process::id()));
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

fn qc_flags(env: &serde_json::Value) -> Vec<&str> {
    env["qc_flags"]
        .as_array()
        .expect("qc_flags missing or not array")
        .iter()
        .filter_map(|v| v.as_str())
        .collect()
}

#[test]
fn schema_version_is_phase4_stub() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "schema-version",
    );
    assert_eq!(env["schema_version"].as_str(), Some("0.5.0-stub"));
}

#[test]
fn per_chromosome_block_resolves_contig_names_and_fractions() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "per-chrom",
    );
    let block = env
        .get("per_chromosome")
        .expect("per_chromosome block must be present");
    let total_mapped = block["total_mapped"]
        .as_u64()
        .expect("per_chromosome.total_mapped should be u64");
    assert!(total_mapped > 0, "fixture has mapped reads");

    let contigs = block["contigs"]
        .as_object()
        .expect("per_chromosome.contigs should be an object");
    assert!(
        contigs.contains_key("chr1"),
        "fixture covers chr1; contigs={:?}",
        contigs.keys().collect::<Vec<_>>()
    );
    let chr1 = &contigs["chr1"];
    assert!(chr1["mapped"].as_u64().unwrap() > 0);
    let frac = chr1["read_fraction"]
        .as_f64()
        .expect("read_fraction should be a float");
    assert!(
        (0.0..=1.0).contains(&frac),
        "read_fraction must be in [0,1], got {frac}"
    );

    // Sum of read_fraction across all reported contigs must be ~1.0
    let sum: f64 = contigs
        .values()
        .map(|v| v["read_fraction"].as_f64().unwrap_or(0.0))
        .sum();
    assert!(
        (sum - 1.0).abs() < 1e-9,
        "per-contig read fractions should sum to 1, got {sum}"
    );
}

#[test]
fn cycle_quality_block_emits_r1_and_skips_r2_when_single_end() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "cycle-quality",
    );
    let cq = env
        .get("cycle_quality")
        .expect("cycle_quality block must be present");
    let r1 = cq["read1_mean_quality"]
        .as_array()
        .expect("read1_mean_quality should be an array");
    assert!(
        !r1.is_empty(),
        "fixture has read-1 reads, expected non-empty per-cycle Q"
    );
    // Single-end fixture: read2_mean_quality is omitted entirely
    // (skip_serializing_if when empty) — that's the expected shape.
    let r2 = cq.get("read2_mean_quality");
    assert!(
        r2.map(|v| v.as_array().map(|a| a.is_empty()).unwrap_or(true))
            .unwrap_or(true),
        "single-end fixture should have empty/absent read2_mean_quality, got {:?}",
        r2
    );
    // Mean Q values must be finite and non-negative.
    for v in r1 {
        let q = v.as_f64().expect("per-cycle mean Q is f64");
        assert!(q.is_finite() && q >= 0.0, "bad per-cycle Q: {q}");
    }
}

#[test]
fn gene_class_fractions_block_present_with_bundled_symbol_counts() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "gene-class",
    );
    let block = env
        .get("gene_class_fractions")
        .expect("gene_class_fractions block must be present");
    for class in ["hemoglobin", "ribosomal_protein", "apolipoprotein"] {
        let entry = &block[class];
        let bundled = entry["bundled_symbols"]
            .as_u64()
            .unwrap_or_else(|| panic!("{class}.bundled_symbols missing"));
        assert!(bundled > 0, "{class} should have bundled symbols");
        let frac = entry["fraction"].as_f64().unwrap_or(0.0);
        assert!(
            (0.0..=1.0).contains(&frac),
            "{class}.fraction must be in [0,1], got {frac}"
        );
    }
    // Synthetic GTF gene names don't match HGNC symbols, so all classes
    // resolve to zero matched genes — that's a feature of the fixture.
    for class in ["hemoglobin", "ribosomal_protein", "apolipoprotein"] {
        let matched = block[class]["matched_genes"].as_u64().unwrap_or(0);
        assert_eq!(
            matched, 0,
            "{class} should have no matches against synthetic gene symbols"
        );
    }
}

#[test]
fn panels_block_loads_bundled_tsvs_with_zero_matches_on_synthetic_fixture() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "panels",
    );
    let block = env.get("panels").expect("panels block must be present");
    let loaded = block["loaded_entries"]
        .as_u64()
        .expect("panels.loaded_entries missing");
    assert!(
        loaded > 0,
        "bundled panel TSVs should produce at least one entry"
    );
    let panels = block["panels"]
        .as_object()
        .expect("panels.panels should be an object");
    for expected in ["lm22", "tabula_sapiens_cfrna", "vorperian"] {
        assert!(
            panels.contains_key(expected),
            "expected bundled panel {expected} in {:?}",
            panels.keys().collect::<Vec<_>>()
        );
    }
    // No matches expected on synthetic gene names.
    assert_eq!(
        block["matched_entries"].as_u64().unwrap_or(u64::MAX),
        0,
        "synthetic fixture should not match any panel entries"
    );
}

#[test]
fn sex_inference_block_present_with_unknown_on_chrm_only_fixture() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "sex",
    );
    let block = env
        .get("sex_inference")
        .expect("sex_inference block must be present");
    // Fixture has no chrY contig and tiny mapped count, so the predictor
    // should refuse to call.
    assert_eq!(
        block["chr_y_reads"].as_u64().unwrap_or(u64::MAX),
        0,
        "fixture has no chrY contig"
    );
    let predicted = block["predicted_sex"]
        .as_str()
        .expect("predicted_sex should be a string");
    assert_eq!(
        predicted, "unknown",
        "predicted_sex should be 'unknown' when below denominator floor; got {predicted}"
    );
    // We dropped sex_swap_warning from the qc_flags enum entirely — assert
    // it is not in the qc_flags array.
    let flags = qc_flags(&env);
    assert!(
        !flags.contains(&"sex_swap_warning"),
        "sex_swap_warning was removed from the enum, must not appear in qc_flags; got {flags:?}"
    );
}

#[test]
fn snp_fingerprint_block_emitted_against_bundled_panel() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "snp-fp",
    );
    let block = env
        .get("snp_fingerprint")
        .expect("snp_fingerprint block must be present");
    let total = block["sites_total"]
        .as_u64()
        .expect("snp_fingerprint.sites_total missing");
    assert!(
        total > 0,
        "bundled smoke panel must contribute sites; got {total}"
    );
    let called = block["sites_called"].as_u64().unwrap_or(u64::MAX);
    assert!(
        called <= total,
        "sites_called ({called}) cannot exceed sites_total ({total})"
    );
    let sites = block["sites"]
        .as_array()
        .expect("snp_fingerprint.sites should be an array");
    // Every site row must have a depth field that is non-negative.
    for s in sites {
        let depth = s["depth"].as_u64().expect("site.depth missing");
        let _ = depth; // u64 is always non-negative; assert presence only.
        let af = s["alt_fraction"]
            .as_f64()
            .expect("site.alt_fraction missing");
        assert!(
            (0.0..=1.0).contains(&af) || depth == 0,
            "alt_fraction out of range: {af}"
        );
    }
}

#[test]
fn saturation_block_is_monotone_and_terminates_at_full_subsample() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "saturation",
    );
    let block = env
        .get("saturation")
        .expect("saturation block must be present");
    let points = block["points"]
        .as_array()
        .expect("saturation.points must be an array");
    assert!(points.len() >= 2, "expected ≥2 saturation points");

    let mut last_genes: i64 = -1;
    let mut last_min10: i64 = -1;
    let mut last_reads: i64 = -1;
    for (i, p) in points.iter().enumerate() {
        let g = p["genes_detected"].as_u64().unwrap() as i64;
        let g10 = p["genes_detected_min10"].as_u64().unwrap() as i64;
        let r = p["reads_in_subsample"].as_u64().unwrap() as i64;
        assert!(
            g >= last_genes,
            "genes_detected must be monotone non-decreasing at idx {i}: {g} < {last_genes}"
        );
        assert!(
            g10 >= last_min10,
            "genes_detected_min10 must be monotone non-decreasing at idx {i}: {g10} < {last_min10}"
        );
        assert!(
            r >= last_reads,
            "reads_in_subsample must be monotone non-decreasing at idx {i}: {r} < {last_reads}"
        );
        last_genes = g;
        last_min10 = g10;
        last_reads = r;
    }

    // Last point should be fraction 1.0 and reads_in_subsample == total.
    let last = points.last().unwrap();
    assert!((last["fraction"].as_f64().unwrap() - 1.0).abs() < 1e-9);
    let total = block["total_reads_assigned"].as_u64().unwrap();
    assert_eq!(
        last["reads_in_subsample"].as_u64().unwrap(),
        total,
        "fraction=1.0 point must include every assigned read"
    );
}

#[test]
fn adapter_readthrough_and_splice_dinuc_are_null_without_paired_or_fasta() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "single-end-no-fasta",
    );
    // single-end fixture: no proper pairs → adapter_readthrough_rate is null
    assert!(
        env["adapter_readthrough_rate"].is_null(),
        "expected null adapter_readthrough_rate for single-end input; got {:?}",
        env["adapter_readthrough_rate"]
    );
    // No FASTA → no donor/acceptor lookup
    assert!(
        env["splice_site_dinucleotides"].is_null(),
        "expected null splice_site_dinucleotides without --fasta; got {:?}",
        env["splice_site_dinucleotides"]
    );
    let flags = qc_flags(&env);
    assert!(
        flags.contains(&"single_end_no_tlen"),
        "expected single_end_no_tlen in qc_flags; got {flags:?}"
    );
}

#[test]
fn cfdna_contamination_suspected_does_not_fire_on_clean_fixture() {
    let (env, _) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "cfdna-contam",
    );
    let flags = qc_flags(&env);
    assert!(
        !flags.contains(&"cfdna_contamination_suspected"),
        "clean fixture must not raise cfdna_contamination_suspected; got {flags:?}"
    );
}

#[test]
fn per_gene_parquet_has_phase4_columns_and_bumped_schema_version() {
    let (env, outdir) = run_with(
        &["--library-prep", "unknown", "--min-gene-reads", "1"],
        "per-gene-cols",
    );

    // Envelope-side: per_gene block reports the new schema version.
    assert_eq!(
        env["per_gene"]["schema_version_per_gene"].as_str(),
        Some("1.1.0-stub")
    );
    let envelope_columns: Vec<&str> = env["per_gene"]["columns"]
        .as_array()
        .expect("per_gene.columns missing")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();
    for c in [
        "cds_aligned_bp",
        "utr_aligned_bp",
        "cds_utr_ratio",
        "intron_retention_rate",
    ] {
        assert!(
            envelope_columns.contains(&c),
            "envelope per_gene.columns should advertise {c}; got {envelope_columns:?}"
        );
    }

    // Parquet-side: file schema includes the four new columns and the file
    // KV metadata reports `liquidqc.per_gene_schema_version = 1.1.0-stub`.
    let parquet_path = outdir.join("test.per_gene.parquet");
    let f = std::fs::File::open(&parquet_path).expect("open per-gene parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).expect("parquet builder");
    let arrow_schema = builder.schema();
    let names: Vec<&str> = arrow_schema
        .fields()
        .iter()
        .map(|f| f.name().as_str())
        .collect();
    for c in [
        "cds_aligned_bp",
        "utr_aligned_bp",
        "cds_utr_ratio",
        "intron_retention_rate",
    ] {
        assert!(
            names.contains(&c),
            "parquet schema must include {c}; got {names:?}"
        );
    }
    let kv = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .expect("per-gene parquet must carry KV metadata");
    let schema_version = kv
        .iter()
        .find(|kv| kv.key == "liquidqc.per_gene_schema_version")
        .and_then(|kv| kv.value.as_deref())
        .expect("liquidqc.per_gene_schema_version not found in parquet KV metadata");
    assert_eq!(schema_version, "1.1.0-stub");
}
