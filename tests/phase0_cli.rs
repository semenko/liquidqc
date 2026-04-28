//! Phase 0 CLI smoke tests for liquidqc.
//!
//! These verify the new subcommands (`dna`, `schema`, `version`) and the
//! `--library-prep` required-flag contract from the v1 plan. They do NOT
//! exercise BAM processing — the inherited integration_test.rs covers that.

use std::path::Path;
use std::process::Command;

/// Path to the liquidqc binary. Cargo sets `CARGO_BIN_EXE_liquidqc` for
/// integration tests; otherwise fall back to building the debug binary.
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

#[test]
fn version_subcommand_prints_version() {
    let out = Command::new(binary())
        .args(["version"])
        .output()
        .expect("Failed to execute liquidqc version");
    assert!(out.status.success(), "liquidqc version failed");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.starts_with("liquidqc "),
        "Expected stdout to start with `liquidqc `, got: {stdout:?}"
    );
    assert!(
        stdout.contains(env!("CARGO_PKG_VERSION")),
        "Expected stdout to contain crate version {}, got: {stdout:?}",
        env!("CARGO_PKG_VERSION"),
    );
}

#[test]
fn version_long_flag_prints_version() {
    // clap auto-handles --version too; sanity-check that it still works
    // alongside the explicit `version` subcommand.
    let out = Command::new(binary())
        .args(["--version"])
        .output()
        .expect("Failed to execute liquidqc --version");
    assert!(out.status.success(), "liquidqc --version failed");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains(env!("CARGO_PKG_VERSION")),
        "Expected --version to contain crate version, got: {stdout:?}"
    );
}

#[test]
fn schema_subcommand_emits_valid_json() {
    let out = Command::new(binary())
        .args(["schema"])
        .output()
        .expect("Failed to execute liquidqc schema");
    assert!(out.status.success(), "liquidqc schema failed");
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Must round-trip through serde_json (proves it's valid JSON, not just
    // a valid byte stream).
    let parsed: serde_json::Value =
        serde_json::from_str(&stdout).expect("liquidqc schema emitted invalid JSON");
    let title = parsed
        .get("title")
        .and_then(|v| v.as_str())
        .expect("schema missing `title`");
    assert!(
        title.contains("liquidqc"),
        "Expected schema title to mention liquidqc, got {title:?}"
    );
    assert_eq!(
        parsed.get("type").and_then(|v| v.as_str()),
        Some("object"),
        "Schema type should be `object`"
    );
    let required = parsed
        .get("required")
        .and_then(|v| v.as_array())
        .expect("schema missing `required` array");
    let required_names: Vec<&str> = required.iter().filter_map(|v| v.as_str()).collect();
    for must_have in [
        "extractor_version",
        "schema_version",
        "library_prep",
        "qc_flags",
    ] {
        assert!(
            required_names.contains(&must_have),
            "Schema `required` should contain `{must_have}`; got {required_names:?}"
        );
    }
}

#[test]
fn live_envelope_satisfies_schema_required_fields() {
    // Phase 1 replaces the static fixture check with a live-build assertion:
    // run the binary against tests/data, then verify the produced envelope
    // contains every field the schema marks as required. This eliminates the
    // drift risk a hand-edited fixture would have.
    let outdir = std::env::temp_dir().join(format!(
        "liquidqc-phase0-live-envelope-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&outdir);
    std::fs::create_dir_all(&outdir).unwrap();

    let out = Command::new(binary())
        .args([
            "rna",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            outdir.to_str().unwrap(),
            "--library-prep",
            "unknown",
            "--skip-dup-check",
        ])
        .output()
        .expect("Failed to execute liquidqc rna");
    assert!(
        out.status.success(),
        "liquidqc rna failed: status={:?}, stderr={}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );

    let envelope_path = outdir.join("test.liquidqc.json");
    let envelope_text = std::fs::read_to_string(&envelope_path)
        .unwrap_or_else(|e| panic!("envelope not found at {envelope_path:?}: {e}"));
    let envelope: serde_json::Value =
        serde_json::from_str(&envelope_text).expect("envelope is not valid JSON");
    let envelope_obj = envelope
        .as_object()
        .expect("envelope must be a JSON object");

    let schema_text =
        std::fs::read_to_string("schema/v1/liquidqc.schema.json").expect("Schema file missing");
    let schema: serde_json::Value =
        serde_json::from_str(&schema_text).expect("Schema is not valid JSON");
    let required: Vec<&str> = schema
        .get("required")
        .and_then(|v| v.as_array())
        .expect("schema missing `required`")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();

    for field in required {
        assert!(
            envelope_obj.contains_key(field),
            "envelope missing schema-required field `{field}` (file: {envelope_path:?})"
        );
    }
}

#[test]
fn dna_subcommand_exits_nonzero_with_clear_message() {
    let out = Command::new(binary())
        .args(["dna", "--bam", "ignored.bam"])
        .output()
        .expect("Failed to execute liquidqc dna");
    assert!(
        !out.status.success(),
        "liquidqc dna should exit non-zero; got status {:?}",
        out.status
    );
    assert_eq!(out.status.code(), Some(1), "Expected exit code 1");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("not implemented"),
        "Expected stderr to contain `not implemented`, got: {stderr:?}"
    );
}

#[test]
fn rna_without_library_prep_fails() {
    // The v1 contract: --library-prep is required, never silently defaulted.
    // clap rejects at parse time with exit code 2.
    let out = Command::new(binary())
        .args([
            "rna",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            "tests/output_phase0_missing_lp",
        ])
        .output()
        .expect("Failed to execute liquidqc rna");
    assert!(
        !out.status.success(),
        "liquidqc rna without --library-prep should fail; got status {:?}",
        out.status,
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("library-prep"),
        "Error should mention `library-prep`, got: {stderr:?}"
    );
}

#[test]
fn rna_help_lists_fragmentomics_flags() {
    let out = Command::new(binary())
        .args(["rna", "--help"])
        .output()
        .expect("Failed to execute liquidqc rna --help");
    assert!(out.status.success(), "liquidqc rna --help should exit 0");
    let stdout = String::from_utf8_lossy(&out.stdout);
    for flag in [
        "--library-prep",
        "--panels",
        "--snp-panel",
        "--min-gene-reads",
        "--fasta",
        "--single-end",
    ] {
        assert!(
            stdout.contains(flag),
            "rna --help should mention {flag}, got:\n{stdout}"
        );
    }
}

#[test]
fn rna_without_gtf_and_empty_cache_errors_clearly() {
    // BAM is synthetic (chr1/chr2 of 20kb) so the genome fingerprint
    // returns None — the error should mention `--gtf` and tell the user
    // that auto-detection failed, surfacing `fetch-references` so they
    // can populate the cache for a real genome.
    let outdir =
        std::env::temp_dir().join(format!("liquidqc-phase0-no-gtf-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&outdir);
    std::fs::create_dir_all(&outdir).unwrap();
    let cache = std::env::temp_dir().join(format!(
        "liquidqc-phase0-empty-cache-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&cache);
    std::fs::create_dir_all(&cache).unwrap();

    let out = Command::new(binary())
        .args([
            "rna",
            "tests/data/test.bam",
            "--outdir",
            outdir.to_str().unwrap(),
            "--library-prep",
            "unknown",
        ])
        .env("LIQUIDQC_CACHE_DIR", cache.to_str().unwrap())
        // Strip any inherited RUSTQC_GTF that would smuggle a default in.
        .env_remove("RUSTQC_GTF")
        .output()
        .expect("Failed to execute liquidqc rna");

    assert!(!out.status.success(), "expected non-zero exit");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--gtf") || stderr.contains("fetch-references"),
        "Error should suggest --gtf or fetch-references, got: {stderr:?}"
    );
}

#[test]
fn fetch_references_list_shows_pinned_genomes() {
    let out = Command::new(binary())
        .args(["fetch-references", "--list"])
        .output()
        .expect("Failed to execute liquidqc fetch-references --list");
    assert!(
        out.status.success(),
        "fetch-references --list should exit 0; stderr={}",
        String::from_utf8_lossy(&out.stderr),
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("GRCh38"),
        "expected GRCh38 in --list output, got:\n{stdout}"
    );
    assert!(
        stdout.contains("https://"),
        "expected an HTTPS URL in --list output, got:\n{stdout}"
    );
}

#[test]
fn fetch_references_unknown_genome_fails_with_message() {
    let out = Command::new(binary())
        .args(["fetch-references", "--genome", "NotAGenome"])
        .output()
        .expect("Failed to execute liquidqc fetch-references");
    assert!(
        !out.status.success(),
        "fetch-references with bogus --genome should fail"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("unknown genome") || stderr.contains("NotAGenome"),
        "Error should mention the unknown genome, got: {stderr:?}"
    );
}
