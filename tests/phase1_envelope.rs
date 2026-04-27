//! Phase 1 envelope smoke tests.
//!
//! These run the binary on `tests/data/test.bam` and verify the generated
//! `<sample_id>.liquidqc.json` envelope has the expected shape. Numerical
//! parity for inherited tools is covered by `tests/integration_test.rs`.

use std::path::Path;
use std::process::Command;

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

/// Run liquidqc on the test BAM and return the parsed envelope JSON.
fn run_and_load_envelope(label: &str) -> serde_json::Value {
    let outdir =
        std::env::temp_dir().join(format!("liquidqc-phase1-{}-{}", label, std::process::id()));
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
        String::from_utf8_lossy(&out.stderr),
    );

    let envelope_path = outdir.join("test.liquidqc.json");
    let text = std::fs::read_to_string(&envelope_path)
        .unwrap_or_else(|e| panic!("envelope missing at {envelope_path:?}: {e}"));
    serde_json::from_str(&text).expect("envelope is not valid JSON")
}

#[test]
fn envelope_has_expected_top_level_shape() {
    let env = run_and_load_envelope("shape");
    let obj = env.as_object().expect("envelope must be a JSON object");

    // Identifiers + provenance.
    assert_eq!(obj["sample_id"].as_str(), Some("test"));
    assert_eq!(obj["library_prep"].as_str(), Some("unknown"));
    assert_eq!(obj["schema_version"].as_str(), Some("0.5.0-stub"));
    assert_eq!(
        obj["extractor_version"].as_str(),
        Some(env!("CARGO_PKG_VERSION"))
    );

    // BAM hash is 32 lowercase hex digits.
    let bam_md5 = obj["bam_md5"].as_str().expect("bam_md5 missing");
    assert_eq!(
        bam_md5.len(),
        32,
        "bam_md5 should be 32 hex chars: {bam_md5}"
    );
    assert!(
        bam_md5
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()),
        "bam_md5 should be lowercase hex: {bam_md5}"
    );
    assert!(obj["bam_size_bytes"].as_u64().unwrap() > 0);

    // Reference FASTA was not supplied, so reference_fasta_md5 is null.
    assert!(obj["reference_fasta_md5"].is_null());

    // Read-length stats derive from bam_stat.
    assert!(obj["read_length_max"].as_u64().unwrap() > 0);
    assert!(obj["read_length_mean"].as_f64().unwrap() > 0.0);

    // Filters block is fully populated with Phase 1 defaults.
    let filters = obj["filters"].as_object().expect("filters block");
    assert_eq!(filters["min_mapq"].as_u64(), Some(30));
    assert_eq!(filters["require_proper_pair"].as_bool(), Some(false));
    assert_eq!(filters["drop_secondary"].as_bool(), Some(true));
    assert_eq!(filters["drop_supplementary"].as_bool(), Some(true));
    assert_eq!(filters["drop_duplicates"].as_bool(), Some(false));

    // Runtime + RSS are non-negative.
    assert!(obj["runtime_seconds"].as_f64().unwrap() >= 0.0);
    assert!(obj["peak_rss_mb"].as_f64().unwrap() >= 0.0);
}

#[test]
fn envelope_qc_flags_include_expected_signals() {
    // The test BAM is single-end and unstranded by default — both deterministic
    // QC-flag triggers in Phase 1.
    let env = run_and_load_envelope("flags");
    let flags: Vec<&str> = env["qc_flags"]
        .as_array()
        .expect("qc_flags array")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();
    assert!(
        flags.contains(&"single_end_no_tlen"),
        "expected single_end_no_tlen in {flags:?}"
    );
    assert!(
        flags.contains(&"unstranded_endbias_unreliable"),
        "expected unstranded_endbias_unreliable in {flags:?}"
    );
}

#[test]
fn envelope_metric_blocks_present_when_tools_run() {
    let env = run_and_load_envelope("blocks");
    // bam_stat is auto-on by default and feeds the schema-required envelope
    // fields, so it must be present.
    assert!(env.get("bam_stat").is_some(), "bam_stat block missing");
    let bs = env["bam_stat"].as_object().unwrap();
    assert!(bs["total_records"].as_u64().unwrap() > 0);

    // featurecounts and dupradar are also default-on in the inherited config.
    assert!(
        env.get("featurecounts").is_some(),
        "featurecounts block missing"
    );

    // bam_stat.total_records should equal the envelope's read_count_total
    // (consistency check between the canonical field and the metric block).
    assert_eq!(
        env["read_count_total"].as_u64(),
        bs["total_records"].as_u64()
    );
}
