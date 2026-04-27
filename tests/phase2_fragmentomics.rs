//! Phase 2 Tier 1 fragmentomics integration tests.
//!
//! Exercises the four new accumulators (fragment-size bins, end-motif 4-mers,
//! soft-clip 4-mers, periodicity FFT) by running the binary on the inherited
//! `tests/data/test.bam` fixture under several invocation flavours, then
//! asserting the resulting envelope JSON has the expected blocks and qc_flags.

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

fn run_with(args: &[&str], label: &str) -> serde_json::Value {
    let outdir =
        std::env::temp_dir().join(format!("liquidqc-phase2-{}-{}", label, std::process::id()));
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
    serde_json::from_str(&text).expect("envelope is not valid JSON")
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
fn schema_version_bumped_to_phase4() {
    let env = run_with(&["--library-prep", "unknown"], "schema-version");
    assert_eq!(env["schema_version"].as_str(), Some("0.5.0-stub"));
}

#[test]
fn fragmentomics_block_emitted_without_fasta_single_end() {
    // Inherited test BAM is treated as single-end without --paired:
    //   - fragment_length and end_motifs blocks are absent (paired-only).
    //   - soft_clips block is present (read-anchored, no pairing required).
    //   - periodicity is absent (depends on TLEN histogram from fragment_size).
    //   - qc_flags includes single_end_no_tlen + end_motifs_skipped_no_fasta.
    let env = run_with(&["--library-prep", "unknown"], "no-fasta-se");
    let flags = qc_flags(&env);

    assert!(
        env.get("soft_clips").is_some(),
        "soft_clips block should be present"
    );
    assert!(
        env["soft_clips"]["primary_reads_observed"]
            .as_u64()
            .unwrap()
            > 0,
        "soft_clips primary_reads_observed should be > 0"
    );
    assert!(
        env.get("end_motifs").is_none(),
        "end_motifs block should be absent"
    );

    assert!(
        flags.contains(&"end_motifs_skipped_no_fasta"),
        "qc_flags should contain end_motifs_skipped_no_fasta; got {:?}",
        flags
    );
    assert!(
        flags.contains(&"single_end_no_tlen"),
        "qc_flags should contain single_end_no_tlen; got {:?}",
        flags
    );
}

#[test]
fn tagmentation_qc_flag_emitted() {
    let env = run_with(
        &[
            "--library-prep",
            "illumina_rna_prep_with_enrichment_tagmentation",
        ],
        "tagmentation",
    );
    let flags = qc_flags(&env);
    assert!(
        flags.contains(&"tagmentation_endmotif_bias"),
        "qc_flags should contain tagmentation_endmotif_bias; got {:?}",
        flags
    );
}

#[test]
fn read_length_caps_qc_flag_emitted() {
    // The inherited test BAM has reads under 100 nt — the threshold for
    // `read_length_caps_long_fragments`. This documents that mechanical
    // detection on the fixture is stable; if the fixture ever grows to
    // ≥100 nt reads this assertion will need updating.
    let env = run_with(&["--library-prep", "unknown"], "read-length-cap");
    let flags = qc_flags(&env);
    let read_length_max = env["read_length_max"].as_u64().unwrap_or(0);
    assert!(
        read_length_max > 0,
        "test BAM should produce a non-zero read_length_max"
    );
    if read_length_max < 100 {
        assert!(
            flags.contains(&"read_length_caps_long_fragments"),
            "read_length_max={} < 100 should set read_length_caps_long_fragments; got {:?}",
            read_length_max,
            flags
        );
    }
}

#[test]
fn schema_subcommand_includes_phase2_blocks_and_flag() {
    // The schema embedded in the binary must mention every new Phase 2 block
    // and the new qc_flag, so callers consuming the schema can validate
    // outputs without reading the source file.
    let out = Command::new(binary())
        .args(["schema"])
        .output()
        .expect("Failed to execute liquidqc schema");
    assert!(out.status.success(), "liquidqc schema failed");
    let parsed: serde_json::Value =
        serde_json::from_slice(&out.stdout).expect("schema is not valid JSON");

    let title = parsed["title"].as_str().unwrap_or("");
    assert!(
        title.contains("Phase 4"),
        "Schema title should mention Phase 4; got {:?}",
        title
    );

    let props = parsed["properties"]
        .as_object()
        .expect("schema properties missing");
    for block in [
        "fragment_length",
        "end_motifs",
        "soft_clips",
        "periodicity",
        "per_gene",
    ] {
        assert!(
            props.contains_key(block),
            "schema properties should define `{block}`; keys = {:?}",
            props.keys().collect::<Vec<_>>()
        );
    }

    // Phase 4 + bug-fix pass: schema_version examples should include 0.5.0-stub.
    let version_examples: Vec<&str> = parsed["properties"]["schema_version"]["examples"]
        .as_array()
        .expect("schema_version examples missing")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();
    assert!(
        version_examples.contains(&"0.5.0-stub"),
        "schema_version examples should include 0.5.0-stub; got {:?}",
        version_examples
    );

    let enum_values: Vec<&str> = parsed["properties"]["qc_flags"]["items"]["enum"]
        .as_array()
        .expect("qc_flags enum missing")
        .iter()
        .filter_map(|v| v.as_str())
        .collect();
    assert!(
        enum_values.contains(&"end_motifs_skipped_no_fasta"),
        "qc_flags enum should contain end_motifs_skipped_no_fasta; got {:?}",
        enum_values
    );
}

#[test]
fn old_reference_flag_is_rejected() {
    // The CLI rename from --reference to --fasta is part of Phase 2's
    // contract. clap should reject the legacy spelling at parse time
    // (exit 2). We do not assert any backward-compat alias.
    let out = Command::new(binary())
        .args([
            "rna",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--library-prep",
            "unknown",
            "--reference",
            "/nonexistent.fa",
        ])
        .output()
        .expect("Failed to execute liquidqc rna --reference");
    assert!(
        !out.status.success(),
        "--reference should be rejected; got status {:?}",
        out.status
    );
}
