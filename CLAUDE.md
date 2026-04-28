# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

**liquidqc** is a hard fork of [seqeralabs/RustQC](https://github.com/seqeralabs/RustQC) extended into a single-pass cfRNA QC + fragmentomics tool. Binary crate (`liquidqc`), Rust edition 2021, GPL-3.0-or-later, MSRV 1.87.

Tier-1 (sample-level fragmentomics), Tier-2 (per-gene sparse Parquet), and Tier-3 (cohort QC + sample-identity add-ons) are all wired and emit into the v1 envelope. The inherited `rna` subcommand runs the full classical RNA-seq QC pipeline (dupRadar, featureCounts, 8 RSeQC tools, preseq, samtools-compatible flagstat/idxstats/stats, Qualimap gene body coverage) plus all fragmentomics features in one BAM pass and writes `<sample_id>.liquidqc.json` per sample. Schema is `1.0.0` (v1 stable).

The directory on disk is `cfqc/` but the crate, binary, and remote are all `liquidqc`. Don't "fix" this.

The `upstream` git remote points at seqeralabs/RustQC. Improvements that make sense upstream are routed back via cherry-picks — keep inherited per-file attribution intact. `AGENTS.md` is inherited from RustQC and still uses the RustQC name in places; treat its style/structure guidance as authoritative for inherited code.

## Commands

```bash
# Build (requires C++ compiler — build.rs compiles cpp/rng_shim.cpp for preseq RNG compat)
cargo build
cargo build --release    # LTO, strip, opt-level 3, codegen-units=1

# Format + lint (CI enforces both, warnings = errors)
cargo fmt --check
cargo clippy -- -D warnings

# Tests (CI runs in --release on Ubuntu + macOS)
cargo test
cargo test --test integration_test          # R-reference numerical comparisons (inherited)
cargo test --test phase0_cli                # CLI smoke tests for liquidqc subcommands
cargo test --test phase1_envelope           # per-sample envelope shape + qc_flags
cargo test --test phase2_fragmentomics      # Tier-1 fragmentomics blocks
cargo test --test phase3_per_gene           # Tier-2 per-gene Parquet
cargo test --test phase4_tier3              # Tier-3 cohort/identity add-ons
cargo test test_<substring>                 # single test by name
cargo test test_x -- --nocapture            # with stdout

# Pre-commit hooks (fmt + clippy + file hygiene)
prek install && prek run --all-files
```

Build script (`build.rs`) embeds `GIT_SHORT_HASH` and `BUILD_TIMESTAMP` as compile-time env vars and compiles the C++ RNG shim. System deps for `rust-htslib`: cmake, zlib, bz2, lzma, curl, ssl, clang, fontconfig.

## Architecture

Pure binary crate, no `lib.rs`. Top-level modules declared in `src/main.rs`; inter-module access uses `crate::` paths.

```
src/
  main.rs            — entry point, subcommand dispatch
  cli.rs             — clap-derive CLI; Commands = { Rna, Dna, Schema, Version }
  config.rs          — YAML config (mirrors CLI subcommand hierarchy under `rna:`)
  cpu.rs             — SIMD CPU compat check (runs before any auto-vectorized code to prevent SIGILL)
  io.rs              — gzip-transparent file reading (magic-byte detection)
  gtf.rs             — GTF parser; auto-detects `gene_type` (GENCODE) vs `gene_biotype` (Ensembl)
  citations.rs       — upstream tool versions for CITATIONS.md
  envelope.rs        — v1 envelope (`SampleEnvelope`, `Filters`, per-tool views, builder/writer)
  qc_flags.rs        — qc_flags rule engine (sample-level diagnostic flags)
  hash.rs            — md5 helpers for `bam_md5` / `gtf_md5` / `reference_fasta_md5` provenance
  runtime_stats.rs   — `peak_rss_mb` via getrusage (Linux/macOS unit handling)
  ui.rs              — terminal UI (verbosity, formatters)
  rna/               — inherited single-pass pipeline + fragmentomics; see AGENTS.md for the submodule map
schema/v1/liquidqc.schema.json   — v1 envelope schema (1.0.0)
cpp/rng_shim.cpp                 — std::mt19937 + binomial shim for preseq bootstrap reproducibility
```

The `rna` subcommand processes BAM files in a single pass, sharing a read dispatcher (`rna::rseqc::accumulators::RseqcAccumulators`) across all RSeQC-equivalent tools. Multiple BAMs run sequentially; `--threads` parallelizes within a sample. See `AGENTS.md` for the full breakdown of `rna/dupradar/`, `rna/featurecounts/`, `rna/qualimap/`, `rna/rseqc/`, and `rna/preseq.rs`.

## liquidqc-specific contracts (v1)

- **`--library-prep` is required and never silently defaults.** Pass `unknown` if truly unknown. Enforced by clap (no default) and part of the v1 contract — do not add a default value. Asserted by tests in `tests/phase0_cli.rs`.
- **`liquidqc dna` is reserved** — wired into the CLI but exits 1 in v1. v1 is cfRNA-only.
- **`liquidqc schema`** prints the embedded `schema/v1/liquidqc.schema.json` to stdout (via `include_str!`).
- **Output filenames:** `<sample_id>.liquidqc.json` (canonical per-sample envelope) and `<sample_id>.per_gene.parquet` (Tier-2 sparse per-gene rows, located/validated by the envelope's `per_gene` block). `sample_id` defaults to the BAM basename without extension; override with `--sample-id`. There is no run-level manifest.
- **Fragmentomics CLI flags** (`--panels`, `--panels-tsv`, `--snp-panel`, `--min-gene-reads`, `--sample-id`, `--saturation-fractions`) are load-bearing — they gate Tier-1/2/3 outputs in the envelope.

## Numerical accuracy

Inherited `rna` outputs are byte- and numerically-validated against upstream tools (R dupRadar, RSeQC 5.0.4, samtools, preseq). Float formatting follows R conventions (15 sig-figs, `NA` for NaN, trailing-zero trimming). `tests/expected/` is regenerated by `tests/create_test_data.R` — do not edit by hand. Preserve this when modifying the inherited pipeline.

## Style notes specific to this repo

- `unwrap()`/`expect()` are restricted to test code. In production, only when a prior guard makes it provably safe (with a comment).
- `anyhow::Result<T>` everywhere; no custom error types. Use `bail!`/`ensure!`/`.context()`.
- `IndexMap` when insertion order matters (gene ordering); `HashMap` otherwise.
- Package-level clippy allows in `Cargo.toml` (`manual_checked_ops`, `collapsible_match`, `unnecessary_sort_by`) exist to keep upstream cherry-picks clean — prefer adding to that list over editing inherited source.

See `AGENTS.md` for the full inherited style guide (formatting, imports, naming, derives, doc-comment conventions).

## Releases

Triggered via `workflow_dispatch` (GitHub Actions UI or `gh workflow run release`). Workflow validates `Cargo.toml` ↔ `Cargo.lock` version match, builds 8 binary targets + Docker images (baseline + 3 SIMD), then tags / publishes / uploads only after all builds pass. To prep: bump `Cargo.toml` version, `cargo update --package liquidqc`, add a `CHANGELOG.md` section, push to `main`, dispatch the workflow.
