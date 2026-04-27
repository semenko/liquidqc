# liquidqc Changelog

## [Unreleased] — Phase 1 (per-sample envelope) — 2026-04-26

### v1 envelope wiring

- The canonical per-sample output is now `<sample_id>.liquidqc.json`,
  written into `--outdir` for every BAM processed. The envelope follows
  `schema/v1/liquidqc.schema.json` and embeds inherited accumulator
  results as optional metric blocks (`bam_stat`, `infer_experiment`,
  `read_distribution`, `read_duplication`, `junction_annotation`,
  `junction_saturation`, `inner_distance`, `tin`, `preseq`, `qualimap`,
  `dupradar`, `featurecounts`).
- Added `src/envelope.rs` with `SampleEnvelope`, `Filters`, the per-tool
  view structs, `SCHEMA_VERSION = "0.1.0-stub"`, and a builder/writer
  pair. View structs decouple the public schema from internal struct
  drift; each has a `from_result(...)` constructor.
- Added `src/qc_flags.rs` (Phase 1 rule engine). Wires
  `single_end_no_tlen`, `low_paired_fraction`, `low_mapping_rate`,
  `high_rrna_fraction`, `tagmentation_endmotif_bias`,
  `unstranded_endbias_unreliable`, `read_length_caps_long_fragments`,
  `low_complexity_library`, and `high_duplication_nonbiological`.
  `cfdna_contamination_suspected` (needs fragmentomics) and
  `sex_swap_warning` (needs user metadata) keep their enum variants but
  are deferred to Phases 2+ / 4+ respectively.
- Added `src/hash.rs` (`md5_and_size`, `md5`, `md5_uncompressed`) for
  envelope provenance fields. `bam_md5` is the BGZF bytes as on disk;
  `gtf_md5` is the **uncompressed** bytes (gzip-transparent);
  `reference_fasta_md5` is the file as supplied. Per-BAM hashes run in
  parallel via rayon before per-BAM processing.
- Added `src/runtime_stats.rs` (`peak_rss_mb`) using `getrusage` with
  Linux/macOS unit-handling. Reported value is process-level — every
  per-sample envelope in a multi-BAM run sees the same final
  high-water mark.
- Refactored `write_rseqc_outputs` to bubble out every accumulator
  result via an extended `RseqcOutputs` struct so the envelope writer
  can read them. No external behavior change for inherited TSV/text/PNG
  outputs.

### Removed (Phase 0 → Phase 1, breaking)

- The legacy per-run `liquidqc_summary.json` writer and the
  `--json-summary` / `RUSTQC_JSON_SUMMARY` flag. Per-sample envelopes
  are the canonical output; no run-level manifest is emitted.
- `tests/snapshot/empty_envelope.json` — replaced by a live-build
  envelope check in `tests/phase0_cli.rs`.
- `src/summary.rs` (RunSummary / InputSummary / CountingSummary /
  DupradarSummary / OutputFile structs) — superseded by `src/envelope.rs`.

### Tests

- Added `tests/phase1_envelope.rs` with smoke tests for envelope shape,
  QC-flag wiring, and metric-block presence on the inherited test BAM.
- `tests/phase0_cli.rs::live_envelope_satisfies_schema_required_fields`
  replaces the static fixture-based check.
- All inherited R-parity tests (`tests/integration_test.rs`) continue
  to pass byte-for-byte.

### Build

- Added direct `md-5 = "0.10"` and `libc = "0.2"` dependencies. Enabled
  the `serde` feature on `indexmap` so biotype-keyed maps serialize
  cleanly.

## [0.1.0-bootstrap (Phase 0)] — 2026-04-26

### Project

- Hard-forked [seqeralabs/RustQC](https://github.com/seqeralabs/RustQC)
  v0.2.1 → renamed to **liquidqc**, version reset to 0.1.0. RustQC's git
  history is preserved; the `upstream` remote points at the original
  project for cherry-picks. License remains GPL-3.0-or-later.
- Updated `Cargo.toml` (crate name, binary name, description, repository,
  authors), `cli.rs` clap command name, `main.rs` doc-comments and the
  output JSON filename (`liquidqc_summary.json`), and `release.yml`
  binary-name and lockfile-pattern references. Inherited per-source-file
  attribution to RustQC remains in docstrings, build script, citations
  module, and integration-test docstrings.
- Added `NOTICE` with full RustQC attribution per GPLv3.
- Rewrote `README.md` for liquidqc identity while preserving upstream
  credit prominently.

### Schema v1 stub

- Added `schema/v1/liquidqc.schema.json` with the always-present envelope
  fields (extractor/schema/git versions, BAM/GTF/FASTA hashes, library
  prep, paired-end, strandedness, read-length and read-count summaries,
  filter parameters, qc_flags, runtime/RSS). Per-metric blocks land in
  Phase 1+. Schema marked `0.1.0-stub`.
- Added `tests/snapshot/empty_envelope.json` reference fixture and a
  parseability smoke test.

### CLI scaffolding (parse-only)

- Added top-level subcommands `dna`, `version`, `schema` alongside the
  inherited `rna`. `dna` prints "not implemented" and exits 1 (the
  v1-cfRNA-only contract). `schema` emits the v1 schema to stdout.
- Added `rna` flag stubs for the new fragmentomics CLI surface
  (`--library-prep`, `--paired-end auto|true|false`, `--panels`,
  `--snp-panel`, `--min-gene-reads`, `--sample-id`, `--out-dir`).
  `--library-prep` is required and never silently defaults — exits
  non-zero if omitted, even though the underlying computation isn't yet
  wired.

## [Version 0.2.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.1) - 2026-04-09

### Bug fixes

- Fix SIMD builds by scoping rustflags to target triple with explicit `--target` (#90, #91)

### Other changes

- Trigger releases via workflow dispatch instead of tag push (#92)

## [Version 0.2.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.0) - 2026-04-09

### Features

- Ship SIMD-optimized binaries with CPU detection and upgrade hints (#81)
- Write `CITATIONS.md` with upstream tool versions (#87)
- Add XDG config discovery and deep-merge support (#88)

### Bug fixes

- Replace header-based duplicate check with flag-based detection (#84)
- Use `.log` extension for junction_annotation output (#80)

### Other changes

- Bump docker/login-action from 4.0.0 to 4.1.0 (#78)
- Bump vite from 7.3.1 to 7.3.2 in docs (#77)
- Bump defu from 6.1.4 to 6.1.6 in docs (#74)

## [Version 0.1.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.1) - 2026-04-02

### Bug fixes

- Fix featureCounts summary to use gene-level stats; add biotype summary (#66)
- Fix inner_distance histogram to include overflow bucket in bulk cutoff loop (#67)

### Other changes

- Add crates.io publishing to release workflow (#62)
- Documentation fixes (#70)

## [Version 0.1.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.0) - 2026-04-01

Initial release of RustQC -- fast quality control tools for sequencing data, written in Rust.

A single `rustqc rna` command runs 15 QC analyses in one pass over the BAM file, producing output that is format- and numerically identical to the upstream tools and fully compatible with [MultiQC](https://multiqc.info/).

### Tools

- **dupRadar** -- PCR duplicate rate vs. expression analysis with density scatter plots, boxplots, and expression histograms. 14-column duplication matrix with logistic-regression model fitting matching the [R dupRadar](https://github.com/ssayols/dupRadar) package.
- **featureCounts** -- Gene-level read counting with assignment summary, compatible with [Subread featureCounts](http://subread.sourceforge.net/). Includes per-biotype read counting and MultiQC integration.
- **RSeQC** (8 tools) -- [RSeQC](https://rseqc.sourceforge.net/)-compatible implementations of bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance, and TIN (Transcript Integrity Number). Includes native plot generation (PNG + SVG) with no R dependency.
- **preseq** -- Library complexity extrapolation (`lc_extrap`) matching [preseq](http://smithlabresearch.org/software/preseq/) v3.2.0, including C++ RNG FFI for reproducible bootstrap sampling.
- **Qualimap rnaseq** -- Gene body coverage profiling, read genomic origin, junction analysis, and transcript coverage bias matching [Qualimap](http://qualimap.conesalab.org/).
- **samtools** -- flagstat, idxstats, and full stats output (SN section + all histogram sections) matching [samtools](http://www.htslib.org/).

### Features

- Single static binary with no runtime dependencies
- SAM, BAM, and CRAM input support (auto-detected)
- Multi-threaded alignment processing with htslib thread pool
- GTF annotation support (gzip-compressed files accepted)
- YAML configuration for output control, chromosome name mapping, and tool parameters
- Multiple BAM file support via positional arguments
- `--sample-name` flag to override BAM-derived sample name in output filenames
- Per-tool seed flags (`--tin-seed`, `--junction-saturation-seed`, `--preseq-seed`) for reproducible results
- Docker container at `ghcr.io/seqeralabs/rustqc`
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
