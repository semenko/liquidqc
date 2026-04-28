<p align="center">
  <img src="docs/public/liquidqc-logo.svg" alt="liquidqc logo" width="180" />
</p>

# liquidqc

**Single-pass cfRNA QC + fragmentomics in Rust.**

Created by [Nick Semenkovich](https://semenko.com).

[![CI](https://github.com/semenko/liquidqc/actions/workflows/ci.yml/badge.svg)](https://github.com/semenko/liquidqc/actions/workflows/ci.yml)
[![Release](https://github.com/semenko/liquidqc/actions/workflows/release.yml/badge.svg)](https://github.com/semenko/liquidqc/actions/workflows/release.yml)
[![Docs](https://github.com/semenko/liquidqc/actions/workflows/build_docs.yml/badge.svg)](https://github.com/semenko/liquidqc/actions/workflows/build_docs.yml)
[![Crates.io](https://img.shields.io/crates/v/liquidqc.svg)](https://crates.io/crates/liquidqc)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)
[![Rust 1.87+](https://img.shields.io/badge/rust-1.87%2B-orange.svg)](https://www.rust-lang.org)
[![Container](https://img.shields.io/badge/container-ghcr.io%2Fsemenko%2Fliquidqc-blue)](https://github.com/semenko/liquidqc/pkgs/container/liquidqc)

`liquidqc` ingests a sorted, duplicate-marked, indexed genome BAM (plus a GTF
and, for end-motif features, a reference FASTA) and emits one versioned JSON
envelope per sample alongside a sparse per-gene Parquet — covering classical
RNA-seq QC (RSeQC + samtools stats + Qualimap + dupRadar + featureCounts +
preseq equivalents) and cfRNA fragmentomics (end-motif k-mers, soft-clip
k-mers, fragment-length periodicity, per-gene Tier-2 features, panels, sex
inference, SNP fingerprint, saturation curves) in **one BAM iteration**.

It is a hard fork of [seqeralabs/RustQC](https://github.com/seqeralabs/RustQC),
and the inherited classical-QC outputs remain byte- and numerically identical
to the upstream tools.

## Features

- **One BAM pass.** A shared read dispatcher feeds dupRadar, featureCounts,
  eight RSeQC tools, Qualimap, samtools stats, preseq, and all fragmentomics
  accumulators. `--threads` parallelizes within a sample; multiple input BAMs
  run sequentially.
- **One versioned output.** A single `<sample_id>.liquidqc.json` envelope per
  sample, with provenance hashes (`bam_md5`, `gtf_md5`, `reference_fasta_md5`),
  filter parameters, and a `qc_flags` rule engine. Companion
  `<sample_id>.per_gene.parquet` holds Tier-2 per-gene rows. Schema published
  as JSON Schema 2020-12 at `schema/v1/liquidqc.schema.json` (`1.0.0`).
- **MultiQC-compatible legacy outputs.** Inherited per-tool TSV/text/PNG files
  remain identical to the upstream tools and parse cleanly in
  [MultiQC](https://multiqc.info/).
- **Static binary.** SAM/BAM/CRAM input via `rust-htslib`; SIMD-optimized
  builds shipped with CPU-feature detection (baseline + 3 SIMD tiers).
- **Reproducible.** Per-tool seeds (`--tin-seed`, `--junction-saturation-seed`,
  `--preseq-seed`); preseq bootstrap RNG bridged to upstream's `std::mt19937`
  via a small C++ shim.

### Tools reimplemented (and numerically validated)

| Tool                | Reimplements                                                                            | Description                                                           |
| ------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| dupRadar            | [dupRadar](https://github.com/ssayols/dupRadar)                                         | PCR duplicate rate vs. expression; density scatter plots              |
| featureCounts       | [featureCounts](http://subread.sourceforge.net/)                                        | Gene-level read counting with biotype summaries                       |
| bam_stat            | [RSeQC](https://rseqc.sourceforge.net/#bam-stat-py) `bam_stat.py`                       | Basic alignment statistics                                            |
| infer_experiment    | [RSeQC](https://rseqc.sourceforge.net/#infer-experiment-py) `infer_experiment.py`       | Library strandedness inference                                        |
| read_duplication    | [RSeQC](https://rseqc.sourceforge.net/#read-duplication-py) `read_duplication.py`       | Position- and sequence-based duplication histograms                   |
| read_distribution   | [RSeQC](https://rseqc.sourceforge.net/#read-distribution-py) `read_distribution.py`     | Read distribution across genomic features                             |
| junction_annotation | [RSeQC](https://rseqc.sourceforge.net/#junction-annotation-py) `junction_annotation.py` | Splice junction classification                                        |
| junction_saturation | [RSeQC](https://rseqc.sourceforge.net/#junction-saturation-py) `junction_saturation.py` | Splice junction saturation analysis                                   |
| inner_distance      | [RSeQC](https://rseqc.sourceforge.net/#inner-distance-py) `inner_distance.py`           | Paired-end inner distance distribution                                |
| TIN                 | [RSeQC](https://rseqc.sourceforge.net/#tin-py) `tin.py`                                 | Transcript Integrity Number                                           |
| preseq              | [preseq](http://smithlabresearch.org/software/preseq/) `lc_extrap`                      | Library complexity extrapolation                                      |
| Qualimap rnaseq     | [Qualimap](http://qualimap.conesalab.org/) `rnaseq`                                     | Gene body coverage, read origin, strand specificity                   |
| flagstat            | [samtools](http://www.htslib.org/) `flagstat`                                           | Alignment flag summary                                                |
| idxstats            | [samtools](http://www.htslib.org/) `idxstats`                                           | Per-chromosome read counts                                            |
| stats               | [samtools](http://www.htslib.org/) `stats`                                              | Full samtools stats output (SN section + all histogram sections)      |

### cfRNA-specific output blocks

- **Sample-level fragmentomics:** end-motif 4-mers (5′ and 3′), soft-clip
  k-mers, fragment-length periodicity FFT, short/long fragment fractions
  with binned histogram, adapter-readthrough rate.
- **Per-gene sparse Parquet:** length histogram, end motifs, mean/median
  fragment length, 5′-3′ coverage bias, soft-clip rates, CDS/UTR/intron
  retention — restricted to genes with ≥ `--min-gene-reads` reads.
- **Cohort QC + sample identity:** per-chromosome read fractions (chrY,
  chrM), hemoglobin and ribosomal-protein gene fractions, marker-panel
  reads (LM22, Tabula Sapiens cfRNA, Vorperian + user TSVs), sex inference
  from XIST/RPS4Y1, per-cycle quality, adapter readthrough, splice-site
  dinucleotide composition, gene-detection saturation curves, ~17.5k-site
  SNP fingerprint, automatic `qc_flags`.

## Install

### From source

```bash
git clone https://github.com/semenko/liquidqc
cd liquidqc
cargo build --release
./target/release/liquidqc --version
```

Requires Rust **1.87+**. `build.rs` compiles a small C++ RNG shim, so a
C++ toolchain is needed alongside the usual `rust-htslib` system deps:

- **macOS:** `brew install cmake xz bzip2`
- **Debian/Ubuntu:** `apt-get install cmake clang pkg-config libssl-dev
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libfontconfig1-dev`

### Docker

```bash
docker pull ghcr.io/semenko/liquidqc:latest
docker run --rm -v "$PWD:/data" ghcr.io/semenko/liquidqc:latest \
  rna /data/sample.bam --gtf /data/genes.gtf.gz \
      --library-prep neb_next_ultra_ii_directional --paired \
      --outdir /data/results
```

Pre-built Linux/macOS binaries (with SIMD tiers) are attached to each
[GitHub Release](https://github.com/semenko/liquidqc/releases).

## Quickstart

```bash
liquidqc rna sample.markdup.bam \
  --gtf gencode.v44.gtf.gz \
  --fasta GRCh38.primary_assembly.genome.fa \
  --library-prep neb_next_ultra_ii_directional \
  --paired \
  --outdir results/
```

This writes:

```
results/
├── sample.liquidqc.json          # canonical v1 envelope (always)
├── sample.per_gene.parquet       # Tier-2 sparse per-gene rows
├── CITATIONS.md                  # upstream tool versions
└── ...                           # MultiQC-compatible legacy outputs
                                  # (dupradar/, featurecounts/, rseqc/,
                                  #  qualimap/, preseq/, samtools/)
```

### Programmatic access to the envelope

```bash
# Pull just the qc_flags section
jq '.qc_flags' results/sample.liquidqc.json

# Inspect the per-gene Parquet without loading the whole file
duckdb -c "SELECT gene_id, n_reads, mean_fragment_length
           FROM 'results/sample.per_gene.parquet'
           ORDER BY n_reads DESC LIMIT 20;"

# Print the schema (embedded at compile time)
liquidqc schema | jq '.properties | keys'
```

### Multiple BAMs, single GTF parse

```bash
liquidqc rna A.bam B.bam C.bam \
  --gtf gencode.v44.gtf.gz \
  --fasta GRCh38.primary_assembly.genome.fa \
  --library-prep unknown \
  --paired \
  --threads 8 \
  --outdir results/
```

The GTF is parsed once; samples are processed sequentially, with
`--threads` parallelizing the htslib read decode within each sample.
Each BAM gets its own `<sample_id>.liquidqc.json` and
`<sample_id>.per_gene.parquet`. `sample_id` defaults to the BAM basename
(without extension); override with `--sample-id`.

### Stranded RNA-seq

```bash
liquidqc rna sample.bam \
  --gtf genes.gtf.gz \
  --fasta GRCh38.fa \
  --stranded reverse \
  --library-prep takara_smarter_v3 \
  --paired \
  --outdir results/
```

`--stranded` accepts `unstranded`, `forward`, or `reverse`. If omitted,
strandedness is inferred from the data and reported in
`infer_experiment.tsv`.

### Custom marker panels

```bash
liquidqc rna sample.bam \
  --gtf genes.gtf.gz \
  --fasta GRCh38.fa \
  --library-prep unknown \
  --panels lm22,tabula_sapiens_cfrna \
  --panels-tsv my_neutrophil_signature.tsv \
  --paired \
  --outdir results/
```

Bundled panels: `lm22`, `tabula_sapiens_cfrna`, `vorperian`. Pass an empty
value (`--panels ""`) to disable bundled panels. User TSVs follow the
columns `gene_id<TAB>gene_symbol<TAB>panel<TAB>cell_type<TAB>weight`.

## CLI reference

| Subcommand          | Purpose                                                    |
| ------------------- | ---------------------------------------------------------- |
| `liquidqc rna`      | Single-pass cfRNA QC + fragmentomics                       |
| `liquidqc dna`      | Reserved for cfDNA; exits 1 in v1                          |
| `liquidqc schema`   | Print the embedded v1 JSON Schema to stdout                |
| `liquidqc version`  | Print version + git short hash + build timestamp           |

The full flag list is available via `liquidqc rna --help`. A few worth
calling out:

| Flag                                | Notes                                                                                                  |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `--library-prep <STRING>`           | **Required.** Never silently defaulted; pass `unknown` if truly unknown. Part of the v1 contract.      |
| `--fasta <FASTA>`                   | Required for CRAM input and for end-motif fragmentomics; other features run without it.                |
| `--paired`                          | Pass for paired-end libraries; affects insert-size, inner-distance, and fragment-size accumulators.    |
| `--min-gene-reads <N>`              | Per-gene Parquet inclusion threshold (default 20).                                                     |
| `--snp-panel <VCF>`                 | Common-SNP fingerprint VCF (default: bundled somalier sites).                                          |
| `--saturation-fractions <CSV>`      | Saturation-curve subsample fractions (default `0.05,0.10,0.25,0.50,0.75,1.00`).                        |
| `--config <FILE>`                   | YAML config (XDG-discoverable) mirroring the CLI hierarchy under `rna:`.                               |

## Output schema

The per-sample envelope is the canonical output. Top-level fields include
provenance (`extractor_version`, `schema_version`, `git_commit`,
`bam_md5`, `gtf_md5`, `reference_fasta_md5`), filter parameters, the
inherited classical-QC blocks (`bam_stat`, `infer_experiment`,
`read_distribution`, `read_duplication`, `junction_annotation`,
`junction_saturation`, `inner_distance`, `tin`, `preseq`, `qualimap`,
`dupradar`, `featurecounts`), the cfRNA blocks (`fragment_length`,
`end_motifs`, `soft_clips`, `periodicity`, `per_gene`,
`adapter_readthrough_rate`, `splice_site_dinucleotides`, `per_chromosome`,
`gene_class_fractions`, `panels`, `sex_inference`, `snp_fingerprint`,
`saturation`, `cycle_quality`), `qc_flags`, and runtime stats
(`runtime_seconds`, `peak_rss_mb`).

Schema: [`schema/v1/liquidqc.schema.json`](schema/v1/liquidqc.schema.json)
(JSON Schema 2020-12). The v1 release will bump the schema version to
`1.0.0`. Per-version notes are in [`CHANGELOG.md`](CHANGELOG.md).

## Numerical accuracy

Inherited classical-QC outputs are validated byte- and numerically against
upstream tools:

- **dupRadar:** R reference fixtures regenerated by `tests/create_test_data.R`
- **RSeQC 5.0.4:** integration tests compare TSV outputs and histograms
- **samtools:** flagstat/idxstats/stats compared against `samtools` of the
  same htslib version
- **preseq:** bootstrap reproducibility via a `std::mt19937` C++ shim

Float formatting follows R conventions (15 sig-figs, `NA` for NaN,
trailing-zero trimming). Do not hand-edit `tests/expected/`; regenerate
via the R script.

## License & attribution

GPL-3.0-or-later, inherited from RustQC. See [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE) for the upstream attribution. Per-source-file
attribution headers are intentionally preserved in inherited modules.
Improvements that make sense upstream are routed back to RustQC via
cherry-picks; the `upstream` git remote points there.

RustQC, and therefore liquidqc, was developed with substantial AI
assistance and validated by output comparison against the upstream tools
on real sequencing data — see the upstream
[AI & Provenance](https://seqeralabs.github.io/RustQC/about/ai-statement/)
notes.

## Citing

If `liquidqc` was useful for your work, please cite both this project
and the upstream tools whose outputs it reproduces. A
`CITATIONS.md` listing the upstream tool versions used in your run is
written to the output directory automatically.
