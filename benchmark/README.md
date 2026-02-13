# RustQC Benchmarks 🧬 🦀

Comparison of [dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor) and RustQC on the same input data.

In addition to dupRadar outputs, RustQC also generates featureCounts-compatible output files (counts TSV, summary, biotype counts, and MultiQC files) in the same pass.

## Small benchmark

A small test BAM file (`test.bam`) with a chr6-only GTF annotation, included in this repository.

### Results

| Metric | dupRadar (R) | RustQC |
| --- | --- | --- |
| **Runtime** | 2.50s | 0.25s (**10x**) |
| **Intercept** | 0.03186 | 0.03186 |
| **Slope** | 1.60189 | 1.60189 |
| **Genes total** | 2,905 | 2,905 |
| **Genes with reads** | 636 | 636 |
| **Genes with duplicates** | 201 | 201 |
| **Total values compared** | 37,765 | 37,765 |
| **Value mismatches** | — | **0** |

#### Count comparison

| Metric | dupRadar (R) | RustQC | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 20,449 | 20,449 | 100% |
| **filteredCounts (unique)** | 17,879 | 17,879 | 100% |
| **allCountsMulti** | 22,812 | 22,812 | **100%** |
| **filteredCountsMulti** | 20,034 | 20,034 | **100%** |

### Replication

```bash
# dupRadar (R)
Rscript benchmark/small/run_dupRadar_R.R

# RustQC (dupRadar + featureCounts outputs)
cargo run --release -- rna benchmark/small/test.bam --gtf benchmark/small/chr6.gtf -p --skip-dup-check -o benchmark/small/RustQC
```

The GTF for the small benchmark uses `gene_biotype` (Ensembl convention), so biotype counts are generated automatically.

## Large benchmark

**GM12878 REP1** — a full-size RNA-seq BAM from the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with picard.
Paired-end, unstranded, aligned to GRCh38 (Ensembl chromosome names).

### Input data

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

### Results

| Metric | dupRadar (R) | RustQC (1 thread) | RustQC (8 threads) | RustQC (10 threads) |
| --- | --- | --- | --- | --- |
| **Runtime** | 29m 56s | 3m 16s (~9x) | 1m 03s (~28x) | 0m 54s (~33x) |
| **Speedup** | — | **~9x** | **~28x** | **~33x** |
| **Max RSS** | N/A (Docker) | 503 MB | 893 MB | 1.3 GB |
| **Intercept** | 0.8245 | 0.8245 | 0.8245 | 0.8245 |
| **Slope** | 1.6774 | 1.6774 | 1.6774 | 1.6774 |
| **Genes total** | 63,086 | 63,086 | 63,086 | 63,086 |
| **Genes with reads (unique)** | 23,597 | 23,597 | 23,597 | 23,597 |
| **Genes with reads (multi)** | 24,719 | 24,719 | 24,719 | 24,719 |
| **Total values compared** | 820,118 | 820,118 | 820,118 | 820,118 |
| **Value mismatches** | — | **0** | **0** | **0** |

### Count comparison

| Metric | dupRadar (R) | RustQC | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 14,654,579 | 14,654,579 | **100%** |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | **100%** |
| **allCountsMulti** | 16,089,488 | 16,089,488 | **100%** |
| **filteredCountsMulti** | 4,503,920 | 4,503,920 | **100%** |

All four count columns match exactly across all 63,086 genes — both unique-mapper and multi-mapper counts. A cell-by-cell comparison of the full duplication matrix (820,118 values) shows **zero mismatches** (relative tolerance 1e-6).

Model fit parameters (**intercept** and **slope**) match to the displayed precision.

### Replication

#### 1. Download input files

```bash
mkdir -p benchmark/large

# Download BAM (~10 GB)
curl -L -o benchmark/large/GM12878_REP1.markdup.sorted.bam \
  "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam"

# Download GTF
curl -L -o benchmark/large/genes.gtf.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
gunzip benchmark/large/genes.gtf.gz

# Index BAM
samtools index benchmark/large/GM12878_REP1.markdup.sorted.bam
```

#### 2. Run dupRadar (R)

Requires R with `dupRadar` and `Rsubread` installed.

```bash
Rscript benchmark/large/run_dupRadar_R.R
```

#### 3. Run RustQC

The alignment file uses Ensembl chromosome names (`1`, `2`, ...) but the GENCODE GTF uses UCSC names (`chr1`, `chr2`, ...).
A config file is used to add the `chr` prefix to alignment chromosome names:

```yaml
# benchmark/large/config.yaml
chromosome_prefix: "chr"
```

```bash
cargo build --release

# Single-threaded
./target/release/rustqc rna \
  benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf \
  -p \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml \
  --biotype-attribute gene_type

# Multi-threaded (8 threads)
./target/release/rustqc rna \
  benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf \
  -p \
  -t 8 \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml \
  --biotype-attribute gene_type
```

> **Note:** The GENCODE GTF uses `gene_type` instead of the Ensembl default `gene_biotype`.
> Use `--biotype-attribute gene_type` to get biotype counts with GENCODE annotations.

### Known differences

Both benchmarks achieve **100% exact match** across all four count columns (unique and multi-mapper) and all genes. A cell-by-cell comparison of the full duplication matrix shows **zero mismatches** across all 820,118 values (37,765 for the small benchmark). Model fit parameters match to at least 10 significant digits.

Runtime may vary depending on hardware — the times above were measured on a single machine (10-core Apple Silicon) for relative comparison. The R dupRadar benchmark was run via Docker (x86 emulation on ARM Mac), so the R timings include container overhead and emulation penalties. Multi-threaded scaling depends on the number of chromosomes with mapped reads and the evenness of their read distribution.
