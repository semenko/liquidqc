# RustQC Benchmarks

Benchmark data and scripts for comparing RustQC against
[dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor),
[Subread featureCounts](http://subread.sourceforge.net/), and
[RSeQC](https://rseqc.sourceforge.net/).

For detailed results, tables, and side-by-side plot comparisons, see the
documentation:

- [dupRadar benchmarks](https://ewels.github.io/RustQC/benchmarks/dupradar/)
- [featureCounts benchmarks](https://ewels.github.io/RustQC/benchmarks/featurecounts/)

## Benchmark data

### Small benchmark

Included in this repository. A test BAM file with chr6 reads, a chr6-only GTF
annotation (2,905 genes), and a BED12 gene model for RSeQC tools.

- `small/test.bam` + `small/chr6.gtf` + `small/chr6.bed`
- `small/dupRadar/` — R dupRadar reference output
- `small/RustQC/` — RustQC output

### Large benchmark

GM12878 REP1 — a full-size RNA-seq BAM (~10 GB) from the
[nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with
Picard. Paired-end, unstranded, aligned to GRCh37 (Ensembl chromosome names).

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

The dupRadar/featureCounts pipeline uses a GENCODE v46 GTF (`genes.gtf`). The
RSeQC tools use a BED12 gene model (`genes.bed`) converted from a
genome-matched GTF with matching Ensembl chromosome names.

### RSeQC reference outputs

Reference outputs for validating RSeQC reimplementations are stored in
`rseqc/small/` and `rseqc/large/`. These were generated with
[RSeQC 5.0.4](https://rseqc.sourceforge.net/) run via Docker
(`--platform linux/amd64`).

**Small BAM** — All 7 tools have reference output:

| Tool | Key metrics |
| ---- | ----------- |
| `bam_stat` | 43,476 total reads, 6,032 duplicates |
| `infer_experiment` | Undetermined (50/50 strand split) |
| `read_duplication` | Position-based and sequence-based histograms |
| `read_distribution` | 68,660 tags, 66,693 assigned (49,589 CDS) |
| `junction_annotation` | 3,261 junctions (2,982 known, 88 partial novel, 191 novel) |
| `junction_saturation` | Saturation curve, 20 sampling points |
| `inner_distance` | 20,861 pairs (mean -38.85, median -72.5) |

**Large BAM** — 6 of 7 tools have reference output (`read_duplication`
pending):

| Tool | Key metrics |
| ---- | ----------- |
| `bam_stat` | 185,718,543 total records, 133,912,519 duplicates, 39,827,099 unique (MAPQ>=30) |
| `infer_experiment` | Reverse stranded (92.2% / 1.2%) |
| `read_distribution` | 55,374,023 tags, 52,400,513 assigned (33.3M CDS, 8.5M 3'UTR) |
| `junction_annotation` | 256,466 junctions (178,797 known, 50,936 partial novel, 26,733 novel) |
| `junction_saturation` | 256,466 junctions at 100% (163,710 known, 92,756 novel) |
| `inner_distance` | 1,000,000 pairs sampled (mean 29.43, median 27.5, SD 32.80) |

## Reproducing benchmarks

### 1. Download large benchmark input files

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

### 2. Generate BED12 for RSeQC tools

Convert the GTF to BED12 format using UCSC tools. This can be done via Docker
with the [wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda ucsc-gtftogenepred ucsc-genepredtobed

# Then run the conversion (substitute the wave image name):
docker run --rm -v $(pwd)/benchmark/large:/data <wave-image> \
  bash -c 'gtfToGenePred /data/genes_matched.gtf /tmp/genes.genePred && \
           genePredToBed /tmp/genes.genePred /data/genes.bed'
```

### 3. Run dupRadar (R)

Requires R with `dupRadar` and `Rsubread` installed.

```bash
# Small
Rscript benchmark/small/run_dupRadar_R.R

# Large
Rscript benchmark/large/run_dupRadar_R.R
```

### 4. Run RustQC

```bash
cargo build --release

# Small
./target/release/rustqc rna benchmark/small/test.bam \
  --gtf benchmark/small/chr6.gtf -p --skip-dup-check \
  -o benchmark/small/RustQC

# Large (the GENCODE GTF uses UCSC chrom names while the BAM uses Ensembl names)
./target/release/rustqc rna benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf -p -t 10 \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml \
  --biotype-attribute gene_type
```

> **Note:** The config file adds a `chr` prefix to alignment chromosome names
> to match the GENCODE GTF. The `--biotype-attribute gene_type` flag is needed
> because GENCODE uses `gene_type` instead of the Ensembl default `gene_biotype`.

### 5. Run RSeQC tools (RustQC reimplementations)

```bash
# Small — all tools
./target/release/rustqc bam-stat benchmark/small/test.bam -o benchmark/small/RustQC
./target/release/rustqc infer-experiment benchmark/small/test.bam -b benchmark/small/chr6.bed -o benchmark/small/RustQC
./target/release/rustqc read-duplication benchmark/small/test.bam -o benchmark/small/RustQC
./target/release/rustqc read-distribution benchmark/small/test.bam -b benchmark/small/chr6.bed -o benchmark/small/RustQC
./target/release/rustqc junction-annotation benchmark/small/test.bam -b benchmark/small/chr6.bed -o benchmark/small/RustQC
./target/release/rustqc junction-saturation benchmark/small/test.bam -b benchmark/small/chr6.bed -o benchmark/small/RustQC
./target/release/rustqc inner-distance benchmark/small/test.bam -b benchmark/small/chr6.bed -o benchmark/small/RustQC

# Large — same commands with large BAM + BED
./target/release/rustqc bam-stat benchmark/large/GM12878_REP1.markdup.sorted.bam -o benchmark/large/RustQC
./target/release/rustqc infer-experiment benchmark/large/GM12878_REP1.markdup.sorted.bam -b benchmark/large/genes.bed -o benchmark/large/RustQC
./target/release/rustqc read-distribution benchmark/large/GM12878_REP1.markdup.sorted.bam -b benchmark/large/genes.bed -o benchmark/large/RustQC
./target/release/rustqc junction-annotation benchmark/large/GM12878_REP1.markdup.sorted.bam -b benchmark/large/genes.bed -o benchmark/large/RustQC
./target/release/rustqc junction-saturation benchmark/large/GM12878_REP1.markdup.sorted.bam -b benchmark/large/genes.bed -o benchmark/large/RustQC
./target/release/rustqc inner-distance benchmark/large/GM12878_REP1.markdup.sorted.bam -b benchmark/large/genes.bed -o benchmark/large/RustQC
```

### 6. Generate RSeQC Python reference outputs

Requires RSeQC 5.0.4 via Docker:

```bash
RSEQC_IMG="wave.seqera.io/wt/ea3e9f972b6e/wave/build:rseqc-5.0.4--14c99cde3bff8d57"

# bam_stat
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data \
  $RSEQC_IMG bam_stat.py -i /data/test.bam > benchmark/rseqc/small/bam_stat.txt 2>&1

# infer_experiment
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data \
  $RSEQC_IMG infer_experiment.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/rseqc/small/infer_experiment.txt 2>&1

# read_duplication
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data -v $(pwd)/benchmark/rseqc/small:/out \
  $RSEQC_IMG read_duplication.py -i /data/test.bam -o /out/read_duplication

# read_distribution
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data \
  $RSEQC_IMG read_distribution.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/rseqc/small/read_distribution.txt 2>&1

# junction_annotation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data -v $(pwd)/benchmark/rseqc/small:/out \
  $RSEQC_IMG junction_annotation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_annotation

# junction_saturation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data -v $(pwd)/benchmark/rseqc/small:/out \
  $RSEQC_IMG junction_saturation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_saturation

# inner_distance
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/small:/data -v $(pwd)/benchmark/rseqc/small:/out \
  $RSEQC_IMG inner_distance.py -i /data/test.bam -r /data/chr6.bed -o /out/inner_distance
```

### 7. Compare results

```bash
# Compare duplication matrices cell-by-cell
python3 -c "
import csv
with open('benchmark/large/dupRadar/dupMatrix.txt') as rf, \
     open('benchmark/large/RustQC/GM12878_REP1.markdup.sorted_dupMatrix.txt') as rustf:
    r = list(csv.reader(rf, delimiter='\t'))
    rust = list(csv.reader(rustf, delimiter='\t'))
    mismatches = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i]))
                     if r[i][j] != 'NA' and rust[i][j] != 'NA'
                     and abs(float(r[i][j]) - float(rust[i][j])) / max(abs(float(r[i][j])), 1e-15) > 1e-6)
    total = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i])))
    print(f'{total} values compared, {mismatches} mismatches')
"
```
