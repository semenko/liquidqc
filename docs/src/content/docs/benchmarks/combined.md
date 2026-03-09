---
title: Combined Benchmark
description: Overall performance comparison between the traditional workflow and RustQC's single-pass approach.
---

The traditional RNA-seq QC workflow requires multiple separate tools:
**featureCounts** for gene-level read counting and biotype quantification,
**dupRadar** for duplication rate analysis,
**RSeQC** for a suite of quality control metrics (including TIN for transcript integrity),
**samtools** for alignment statistics,
**preseq** for library complexity estimation, and
**Qualimap** for gene body coverage.
RustQC replaces all of these in a single pass, producing every output together.

## Traditional workflow vs. RustQC

Benchmarks were run on AWS (2026-03-07) using the nf-core/rnaseq pipeline for the
traditional tools and `rustqc rna --gtf` for RustQC.

### Large dataset: GM12878 REP1 (~186M reads)

A real-world large paired-end RNA-seq BAM aligned to GRCh38 (63,677 genes).

| Step                                    |     Traditional workflow |      RustQC |
| --------------------------------------- | -----------------------: | ----------: |
| samtools index                          |                   1m 36s |          -- |
| samtools idxstats                       |                      19s |          -- |
| samtools flagstat                       |                    2m 8s |          -- |
| samtools stats                          |                   4m 44s |          -- |
| samtools sort (for Qualimap)            |                   7m 57s |          -- |
| GTF → BED conversion                    |                      22s |          -- |
| featureCounts                           |                      56s |          -- |
| infer_experiment (RSeQC)                |                     5.7s |          -- |
| inner_distance (RSeQC)                  |                   1m 11s |          -- |
| bam_stat (RSeQC)                        |                    7m 8s |          -- |
| junction_annotation (RSeQC)             |                   6m 29s |          -- |
| junction_saturation (RSeQC)             |                   8m 19s |          -- |
| read_distribution (RSeQC)               |                   7m 35s |          -- |
| read_duplication (RSeQC)                |                  21m 53s |          -- |
| TIN (RSeQC tin.py)                      | **FAILED** (6h+ timeout) |          -- |
| Qualimap rnaseq                         |                  38m 31s |          -- |
| dupRadar                                |               1h 24m 56s |          -- |
| preseq lc_extrap                        |                  23m 26s |          -- |
| **All outputs, single pass**            |                       -- | **18m 52s** |
| **Sequential total** (excl. failed TIN) |              **~2h 58m** | **18m 52s** |

RustQC ran with 224.8% average CPU utilisation (multi-threaded), processing the
entire BAM in under 19 minutes while the sum of all individual tool runtimes
exceeds 2 hours and 58 minutes.

> **Note:** In a real pipeline, some tools can run in parallel after the BAM is
> ready, so wall-clock time for the traditional workflow will be less than the
> sequential sum. However, tools like Qualimap require a name-sorted BAM (adding
> ~8m of preprocessing), and dupRadar alone (1h 25m) sets a hard lower bound on
> wall-clock time for any parallelised run. RustQC eliminates all of this.

Key highlights from the large dataset run:

- **DupRadar** alone takes **1h 24m 56s** — more than 4× RustQC's total runtime
- **Qualimap pipeline** (sort + rnaseq) takes **46m 28s** — more than 2× RustQC's total runtime
- **read_duplication** (RSeQC) takes **21m 53s** — longer than all of RustQC
- **TIN failed entirely** in the upstream run (6+ hour timeout with no result);
  RustQC completes TIN successfully as part of its normal run
- **Peak memory:** RustQC uses 13.8 GB in a single process vs. up to 10.2 GB
  for read_duplication alone, with many other tools also requiring several GB
  each (Qualimap: 6.0 GB, GTF2BED: 2.7 GB, preseq: 2.0 GB, etc.)

### Small dataset: `test` sample (~52K reads, chr6)

A small test BAM (chromosome 6 only) used to verify correctness.

| Step                         | Traditional workflow |    RustQC |
| ---------------------------- | -------------------: | --------: |
| GTF decompression            |                479ms |        -- |
| samtools flagstat            |                  <1s |        -- |
| samtools idxstats            |                  <1s |        -- |
| samtools stats               |                  <1s |        -- |
| samtools sort (for Qualimap) |                  <1s |        -- |
| featureCounts                |                 1.0s |        -- |
| bam_stat (RSeQC)             |                 0.6s |        -- |
| infer_experiment (RSeQC)     |                 0.8s |        -- |
| read_distribution (RSeQC)    |                 1.9s |        -- |
| inner_distance (RSeQC)       |                 1.9s |        -- |
| junction_annotation (RSeQC)  |                 1.4s |        -- |
| junction_saturation (RSeQC)  |                 1.1s |        -- |
| read_duplication (RSeQC)     |                 1.3s |        -- |
| TIN (RSeQC tin.py)           |               1m 35s |        -- |
| Qualimap rnaseq              |                 5.0s |        -- |
| dupRadar                     |                 7.2s |        -- |
| preseq lc_extrap             |                 1.0s |        -- |
| **All outputs, single pass** |                   -- | **25.9s** |

On this tiny dataset the per-tool runtimes are dominated by startup overhead;
TIN is still the slowest individual tool at 1m 35s. RustQC completes all
analyses in 25.9s.

## RustQC outputs all produce identical results

Every output file produced by RustQC matches the original tools:

- **dupRadar duplication matrix** — all 14 columns match across 63,677 genes
- **Gene-level read counts** identical across all 63,677 genes
- **Model fit parameters** (intercept and slope) match to 10+ significant digits
- **Assignment statistics** (Assigned, NoFeatures, Ambiguous) match exactly
- **flagstat** all 16 metrics identical to samtools flagstat
- **idxstats** all per-chromosome counts identical to samtools idxstats
- **stats** full samtools stats format including all histogram sections (SN, FFQ, LFQ, GCF, GCL, IS, RL, etc.) matches samtools stats
- **preseq** extrapolation curve within <0.1% of preseq v3.2.0 across entire range
- **RSeQC tools** all data values and R plotting scripts match Python RSeQC.
  Junction saturation intermediate sampling points show expected stochastic
  variation from random subsampling (see [RSeQC benchmark](rseqc/))
- **Qualimap** all read counts, genomic origin percentages, and coverage bias metrics
  match Qualimap Java 2.3 exactly (see [Qualimap benchmark](qualimap/))

See the individual benchmark pages for detailed per-tool comparisons:
[dupRadar](dupradar/), [featureCounts](featurecounts/), [RSeQC](rseqc/),
[preseq](preseq/), [Samtools](samtools/), [Qualimap](qualimap/).

## Where the speedup comes from

1. **Single-pass architecture:** The traditional workflow reads the BAM file
   many times -- once per tool. RustQC reads it once and produces everything
   in that single pass.
2. **Compiled Rust vs. interpreted R/Python:** Compiled code with zero-cost
   abstractions and efficient memory management.
3. **Multi-threaded parallelism:** RustQC parallelizes across chromosomes.
   Scaling depends on the number of chromosomes with mapped reads and the
   evenness of their read distribution.

## Benchmark conditions

- **Hardware:** AWS cloud (2026-03-07), run via nf-core/rnaseq pipeline.
- **RustQC:** `--gtf` mode (all tools enabled), multi-threaded (224.8% CPU on large dataset).
- **Traditional tools:** Each tool run in its own container as part of the nf-core/rnaseq pipeline.
  Timings are wall-clock realtime from the Nextflow execution trace.
- **Reproducibility:** Benchmark scripts and data are in the
  [benchmark directory](https://github.com/seqeralabs/RustQC/tree/main/benchmark).
