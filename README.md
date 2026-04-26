# liquidqc

**Single-pass cfRNA QC + fragmentomics in Rust.**

liquidqc is a hard fork of [seqeralabs/RustQC](https://github.com/seqeralabs/RustQC)
extended with cfRNA-specific fragmentomics features. It ingests a sorted/indexed
genome BAM, a GTF annotation, and a reference FASTA, and emits one canonical
versioned JSON per sample plus a sparse per-gene Parquet — covering both
classical RNA-seq QC (RSeQC + samtools stats + Picard CollectRnaSeqMetrics
equivalents) and fragmentomics (end-motif k-mers, soft-clip k-mers,
fragment-length periodicity, per-gene Tier-2 features) in one BAM iteration.

## Status

**Phase 0 (bootstrap, 2026-04-26):** repo forked, renamed, schema v1 stub
committed, CLI scaffolds the new fragmentomics flags but does not yet compute
fragmentomics features. The inherited RustQC `rna` subcommand still runs end-
to-end and produces all classical QC outputs.

Subsequent phases (Phase 1 → 6, ~4–6 weeks total) wire RustQC accumulators
into a unified per-sample JSON writer, add Tier-1 fragmentomics (end-motif
k-mers, soft-clip k-mers, fragment-length periodicity FFT, fragment-fraction
bins), Tier-2 per-gene Parquet, Tier-3 cohort QC + sample-identity add-ons,
and a static-binary release.

## Usage

```bash
# Classical RNA-seq QC (inherited from RustQC, fully functional today)
liquidqc rna sample.markdup.bam \
  --gtf genes.gtf \
  --paired \
  --outdir results/

# Same, with the new --library-prep flag required by the v1 fragmentomics
# contract. In Phase 0 the flag is parsed but not yet acted on; from Phase 1
# onward it gates the fragmentomics output schema.
liquidqc rna sample.markdup.bam \
  --gtf genes.gtf \
  --library-prep neb_next_ultra_ii_directional \
  --paired \
  --outdir results/

# Print the schema
liquidqc schema

# Print the version (with git short hash and build timestamp)
liquidqc --version
```

`liquidqc dna` is reserved for a future cfDNA subcommand and currently exits 1.

## Inherited RustQC tooling

The classical `rna` pipeline reimplements (and is numerically validated against)
the following upstream tools, in a single BAM pass:

| Tool                | Reimplements                                                                            | Description                                                           |
| ------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| dupRadar            | [dupRadar](https://github.com/ssayols/dupRadar)                                         | PCR duplicate rate vs. expression analysis with density scatter plots |
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
| stats               | [samtools](http://www.htslib.org/) `stats`                                              | Full samtools stats output including all histogram sections           |

These outputs are format- and numerically identical to the upstream tools and
compatible with [MultiQC](https://multiqc.info/). RustQC (and therefore
liquidqc) was developed with substantial AI assistance and validated by
output comparison against the upstream tools on real sequencing data — see
the upstream
[AI & Provenance](https://seqeralabs.github.io/RustQC/about/ai-statement/)
notes.

## liquidqc additions (planned)

Per the phased roadmap:

- **Tier 1 (sample-level fragmentomics):** end-motif 4-mers (5' and 3'),
  soft-clip k-mers, fragment-length periodicity FFT, short/long fragment
  fractions and binned histogram.
- **Tier 2 (per-gene sparse Parquet):** length histogram, end motifs, mean/
  median fragment length, 5'-3' coverage bias, soft-clip rates, intron-exon
  ratio, all restricted to genes with ≥ N reads.
- **Tier 3 (cohort QC + sample identity):** per-chromosome fractions
  (chrY, chrM), hemoglobin and ribosomal-protein gene fractions, marker-
  panel reads (LM22, Tabula Sapiens, Vorperian), per-cycle quality, sex
  inference, ~17.5k-site SNP fingerprint, rRNA fraction estimate, adapter
  readthrough, per-chromosome insert-size mean, splice-site dinucleotide
  composition, gene-detection saturation curves, automatic `qc_flags`.

## Build

```bash
git clone https://github.com/semenko/liquidqc
cd liquidqc
cargo build --release
```

Requires Rust 1.87+. macOS (`brew install bzip2 xz`) and Linux
(`cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
libfontconfig1-dev pkg-config clang`) build dependencies are inherited from
RustQC.

## License

GPL-3.0-or-later (inherited from RustQC). See [LICENSE](LICENSE) and
[NOTICE](NOTICE) for the upstream attribution. Improvements that make sense
upstream are routed back to RustQC via cherry-picks; the `upstream` git
remote points there.
