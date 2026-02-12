# dupRust Changelog  🧬 🦀

## [Unreleased]

Initial release of dupRust, a fast Rust reimplementation of dupRadar.

- Single-pass BAM counting engine with CIGAR-aware alignment block extraction
- Support for single-end and paired-end libraries
- Strand-aware counting (unstranded, forward, reverse-stranded)
- Multimapper and duplicate tracking across four simultaneous counting modes
- GTF annotation parser with non-overlapping exon length calculation
- 14-column duplication matrix output matching R dupRadar format exactly
- Logistic regression fitting via IRLS (matches R's `glm()` defaults)
- Density scatter plot matching R's `densCols(nbin=500)` behavior
- Boxplot of duplication rate by expression quantile
- Expression histogram with log10(RPK) distribution
- MultiQC-compatible output files (intercept stats + fitted curve)
- YAML configuration file support for chromosome name mapping
- Cross-platform release builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container published to `ghcr.io/ewels/duprust`
- CI pipeline with tests, formatting checks, and clippy linting
- Benchmark suite comparing against R dupRadar on GM12878 dataset (~7x faster)
