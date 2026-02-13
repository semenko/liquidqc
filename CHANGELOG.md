# dupRust Changelog 🧬 🦀

## [Version 0.1.0](https://github.com/ewels/dupRust/releases/tag/v0.1.0) - 2026-02-13

Initial release of dupRust — a fast Rust reimplementation of [dupRadar](https://github.com/ssayols/dupRadar) for assessing PCR duplicate rates in RNA-Seq datasets.

### Features

- Drop-in replacement for R dupRadar with identical numerical output
- Single-end and paired-end library support
- Strand-aware counting (unstranded, forward, reverse-stranded)
- Multi-threaded BAM processing across chromosomes (`--threads`)
- 14-column duplication matrix, density scatter plot, boxplot, and expression histogram
- MultiQC-compatible output files
- YAML configuration for chromosome name mapping
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container at `ghcr.io/ewels/duprust`

### Performance

- ~7x faster than R dupRadar single-threaded, ~27x faster with 10 threads
- Parallel BAM processing using rayon thread pool
- Cache-oblivious interval trees (coitrees) for fast overlap queries
- Interned gene IDs and reusable buffers to minimise allocations
