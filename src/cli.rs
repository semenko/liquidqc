//! Command-line interface definition for liquidqc.
//!
//! Forked from seqeralabs/RustQC; the doc-comment summary below describes
//! the inherited `rna` subcommand. Future fragmentomics flags will be
//! added in subsequent phases.
//!
//! Provides a subcommand-based CLI. The `rna` subcommand runs all RNA-Seq QC
//! analyses in a single pass: dupRadar duplication rate analysis, featureCounts-
//! compatible output, RSeQC-equivalent metrics (bam_stat, infer_experiment,
//! read_duplication, read_distribution, junction_annotation, junction_saturation,
//! inner_distance), TIN (Transcript Integrity Number), preseq library complexity
//! extrapolation, samtools-compatible outputs (flagstat, idxstats, stats), and
//! Qualimap RNA-seq QC. Individual tools can be disabled via the YAML config file.
//!
//! A GTF gene annotation file is required for all analyses.

use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use serde::Deserialize;

/// Library strandedness protocol.
///
/// Determines how read strand is interpreted relative to the gene annotation
/// strand during counting. Accepted CLI values: `unstranded`, `forward`, `reverse`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Strandedness {
    /// Count reads on either strand (library is not strand-specific).
    #[default]
    Unstranded,
    /// Forward stranded: read 1 maps to the transcript strand.
    Forward,
    /// Reverse stranded: read 2 maps to the transcript strand (e.g. dUTP).
    Reverse,
}

impl std::fmt::Display for Strandedness {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strandedness::Unstranded => write!(f, "unstranded"),
            Strandedness::Forward => write!(f, "forward"),
            Strandedness::Reverse => write!(f, "reverse"),
        }
    }
}

/// Fast quality control tools for sequencing data, written in Rust.
#[derive(Parser, Debug)]
#[command(name = env!("CARGO_PKG_NAME"), version, about, long_about = None)]
pub struct Cli {
    /// The analysis subcommand to run.
    #[command(subcommand)]
    pub command: Commands,
}

/// Available analysis subcommands.
//
// `RnaArgs` is large; the other variants are unit structs. The size delta
// is fine here — this enum only exists once per process at startup.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// RNA-Seq QC — single-pass analysis of BAM/SAM/CRAM files.
    ///
    /// Runs featureCounts, dupRadar, Qualimap, samtools stats, and RSeQC
    /// analyses in one pass. Requires a GTF annotation and duplicate-marked
    /// (not removed) alignments.
    Rna(RnaArgs),

    /// cfDNA QC — placeholder; not implemented in v1.
    ///
    /// Reserved for a future single-pass cfDNA fragmentomics path. v1 of
    /// liquidqc is cfRNA-only; invoking this subcommand exits 1.
    Dna(DnaArgs),

    /// Download pinned GTF annotations into the per-user cache.
    ///
    /// With no arguments, fetches every supported genome. Use `--list` to
    /// see what would be fetched, or `--genome <NAME>` for a single
    /// assembly. Already-cached files are skipped unless `--force` is
    /// passed.
    FetchReferences(FetchReferencesArgs),

    /// Print the v1 JSON schema (`schema/v1/liquidqc.schema.json`) to stdout.
    Schema,

    /// Print the version (with git short hash and build timestamp).
    Version,
}

/// Arguments for the `fetch-references` subcommand.
#[derive(Parser, Debug)]
pub struct FetchReferencesArgs {
    /// Genome cache name to fetch (e.g. `GRCh38`, `GRCh37`). Omit to
    /// fetch every supported genome.
    #[arg(long, value_name = "NAME")]
    pub genome: Option<String>,

    /// Re-download even when the cache file already exists.
    #[arg(long, default_value_t = false)]
    pub force: bool,

    /// Print the table of pinned downloads and exit without fetching.
    #[arg(long, default_value_t = false)]
    pub list: bool,

    /// Suppress progress messages.
    #[arg(short = 'q', long, default_value_t = false)]
    pub quiet: bool,
}

/// Arguments for the `dna` subcommand. Reserved; not implemented in v1.
#[derive(Parser, Debug)]
pub struct DnaArgs {
    /// Input BAM (parsed but ignored — `dna` exits 1 in v1).
    #[arg(long, value_name = "BAM")]
    pub bam: Option<String>,
}

/// Arguments for the `rna` subcommand.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct RnaArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Duplicate-marked alignment file(s)
    #[arg(value_name = "INPUT", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// GTF gene annotation (plain or .gz). If omitted, liquidqc
    /// fingerprints the BAM's `@SQ` header and loads a cached GTF from
    /// `${XDG_CACHE_HOME:-~/.cache}/liquidqc/gtf/` — see
    /// `liquidqc fetch-references`.
    #[arg(
        short,
        long,
        value_name = "GTF",
        env = "RUSTQC_GTF",
        help_heading = "Input / Output"
    )]
    pub gtf: Option<String>,

    /// Reference FASTA. Required for CRAM input; required for end-motif
    /// fragmentomics (Phase 2). Other QC features run without it.
    #[arg(
        long,
        value_name = "FASTA",
        env = "LIQUIDQC_FASTA",
        help_heading = "Input / Output"
    )]
    pub fasta: Option<String>,

    /// Output directory [default: .]
    #[arg(
        short,
        long,
        default_value = ".",
        hide_default_value = true,
        env = "RUSTQC_OUTDIR",
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Override sample name for output filenames (default: derived from BAM filename)
    #[arg(
        long,
        value_name = "NAME",
        env = "RUSTQC_SAMPLE_NAME",
        help_heading = "Input / Output"
    )]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_FLAT_OUTPUT",
        help_heading = "Input / Output"
    )]
    pub flat_output: bool,

    /// YAML configuration file (see also: RUSTQC_CONFIG env var)
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    // ── Library ─────────────────────────────────────────────────────────
    /// Strandedness: unstranded, forward, reverse
    #[arg(
        short,
        long,
        value_enum,
        env = "RUSTQC_STRANDED",
        help_heading = "Library"
    )]
    pub stranded: Option<Strandedness>,

    /// Force paired-end interpretation (otherwise auto-detected).
    #[arg(
        short,
        long,
        env = "RUSTQC_PAIRED",
        conflicts_with = "single_end",
        help_heading = "Library"
    )]
    pub paired: bool,

    /// Force single-end interpretation.
    #[arg(long = "single-end", help_heading = "Library")]
    pub single_end: bool,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[arg(
        short,
        long,
        default_value_t = 1,
        hide_default_value = true,
        env = "RUSTQC_THREADS",
        help_heading = "General"
    )]
    pub threads: usize,

    /// MAPQ cutoff for quality filtering [default: 30]
    #[arg(
        short = 'Q',
        long = "mapq",
        default_value_t = 30,
        hide_default_value = true,
        env = "RUSTQC_MAPQ",
        help_heading = "General"
    )]
    pub mapq_cut: u8,

    /// GTF attribute for biotype grouping
    #[arg(
        long,
        value_name = "ATTR",
        env = "RUSTQC_BIOTYPE_ATTRIBUTE",
        help_heading = "General"
    )]
    pub biotype_attribute: Option<String>,

    /// Skip duplicate-marking check
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_DUP_CHECK",
        help_heading = "General"
    )]
    pub skip_dup_check: bool,

    /// Suppress output except warnings/errors
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        env = "RUSTQC_QUIET",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[arg(
        short = 'v',
        long,
        conflicts_with = "quiet",
        env = "RUSTQC_VERBOSE",
        help_heading = "General"
    )]
    pub verbose: bool,

    // ── Tool parameters ─────────────────────────────────────────────────
    /// infer_experiment: sample size [default: 200000]
    #[arg(
        long = "infer-experiment-sample-size",
        value_name = "N",
        env = "RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub infer_experiment_sample_size: Option<u64>,

    /// junction_annotation: min intron size [default: 50]
    #[arg(
        long = "min-intron",
        value_name = "N",
        env = "RUSTQC_MIN_INTRON",
        help_heading = "Tool parameters"
    )]
    pub min_intron: Option<u64>,

    /// junction_saturation: random seed for reproducible results
    #[arg(
        long = "junction-saturation-seed",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_SEED",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_seed: Option<u64>,

    /// junction_saturation: min coverage [default: 1]
    #[arg(
        long = "junction-saturation-min-coverage",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_MIN_COVERAGE",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_min_coverage: Option<u64>,

    /// junction_saturation: start % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-floor",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_FLOOR",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_floor: Option<u64>,

    /// junction_saturation: end % [default: 100]
    #[arg(
        long = "junction-saturation-percentile-ceiling",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_CEILING",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_ceiling: Option<u64>,

    /// junction_saturation: step % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-step",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_STEP",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_step: Option<u64>,

    /// inner_distance: sample size [default: 1000000]
    #[arg(
        long = "inner-distance-sample-size",
        value_name = "N",
        env = "RUSTQC_INNER_DISTANCE_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_sample_size: Option<u64>,

    /// inner_distance: lower bound [default: -250]
    #[arg(
        long = "inner-distance-lower-bound",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_LOWER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_lower_bound: Option<i64>,

    /// inner_distance: upper bound [default: 250]
    #[arg(
        long = "inner-distance-upper-bound",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_UPPER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_upper_bound: Option<i64>,

    /// inner_distance: bin width [default: 5]
    #[arg(
        long = "inner-distance-step",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_STEP",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_step: Option<i64>,

    /// TIN: random seed for reproducible results
    #[arg(
        long = "tin-seed",
        value_name = "N",
        env = "RUSTQC_TIN_SEED",
        help_heading = "Tool parameters"
    )]
    pub tin_seed: Option<u64>,

    /// Skip TIN analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_TIN",
        help_heading = "Tool parameters"
    )]
    pub skip_tin: bool,

    /// Skip read duplication analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_READ_DUPLICATION",
        help_heading = "Tool parameters"
    )]
    pub skip_read_duplication: bool,

    /// Skip preseq library complexity analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_PRESEQ",
        help_heading = "Tool parameters"
    )]
    pub skip_preseq: bool,

    /// preseq: random seed for bootstrap CIs
    #[arg(
        long = "preseq-seed",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEED",
        help_heading = "Tool parameters"
    )]
    pub preseq_seed: Option<u64>,

    /// preseq: max extrapolation depth
    #[arg(
        long = "preseq-max-extrap",
        value_name = "N",
        env = "RUSTQC_PRESEQ_MAX_EXTRAP",
        help_heading = "Tool parameters"
    )]
    pub preseq_max_extrap: Option<f64>,

    /// preseq: step size between points
    #[arg(
        long = "preseq-step-size",
        value_name = "N",
        env = "RUSTQC_PRESEQ_STEP_SIZE",
        help_heading = "Tool parameters"
    )]
    pub preseq_step_size: Option<f64>,

    /// preseq: bootstrap replicates for CIs
    #[arg(
        long = "preseq-n-bootstraps",
        value_name = "N",
        env = "RUSTQC_PRESEQ_N_BOOTSTRAPS",
        help_heading = "Tool parameters"
    )]
    pub preseq_n_bootstraps: Option<u32>,

    /// preseq: max segment length for PE merging
    #[arg(
        long = "preseq-seg-len",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEG_LEN",
        help_heading = "Tool parameters"
    )]
    pub preseq_seg_len: Option<i64>,

    // ── Fragmentomics (liquidqc v1) ─────────────────────────────────────
    //
    // Phase 0: parsed but not yet acted on. `--library-prep` is required
    // and never silently defaulted (per the v1 contract). Other flags
    // here become load-bearing in Phases 1–4.
    /// Library prep label (e.g. neb_next_ultra_ii_directional). Required.
    /// Never silently defaulted — pass `unknown` if truly unknown.
    #[arg(
        long = "library-prep",
        value_name = "STRING",
        env = "LIQUIDQC_LIBRARY_PREP",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub library_prep: String,

    /// Bundled marker panels to load (CSV; default: all bundled panels).
    /// Names: `lm22`, `tabula_sapiens_cfrna`, `vorperian`. Pass an empty
    /// value to disable bundled panels.
    #[arg(
        long = "panels",
        value_name = "CSV",
        env = "LIQUIDQC_PANELS",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub panels: Option<String>,

    /// Additional marker-panel TSVs to load alongside the bundled ones.
    /// Format: `gene_id<TAB>gene_symbol<TAB>panel<TAB>cell_type<TAB>weight`.
    /// Repeat the flag for multiple files.
    #[arg(
        long = "panels-tsv",
        value_name = "TSV",
        env = "LIQUIDQC_PANELS_TSV",
        help_heading = "Fragmentomics (liquidqc v1)",
        action = clap::ArgAction::Append
    )]
    pub panels_tsv: Vec<String>,

    /// Common-SNP fingerprint VCF (default: bundled somalier sites).
    /// Phase 4.
    #[arg(
        long = "snp-panel",
        value_name = "VCF",
        env = "LIQUIDQC_SNP_PANEL",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub snp_panel: Option<String>,

    /// Minimum primary read count for a gene to appear as a row in the
    /// per-gene Tier-2 Parquet sibling (default 20).
    #[arg(
        long = "min-gene-reads",
        value_name = "N",
        default_value_t = 20,
        env = "LIQUIDQC_MIN_GENE_READS",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub min_gene_reads: u32,

    /// Sample identifier for output filenames and JSON envelope.
    /// Defaults to BAM basename (without extension).
    #[arg(
        long = "sample-id",
        value_name = "ID",
        env = "LIQUIDQC_SAMPLE_ID",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub sample_id: Option<String>,

    /// Saturation curve fractions (CSV of floats in (0, 1]).
    /// Default: `0.05,0.10,0.25,0.50,0.75,1.00`.
    #[arg(
        long = "saturation-fractions",
        value_name = "CSV",
        env = "LIQUIDQC_SATURATION_FRACTIONS",
        help_heading = "Fragmentomics (liquidqc v1)"
    )]
    pub saturation_fractions: Option<String>,
}

/// Parse command-line arguments and return the Cli struct.
///
/// Sets a `long_version` that includes the git commit, build timestamp,
/// and CPU info line, shown when the user runs `liquidqc --version`.
pub fn parse_args() -> Cli {
    use clap::FromArgMatches;

    let long_version: &'static str = Box::leak(
        format!(
            "{} ({}, built {})\n{}",
            env!("CARGO_PKG_VERSION"),
            env!("GIT_SHORT_HASH"),
            env!("BUILD_TIMESTAMP"),
            crate::cpu::cpu_info_line(),
        )
        .into_boxed_str(),
    );
    let cmd = Cli::command().long_version(long_version);
    let matches = cmd.get_matches();
    Cli::from_arg_matches(&matches).expect("clap arg matching failed")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rna_default_args_gtf() {
        // Test that defaults are sensible with a GTF annotation
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert_eq!(args.gtf.as_deref(), Some("genes.gtf"));
                assert_eq!(args.stranded, None);
                assert!(!args.paired);
                assert!(!args.single_end);
                assert_eq!(args.threads, 1);
                assert_eq!(args.outdir, ".");
                assert!(args.biotype_attribute.is_none());
                assert_eq!(args.mapq_cut, 30);
                assert_eq!(args.infer_experiment_sample_size, None);
                assert_eq!(args.min_intron, None);
                assert_eq!(args.inner_distance_step, None);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_multiple_bams() {
        // Test that multiple BAM files are accepted
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "a.bam",
            "b.bam",
            "c.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["a.bam", "b.bam", "c.bam"]);
                assert_eq!(args.gtf.as_deref(), Some("genes.gtf"));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_gtf_all_args() {
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--stranded",
            "reverse",
            "--paired",
            "--threads",
            "4",
            "--outdir",
            "/tmp/out",
            "--fasta",
            "genome.fa",
            "-Q",
            "20",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.gtf.as_deref(), Some("genes.gtf"));
                assert_eq!(args.stranded, Some(Strandedness::Reverse));
                assert!(args.paired);
                assert_eq!(args.threads, 4);
                assert_eq!(args.outdir, "/tmp/out");
                assert_eq!(args.fasta, Some("genome.fa".to_string()));
                assert_eq!(args.mapq_cut, 20);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_omitting_gtf_is_allowed_at_parse_time() {
        // --gtf is now optional; resolution to a cache file (or a clear
        // missing-cache error) happens at run time. Parse must succeed
        // here; tests/phase0_cli covers the runtime behaviour.
        let result = Cli::try_parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--library-prep",
            "unknown",
        ]);
        assert!(result.is_ok(), "Parse should succeed without --gtf");
    }

    #[test]
    fn test_rna_paired_and_single_end_conflict() {
        let result = Cli::try_parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--paired",
            "--single-end",
        ]);
        assert!(
            result.is_err(),
            "Expected error when both --paired and --single-end are passed"
        );
    }

    #[test]
    fn test_rna_missing_library_prep() {
        // --library-prep is required by the v1 contract; never silently defaulted
        let result = Cli::try_parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
        ]);
        assert!(
            result.is_err(),
            "Expected error when --library-prep is not provided"
        );
    }

    #[test]
    fn test_rna_rseqc_params() {
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--infer-experiment-sample-size",
            "500000",
            "--min-intron",
            "100",
            "--junction-saturation-min-coverage",
            "5",
            "--inner-distance-lower-bound",
            "-500",
            "--inner-distance-upper-bound",
            "500",
            "--inner-distance-step",
            "10",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.infer_experiment_sample_size, Some(500_000));
                assert_eq!(args.min_intron, Some(100));
                assert_eq!(args.junction_saturation_min_coverage, Some(5));
                assert_eq!(args.inner_distance_lower_bound, Some(-500));
                assert_eq!(args.inner_distance_upper_bound, Some(500));
                assert_eq!(args.inner_distance_step, Some(10));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_preseq_params() {
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--preseq-max-extrap",
            "5000000000",
            "--preseq-step-size",
            "500000",
            "--preseq-n-bootstraps",
            "200",
            "--preseq-seg-len",
            "100000000",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(!args.skip_preseq);
                assert_eq!(args.preseq_max_extrap, Some(5_000_000_000.0));
                assert_eq!(args.preseq_step_size, Some(500_000.0));
                assert_eq!(args.preseq_n_bootstraps, Some(200));
                assert_eq!(args.preseq_seg_len, Some(100_000_000));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_tool_seeds() {
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--preseq-seed",
            "1",
            "--tin-seed",
            "2",
            "--junction-saturation-seed",
            "3",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.preseq_seed, Some(1));
                assert_eq!(args.tin_seed, Some(2));
                assert_eq!(args.junction_saturation_seed, Some(3));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_skip_preseq() {
        let cli = Cli::parse_from([
            env!("CARGO_PKG_NAME"),
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--library-prep",
            "unknown",
            "--skip-preseq",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(args.skip_preseq);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }
}
