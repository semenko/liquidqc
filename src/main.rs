//! liquidqc — single-pass cfRNA QC + fragmentomics
//!
//! Forked from seqeralabs/RustQC (GPLv3+). The inherited `rna` subcommand
//! still runs the full RNA-Seq QC analysis suite in a single BAM pass:
//! dupRadar duplication rate analysis, featureCounts-compatible gene counting,
//! 8 RSeQC-equivalent tools (bam_stat, infer_experiment, read_duplication,
//! read_distribution, junction_annotation, junction_saturation, inner_distance, TIN),
//! preseq library complexity extrapolation, samtools-compatible outputs
//! (flagstat, idxstats, stats), and Qualimap gene body coverage profiling.
//!
//! liquidqc extends this base with cfRNA-specific fragmentomics features
//! (end-motif k-mers, soft-clip k-mers, fragment-length periodicity, per-gene
//! Tier-2 features, panels, sex inference, SNP fingerprinting, saturation
//! curves) and a versioned JSON+Parquet output schema (`schema/v1/`).
//! Individual tools can be disabled via the YAML config file.

mod citations;
mod cli;
mod config;
mod cpu;
mod envelope;
mod gtf;
mod hash;
mod io;
mod qc_flags;
mod rna;
mod runtime_stats;
mod ui;

use anyhow::{ensure, Context, Result};
use indexmap::IndexMap;
use log::debug;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

use ui::{format_count, format_duration, format_pct, Ui, Verbosity};

use rust_htslib::bam::Read as BamRead;

use rna::rseqc::accumulators::{RseqcAccumulators, RseqcAnnotations, RseqcConfig};

fn main() -> Result<()> {
    // Guard against running a SIMD-optimized binary on incompatible hardware.
    // Must run before any auto-vectorized code to prevent SIGILL.
    cpu::check_cpu_compat()?;

    let cli = cli::parse_args();

    // Subcommands that don't run the full RNA pipeline are handled before
    // env_logger init so they have minimal output.
    match &cli.command {
        cli::Commands::Schema => {
            const SCHEMA: &str = include_str!("../schema/v1/liquidqc.schema.json");
            print!("{SCHEMA}");
            return Ok(());
        }
        cli::Commands::Version => {
            println!(
                "liquidqc {} ({}, built {})",
                env!("CARGO_PKG_VERSION"),
                env!("GIT_SHORT_HASH"),
                env!("BUILD_TIMESTAMP"),
            );
            return Ok(());
        }
        cli::Commands::Dna(_) => {
            eprintln!("liquidqc dna: not implemented in v1; planned");
            std::process::exit(1);
        }
        cli::Commands::Rna(_) => {}
    }

    // Determine verbosity from CLI flags
    let verbosity = match &cli.command {
        cli::Commands::Rna(args) if args.quiet => Verbosity::Quiet,
        cli::Commands::Rna(args) if args.verbose => Verbosity::Verbose,
        _ => Verbosity::Normal,
    };

    // Initialize env_logger: only for debug/trace (user-facing output goes through Ui).
    // In verbose mode, lower the threshold so debug!() messages are visible too.
    let log_level = match verbosity {
        Verbosity::Quiet => "warn",
        Verbosity::Normal => "warn",
        Verbosity::Verbose => "debug",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level))
        .format_timestamp(None)
        .init();

    let ui = Ui::new(verbosity);

    match cli.command {
        cli::Commands::Rna(args) => run_rna(args, &ui),
        // Schema/Version/Dna are handled before env_logger init above and
        // early-return out of main(); reaching them here would be a bug.
        cli::Commands::Schema | cli::Commands::Version | cli::Commands::Dna(_) => {
            unreachable!("non-rna subcommands handled before this match")
        }
    }
}

/// Reconstruct the command line for the featureCounts-compatible header comment.
fn reconstruct_command_line(args: &cli::RnaArgs) -> String {
    let mut parts = vec![format!(
        "liquidqc rna {}",
        args.input
            .iter()
            .map(|s| shell_escape(s))
            .collect::<Vec<_>>()
            .join(" "),
    )];
    parts.push(format!("--gtf {}", shell_escape(&args.gtf)));
    if let Some(s) = args.stranded {
        parts.push(format!("-s {}", s));
    }
    if args.paired {
        parts.push("-p".to_string());
    }
    if args.threads != 1 {
        parts.push(format!("-t {}", args.threads));
    }
    if args.outdir != "." {
        parts.push(format!("-o {}", shell_escape(&args.outdir)));
    }
    if let Some(ref config_path) = args.config {
        parts.push(format!("-c {}", shell_escape(config_path)));
    }
    if let Some(ref biotype) = args.biotype_attribute {
        parts.push(format!("--biotype-attribute {}", shell_escape(biotype)));
    }
    if let Some(ref fasta) = args.fasta {
        parts.push(format!("--fasta {}", shell_escape(fasta)));
    }
    if args.flat_output {
        parts.push("--flat-output".to_string());
    }
    if args.skip_dup_check {
        parts.push("--skip-dup-check".to_string());
    }
    parts.join(" ")
}

/// Shell-escaping: wrap in single quotes if the string contains shell metacharacters.
///
/// Single quotes prevent all shell interpretation. Any embedded single quotes
/// are escaped using the `'\''` pattern (end quote, escaped quote, restart quote).
fn shell_escape(s: &str) -> String {
    if s.contains(|c: char| c.is_whitespace() || "\"'\\$`!#&|;(){}[]<>?*~".contains(c)) {
        format!("'{}'", s.replace('\'', "'\\''"))
    } else {
        s.to_string()
    }
}

/// Run the full RNA-Seq QC pipeline: dupRadar + featureCounts + RSeQC analyses.
///
/// When multiple BAM files are provided, they are processed in parallel using
/// rayon. The GTF annotation is parsed once and shared across all BAM files.
/// Available threads are distributed across the parallel BAM processing jobs.
fn run_rna(args: cli::RnaArgs, ui: &Ui) -> Result<()> {
    // Validate thread count
    ensure!(
        args.threads >= 1,
        "--threads must be at least 1 (got {})",
        args.threads
    );

    // Load and merge configuration files from all sources:
    // XDG system → XDG user → RUSTQC_CONFIG env → explicit -c flag.
    // Each layer overrides only the leaf fields it sets.
    let (cfg, config_sources) = config::load_merged_config(args.config.as_deref())?;
    let mut config = cfg.rna;

    if !config_sources.is_empty() {
        for (path, source) in &config_sources {
            ui.detail(&format!("Loaded config: {} ({})", path.display(), source));
        }
        if config.has_chromosome_mapping() {
            ui.detail(&format!(
                "Chromosome name mapping: {} entries",
                config.chromosome_mapping.len()
            ));
        }
    }

    // Apply CLI overrides to skip flags
    if args.skip_tin {
        config.tin.enabled = false;
    }
    if args.skip_read_duplication {
        config.read_duplication.enabled = false;
    }
    if args.skip_preseq {
        config.preseq.enabled = false;
    }
    if let Some(val) = args.preseq_max_extrap {
        config.preseq.max_extrap = val;
    }
    if let Some(val) = args.preseq_step_size {
        config.preseq.step_size = val;
    }
    if let Some(val) = args.preseq_n_bootstraps {
        config.preseq.n_bootstraps = val;
    }
    if let Some(val) = args.preseq_seg_len {
        config.preseq.max_segment_length = val;
    }

    // Apply CLI overrides for RSeQC tool parameters
    if let Some(val) = args.infer_experiment_sample_size {
        config.infer_experiment.sample_size = Some(val);
    }
    if let Some(val) = args.min_intron {
        config.junction_annotation.min_intron = Some(val);
    }
    if let Some(val) = args.junction_saturation_min_coverage {
        config.junction_saturation.min_coverage = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_floor {
        config.junction_saturation.percentile_floor = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_ceiling {
        config.junction_saturation.percentile_ceiling = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_step {
        config.junction_saturation.percentile_step = Some(val);
    }
    if let Some(val) = args.inner_distance_sample_size {
        config.inner_distance.sample_size = Some(val);
    }
    if let Some(val) = args.inner_distance_lower_bound {
        config.inner_distance.lower_bound = Some(val);
    }
    if let Some(val) = args.inner_distance_upper_bound {
        config.inner_distance.upper_bound = Some(val);
    }
    if let Some(val) = args.inner_distance_step {
        config.inner_distance.step = Some(val);
    }

    // Apply per-tool seed overrides from CLI flags
    if let Some(seed) = args.preseq_seed {
        config.preseq.seed = seed;
    }
    if let Some(seed) = args.tin_seed {
        config.tin.seed = Some(seed);
    }
    if let Some(seed) = args.junction_saturation_seed {
        config.junction_saturation.seed = Some(seed);
    }

    // Warn early if CRAM input is likely but no FASTA is provided
    if args.input.iter().any(|f| f.ends_with(".cram"))
        && args.fasta.is_none()
        && std::env::var("REF_PATH").is_err()
        && std::env::var("REF_CACHE").is_err()
    {
        anyhow::bail!(
            "CRAM input requires a reference FASTA. \
             Pass --fasta <FASTA> or set the REF_PATH environment variable."
        );
    }

    // Validate all input alignment files before expensive GTF parsing
    for bam_path in &args.input {
        let mut reader = rust_htslib::bam::Reader::from_path(bam_path)
            .with_context(|| format!("Cannot open alignment file '{}'", bam_path))?;
        if let Some(ref fasta) = args.fasta {
            reader
                .set_reference(fasta)
                .with_context(|| format!("Cannot set reference for CRAM file '{}'", bam_path))?;
        }
        let _header = reader.header().clone();
        // Reader dropped — just validating the file is openable
    }

    let start = Instant::now();
    let n_bams = args.input.len();

    // Hash inputs for the v1 envelope provenance fields. GTF + FASTA hash
    // once and are shared across all BAMs; per-BAM MD5s are computed in
    // parallel before per-BAM processing starts.
    let gtf_md5 =
        hash::md5_uncompressed(&args.gtf).with_context(|| format!("hashing GTF: {}", args.gtf))?;
    let gtf_path_abs = absolutize(&args.gtf);
    let reference_fasta_md5 = match args.fasta {
        Some(ref p) => {
            Some(hash::md5(p).with_context(|| format!("hashing reference FASTA: {}", p))?)
        }
        None => None,
    };
    let bam_hashes: HashMap<String, (String, u64)> = args
        .input
        .par_iter()
        .map(|p| {
            hash::md5_and_size(p)
                .map(|h| (p.clone(), h))
                .with_context(|| format!("hashing BAM: {}", p))
        })
        .collect::<Result<_>>()?;

    // --sample-name / --sample-id / config sample_name only make sense for a
    // single BAM — using them with multi-BAM would collide output filenames.
    let effective_sample_name = args
        .sample_name
        .as_deref()
        .or(config.sample_name.as_deref());
    if n_bams > 1 && effective_sample_name.is_some() {
        let source = if args.sample_name.is_some() {
            "--sample-name flag"
        } else {
            "config file sample_name"
        };
        anyhow::bail!(
            "{source} cannot be used with multiple input files \
             (would produce identical output filenames)"
        );
    }
    if n_bams > 1 && args.sample_id.is_some() {
        anyhow::bail!(
            "--sample-id cannot be used with multiple input files \
             (would produce identical envelope filenames)"
        );
    }

    let cpu_info = cpu::cpu_info_line();
    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&cpu_info),
    );

    // Resolve effective stranded/paired early for display (config loaded above)
    let effective_stranded = args
        .stranded
        .or(config.stranded)
        .unwrap_or(cli::Strandedness::Unstranded);
    let effective_paired = args.paired || config.paired.unwrap_or(false);

    if n_bams == 1 {
        ui.config("Input", &args.input[0]);
    } else {
        ui.config("Input", &format!("{} BAM files", n_bams));
        for f in &args.input {
            ui.detail(&format!("  {f}"));
        }
    }
    ui.config("Annotation", &args.gtf);
    if let Some(ref fasta) = args.fasta {
        ui.config("Reference", fasta);
    }
    ui.config("Stranded", &effective_stranded.to_string());
    ui.config("Paired", &effective_paired.to_string());
    ui.config("CPU Threads", &args.threads.to_string());
    let outdir_display = if std::path::Path::new(&args.outdir).is_relative() {
        format!("./{}", args.outdir)
    } else {
        args.outdir.clone()
    };
    ui.config("Output dir", &outdir_display);

    // Determine biotype attribute name (CLI overrides config, with auto-detection fallback)
    let configured_biotype = args
        .biotype_attribute
        .clone()
        .unwrap_or_else(|| config.featurecounts.biotype_attribute.clone());

    // Determine which extra GTF attributes we need, and parse GTF if provided
    let mut extra_attributes: Vec<String> = Vec::new();
    let mut biotype_attribute = configured_biotype.clone();
    let need_biotype = config.any_biotype_output();

    // Detect biotype attributes in GTF
    let gtf_path = &args.gtf;
    if need_biotype {
        // Check if the configured biotype attribute exists in the GTF.
        // If not found and the user didn't explicitly set it, try common alternatives:
        // Ensembl GTFs use "gene_biotype", GENCODE GTFs use "gene_type".
        let user_explicit = args.biotype_attribute.is_some();
        if gtf::attribute_exists_in_gtf(gtf_path, &biotype_attribute, 1000) {
            extra_attributes.push(biotype_attribute.clone());
            ui.detail(&format!("Biotype attribute: {}", biotype_attribute));
        } else if !user_explicit {
            // Auto-detect: try known alternatives
            let alternatives = if biotype_attribute == "gene_biotype" {
                vec!["gene_type"]
            } else if biotype_attribute == "gene_type" {
                vec!["gene_biotype"]
            } else {
                vec!["gene_biotype", "gene_type"]
            };
            let mut found = false;
            for alt in &alternatives {
                if gtf::attribute_exists_in_gtf(gtf_path, alt, 1000) {
                    ui.detail(&format!(
                        "Biotype attribute '{}' not found, using '{}'",
                        biotype_attribute, alt
                    ));
                    biotype_attribute = alt.to_string();
                    extra_attributes.push(biotype_attribute.clone());
                    found = true;
                    break;
                }
            }
            if !found {
                let tried: Vec<_> = std::iter::once(configured_biotype.as_str())
                    .chain(alternatives.iter().copied())
                    .collect();
                let names = tried
                    .iter()
                    .map(|a| format!("'{a}'"))
                    .collect::<Vec<_>>()
                    .join(" and ");
                ui.warn(&format!(
                    "Biotype attributes {} not found in GTF, skipping biotype outputs \
                         (use --biotype-attribute to specify)",
                    names
                ));
            }
        } else {
            ui.warn(&format!(
                "Biotype attribute '{}' not found in GTF, skipping biotype outputs",
                biotype_attribute
            ));
        }
    }

    // Step 1: Parse GTF annotation (shared across all BAM files)
    ui.blank();
    ui.step("Parsing GTF annotation...");
    let gtf_start = Instant::now();
    // Phase 4 gene-class fractions resolve genes by `gene_name`; ensure the
    // attribute is parsed even when the user didn't ask for biotype outputs.
    if !extra_attributes.iter().any(|a| a == "gene_name") {
        extra_attributes.push("gene_name".to_string());
    }
    let genes = gtf::parse_gtf(gtf_path, &extra_attributes)?;
    ui.detail(&format!(
        "Parsed {} genes in {}",
        format_count(genes.len() as u64),
        format_duration(gtf_start.elapsed()),
    ));

    // Build RSeQC data structures from GTF annotation.
    // These are built once and shared across all BAM files.
    // Each tool's data is only built when enabled in the config.
    let gene_model = if config.infer_experiment.enabled {
        ui.detail("Building gene model for infer_experiment...");
        Some(rna::rseqc::infer_experiment::GeneModel::from_genes(&genes))
    } else {
        None
    };

    let ref_junctions = if config.junction_annotation.enabled {
        ui.detail("Building reference junctions...");
        Some(rna::rseqc::common::build_reference_junctions_from_genes(
            &genes,
        ))
    } else {
        None
    };

    let known_junctions = if config.junction_saturation.enabled {
        ui.detail("Building known junction set...");
        Some(rna::rseqc::common::build_known_junctions_from_genes(&genes))
    } else {
        None
    };

    let rd_regions = if config.read_distribution.enabled {
        ui.detail("Building genomic region sets...");
        Some(rna::rseqc::read_distribution::build_regions_from_genes(
            &genes,
        ))
    } else {
        None
    };

    let exon_bitset = if config.inner_distance.enabled {
        ui.detail("Building exon bitset...");
        Some(rna::rseqc::inner_distance::ExonBitset::from_genes(&genes))
    } else {
        None
    };

    let transcript_tree = if config.inner_distance.enabled {
        ui.detail("Building transcript tree...");
        Some(rna::rseqc::inner_distance::TranscriptTree::from_genes(
            &genes,
        ))
    } else {
        None
    };

    let tin_sample_size = config.tin.sample_size.unwrap_or(100) as usize;
    let tin_index = if config.tin.enabled {
        ui.detail("Building TIN index...");
        Some(rna::rseqc::tin::TinIndex::from_genes(
            &genes,
            tin_sample_size,
        ))
    } else {
        None
    };

    let chrom_mapping = config.alignment_to_gtf_mapping();
    let chrom_prefix = config.chromosome_prefix().map(|s| s.to_owned());

    // Reconstruct command line for featureCounts-compatible header
    let command_line = reconstruct_command_line(&args);

    // Validate that input file stems are unique (otherwise outputs would collide).
    if n_bams > 1 {
        let mut seen_stems = HashSet::new();
        for bam_path in &args.input {
            let stem = Path::new(bam_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(bam_path);
            anyhow::ensure!(
                seen_stems.insert(stem.to_owned()),
                "Duplicate BAM file stem '{}': multiple BAM files with the same \
                 filename would produce conflicting output files. Rename or \
                 reorganise input files so each has a unique filename.",
                stem
            );
        }
    }

    // Determine thread allocation for parallel BAM processing.
    // When processing multiple BAMs, we run BAMs in parallel and divide threads
    // among them. Each BAM's count_reads() creates its own rayon pool internally.
    // Note: the outer pool threads are mostly blocked waiting on inner pools, so
    // actual CPU-active threads stay close to `args.threads`. However, the total
    // OS thread count may briefly exceed `--threads` due to the outer pool threads
    // and temporary plot-generation threads (3 per BAM via std::thread::scope).
    // Integer division may leave up to `n_parallel - 1` threads unused.
    let n_parallel = n_bams.min(args.threads).max(1);
    let threads_per_bam = (args.threads / n_parallel).max(1);
    if n_bams > 1 {
        ui.detail(&format!(
            "Processing {} BAM files ({} in parallel, {} threads each)",
            n_bams, n_parallel, threads_per_bam
        ));
    }

    // Create output directory
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    // Determine if biotype attribute was found in the GTF
    let biotype_in_gtf = extra_attributes.contains(&biotype_attribute);

    // Effective flat_output: true if enabled by either CLI flag or config file
    let flat_output = args.flat_output || config.flat_output;

    // Build the shared parameters struct for process_single_bam
    let shared = SharedParams {
        ui,
        stranded: effective_stranded,
        paired: effective_paired,
        chrom_mapping: &chrom_mapping,
        chrom_prefix: chrom_prefix.as_deref(),
        outdir,
        flat_output,
        fasta: args.fasta.as_deref(),
        skip_dup_check: args.skip_dup_check,
        config: &config,
        biotype_attribute: &biotype_attribute,
        biotype_in_gtf,
        command_line: &command_line,
        gene_model: gene_model.as_ref(),
        ref_junctions: ref_junctions.as_ref(),
        known_junctions: known_junctions.as_ref(),
        rd_regions: rd_regions.as_ref(),
        exon_bitset: exon_bitset.as_ref(),
        transcript_tree: transcript_tree.as_ref(),
        mapq_cut: args.mapq_cut,
        infer_experiment_sample_size: config.infer_experiment.sample_size.unwrap_or(200_000),
        min_intron: config.junction_annotation.min_intron.unwrap_or(50),
        junction_saturation_min_coverage: config.junction_saturation.min_coverage.unwrap_or(1),
        junction_saturation_percentile_floor: config
            .junction_saturation
            .percentile_floor
            .unwrap_or(5),
        junction_saturation_percentile_ceiling: config
            .junction_saturation
            .percentile_ceiling
            .unwrap_or(100),
        junction_saturation_percentile_step: config
            .junction_saturation
            .percentile_step
            .unwrap_or(5),
        junction_saturation_seed: config.junction_saturation.seed.unwrap_or(42),
        inner_distance_sample_size: config.inner_distance.sample_size.unwrap_or(1_000_000),
        inner_distance_lower_bound: config.inner_distance.lower_bound.unwrap_or(-250),
        inner_distance_upper_bound: config.inner_distance.upper_bound.unwrap_or(250),
        inner_distance_step: config.inner_distance.step.unwrap_or(5),
        tin_index: tin_index.as_ref(),
        tin_sample_size,
        tin_min_coverage: config.tin.min_coverage.unwrap_or(10),
        gtf_path: &args.gtf,
        gtf_path_abs: &gtf_path_abs,
        gtf_md5: &gtf_md5,
        reference_fasta_md5: reference_fasta_md5.as_deref(),
        bam_hashes: &bam_hashes,
        library_prep: &args.library_prep,
        sample_id_override: args.sample_id.as_deref(),
        sample_name_override: effective_sample_name,
        min_gene_reads: args.min_gene_reads,
        panels: args.panels.as_deref(),
        panels_tsv: args.panels_tsv.clone(),
        saturation_fractions: args.saturation_fractions.as_deref(),
        snp_panel: args.snp_panel.as_deref(),
    };

    // Step 2: Process all alignment files (in parallel when multiple)
    let bam_results: Vec<(String, Result<BamResult>)> = if n_bams == 1 {
        // Single file: use all threads directly, no outer rayon pool needed
        vec![(
            args.input[0].clone(),
            process_single_bam(&args.input[0], &genes, args.threads, &shared),
        )]
    } else {
        // Multiple files: process in parallel with a dedicated rayon pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_parallel)
            .build()
            .context("Failed to create rayon thread pool for parallel BAM processing")?;

        pool.install(|| {
            args.input
                .par_iter()
                .map(|bam_path| {
                    (
                        bam_path.clone(),
                        process_single_bam(bam_path, &genes, threads_per_bam, &shared),
                    )
                })
                .collect()
        })
    };

    // Collect results
    let mut n_err = 0;
    for (bam_path, result) in &bam_results {
        if let Err(e) = result {
            n_err += 1;
            ui.error(&format!("Failed to process {}: {:#}", bam_path, e));
        }
    }

    // Multi-BAM summary
    if n_bams > 1 {
        ui.blank();
        ui.section(&format!(
            "Processed {} files in {}",
            n_bams,
            format_duration(start.elapsed()),
        ));
        for (bam_path, result) in &bam_results {
            let name = Path::new(bam_path)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or(bam_path);
            match result {
                Ok(bam_result) => {
                    ui.bam_result_ok(name, bam_result.duration);
                }
                Err(e) => {
                    let msg = format!("{:#}", e);
                    // Truncate long error messages for the summary line
                    let short = console::truncate_str(&msg, 60, "…");
                    ui.bam_result_err(name, &short);
                }
            }
        }
    }

    let elapsed = start.elapsed();

    // Write citations file
    let citations_path = outdir.join("CITATIONS.md");
    citations::write_citations(
        &citations_path,
        &config,
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
    )?;
    ui.output_item("Citations", &citations_path.display().to_string());

    // Check for strandedness mismatch between user-specified and inferred values
    for (bam_path, result) in &bam_results {
        if let Ok(bam_result) = result {
            if let Some(ref ie_result) = bam_result.infer_experiment {
                if let Some((inferred, suggestion)) =
                    rna::rseqc::infer_experiment::check_strandedness_mismatch(
                        ie_result,
                        effective_stranded,
                    )
                {
                    let bam_name = Path::new(bam_path)
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or(bam_path);
                    let line1 = "Strandedness mismatch detected!".to_string();
                    let line2 = format!(
                        "{} - you specified '--stranded {}' but infer_experiment suggests '{}'",
                        bam_name, effective_stranded, inferred,
                    );
                    let line3 = format!(
                        "(forward fraction: {:.4}, reverse fraction: {:.4})",
                        ie_result.frac_protocol1, ie_result.frac_protocol2,
                    );
                    let line4 = format!("Consider re-running with '--stranded {}'", suggestion,);
                    ui.blank();
                    ui.warn_box(&[&line1, &line2, &line3, &line4]);
                }
            }
        }
    }

    ui.finish("liquidqc run finished", elapsed);

    if n_err > 0 {
        anyhow::bail!("{} file(s) failed to process", n_err);
    }

    Ok(())
}

// ============================================================================
// Shared parameters
// ============================================================================

/// Parameters shared across all BAM files in a run.
///
/// Bundles the read-only configuration, annotation data, and tool parameters
/// that are computed once in `run_rna()` and passed to each `process_single_bam()`.
/// This avoids a long parameter list on the processing function.
struct SharedParams<'a> {
    /// Terminal UI handle.
    ui: &'a Ui,
    /// Library strandedness.
    stranded: cli::Strandedness,
    /// Whether the library is paired-end.
    paired: bool,
    /// Alignment-to-GTF chromosome name mapping.
    chrom_mapping: &'a HashMap<String, String>,
    /// Optional chromosome name prefix.
    chrom_prefix: Option<&'a str>,
    /// Output directory for results.
    outdir: &'a Path,
    /// When true, write all files directly to outdir (no subdirectories).
    flat_output: bool,
    /// Optional override for the sample name used in output filenames.
    sample_name_override: Option<&'a str>,
    /// Optional reference FASTA. Used for CRAM decoding and end-motif
    /// fragmentomics (Phase 2).
    fasta: Option<&'a str>,
    /// Whether to skip duplicate-marking validation.
    skip_dup_check: bool,
    /// Configuration for conditional outputs.
    config: &'a config::RnaConfig,
    /// GTF attribute name for biotype counting.
    biotype_attribute: &'a str,
    /// Whether the biotype attribute was found in the GTF.
    biotype_in_gtf: bool,
    /// Reconstructed command line for featureCounts header.
    command_line: &'a str,
    /// Pre-built gene model for infer_experiment (from GTF).
    gene_model: Option<&'a rna::rseqc::infer_experiment::GeneModel>,
    /// Pre-built reference junctions for junction_annotation (from GTF).
    ref_junctions: Option<&'a rna::rseqc::common::ReferenceJunctions>,
    /// Pre-built known junction set for junction_saturation (from GTF).
    known_junctions: Option<&'a rna::rseqc::common::KnownJunctionSet>,
    /// Pre-built genomic region sets for read_distribution (from GTF).
    rd_regions: Option<&'a rna::rseqc::read_distribution::RegionSets>,
    /// Pre-built exon bitset for inner_distance (from GTF).
    exon_bitset: Option<&'a rna::rseqc::inner_distance::ExonBitset>,
    /// Pre-built transcript tree for inner_distance (from GTF).
    transcript_tree: Option<&'a rna::rseqc::inner_distance::TranscriptTree>,
    /// MAPQ cutoff for read quality filtering.
    mapq_cut: u8,
    /// Maximum reads to sample for strandedness inference.
    infer_experiment_sample_size: u64,
    /// Minimum intron size for junction filtering.
    min_intron: u64,
    /// Minimum coverage for junction saturation.
    junction_saturation_min_coverage: u64,
    /// Sampling start percentage for junction saturation.
    junction_saturation_percentile_floor: u64,
    /// Sampling end percentage for junction saturation.
    junction_saturation_percentile_ceiling: u64,
    /// Sampling step percentage for junction saturation.
    junction_saturation_percentile_step: u64,
    /// Random seed for junction saturation shuffle.
    junction_saturation_seed: u64,
    /// Maximum read pairs to sample for inner distance.
    inner_distance_sample_size: u64,
    /// Lower bound of inner distance histogram.
    inner_distance_lower_bound: i64,
    /// Upper bound of inner distance histogram.
    inner_distance_upper_bound: i64,
    /// Bin width for inner distance histogram.
    inner_distance_step: i64,
    /// Pre-built TIN index for transcript integrity analysis (from GTF).
    tin_index: Option<&'a rna::rseqc::tin::TinIndex>,
    /// Number of equally-spaced positions to sample per transcript for TIN.
    tin_sample_size: usize,
    /// Minimum read-start count per transcript to compute TIN.
    tin_min_coverage: u32,
    /// Path to GTF file (for Qualimap report output).
    gtf_path: &'a str,
    /// Absolute GTF path string for the envelope (computed once in `run_rna`).
    gtf_path_abs: &'a str,
    /// MD5 of the (uncompressed) GTF bytes — schema-required envelope field.
    gtf_md5: &'a str,
    /// MD5 of the reference FASTA file as supplied, or `None` if no
    /// `--fasta` was provided.
    reference_fasta_md5: Option<&'a str>,
    /// Per-BAM (md5, size_bytes) keyed by the path string passed on the CLI.
    bam_hashes: &'a HashMap<String, (String, u64)>,
    /// `--library-prep` value (required, never silently defaulted).
    library_prep: &'a str,
    /// Optional `--sample-id` override (single-BAM only); otherwise
    /// `sample_name_override` and BAM stem are used.
    sample_id_override: Option<&'a str>,
    /// Minimum primary read count for a gene to appear as a row in the
    /// per-gene Tier-2 Parquet (`--min-gene-reads`).
    min_gene_reads: u32,
    /// CSV of bundled panel names selected via `--panels`. `None` =
    /// all bundled panels load (default). Empty string disables bundled.
    panels: Option<&'a str>,
    /// User-supplied panel TSV paths via `--panels-tsv`.
    panels_tsv: Vec<String>,
    /// CSV of saturation-curve fractions in (0, 1]. `None` = use default
    /// `[0.05, 0.10, 0.25, 0.50, 0.75, 1.00]`.
    saturation_fractions: Option<&'a str>,
    /// User SNP-site panel (TSV). `None` loads the bundled smoke panel.
    snp_panel: Option<&'a str>,
}

// ============================================================================
// Per-BAM result
// ============================================================================

/// Per-BAM result returned from `process_single_bam`.
///
/// Carries only the fields used by the run-level loop (timing for the
/// summary box, infer_experiment result for the strandedness-mismatch
/// warning). The canonical per-sample output is the `<sample_id>.liquidqc.json`
/// envelope written from inside `process_single_bam`.
#[derive(Debug)]
struct BamResult {
    /// Wall-clock processing time.
    duration: std::time::Duration,
    /// infer_experiment result (if the tool was enabled).
    infer_experiment: Option<rna::rseqc::infer_experiment::InferExperimentResult>,
}

// ============================================================================
// Per-BAM processing
// ============================================================================

/// Process a single alignment file through the full analysis pipeline.
///
/// This runs the complete analysis for one file: counting, featureCounts output,
/// biotype counting, duplication matrix, model fitting, plotting, MultiQC output,
/// and all enabled RSeQC analyses.
///
/// # Arguments
///
/// * `bam_path` - Path to the duplicate-marked alignment file (SAM/BAM/CRAM)
/// * `genes` - Parsed GTF gene annotations
/// * `threads` - Number of threads for this file's read counting
/// * `params` - Shared parameters (config, annotations, tool settings)
fn process_single_bam(
    bam_path: &str,
    genes: &IndexMap<String, gtf::Gene>,
    threads: usize,
    params: &SharedParams,
) -> Result<BamResult> {
    let ui = params.ui;
    let bam_stem = Path::new(bam_path)
        .file_stem()
        .context("Input path has no filename")?
        .to_str()
        .context("Input filename is not valid UTF-8")?;
    // sample_id resolution: --sample-id wins, then --sample-name (inherited),
    // then BAM stem with QC suffixes stripped.
    let sample_name = params
        .sample_id_override
        .or(params.sample_name_override)
        .map(|s| s.to_string())
        .unwrap_or_else(|| bam_stem.to_string());

    let bam_start = Instant::now();
    let config = params.config;
    let outdir = params.outdir;

    // Track output files for the summary
    let mut written_outputs: Vec<(String, String)> = Vec::new();

    // === Build RSeQC config and annotations ===
    let rseqc_config = RseqcConfig {
        mapq_cut: params.mapq_cut,
        infer_experiment_sample_size: params.infer_experiment_sample_size,
        min_intron: params.min_intron,
        junction_saturation_min_coverage: params.junction_saturation_min_coverage as u32,
        junction_saturation_sample_start: params.junction_saturation_percentile_floor as u32,
        junction_saturation_sample_end: params.junction_saturation_percentile_ceiling as u32,
        junction_saturation_sample_step: params.junction_saturation_percentile_step as u32,
        inner_distance_sample_size: params.inner_distance_sample_size,
        inner_distance_lower_bound: params.inner_distance_lower_bound,
        inner_distance_upper_bound: params.inner_distance_upper_bound,
        inner_distance_step: params.inner_distance_step,
        bam_stat_enabled: config.bam_stat.enabled
            || config.flagstat.enabled
            || config.idxstats.enabled
            || config.samtools_stats.enabled,
        infer_experiment_enabled: config.infer_experiment.enabled && params.gene_model.is_some(),
        read_duplication_enabled: config.read_duplication.enabled,
        read_distribution_enabled: config.read_distribution.enabled && params.rd_regions.is_some(),
        junction_annotation_enabled: config.junction_annotation.enabled
            && params.ref_junctions.is_some(),
        junction_saturation_enabled: config.junction_saturation.enabled
            && params.known_junctions.is_some(),
        inner_distance_enabled: config.inner_distance.enabled
            && params.exon_bitset.is_some()
            && params.transcript_tree.is_some(),
        tin_enabled: config.tin.enabled && params.tin_index.is_some(),
        tin_sample_size: params.tin_sample_size,
        tin_min_coverage: params.tin_min_coverage,
        tin_seed: config.tin.seed,
        junction_saturation_seed: config.junction_saturation.seed.unwrap_or(42),
        preseq_enabled: config.preseq.enabled,
        preseq_max_segment_length: config.preseq.max_segment_length,
    };

    let rseqc_annotations = RseqcAnnotations {
        gene_model: params.gene_model,
        ref_junctions: params.ref_junctions,
        rd_regions: params.rd_regions,
        exon_bitset: params.exon_bitset,
        transcript_tree: params.transcript_tree,
        tin_index: params.tin_index,
    };

    let any_rseqc_enabled = rseqc_config.bam_stat_enabled
        || rseqc_config.infer_experiment_enabled
        || rseqc_config.read_duplication_enabled
        || rseqc_config.read_distribution_enabled
        || rseqc_config.junction_annotation_enabled
        || rseqc_config.junction_saturation_enabled
        || rseqc_config.inner_distance_enabled
        || rseqc_config.preseq_enabled
        || rseqc_config.tin_enabled;

    // Phase 2 fragmentomics: end-motifs additionally require --fasta; that
    // gating happens inside FragmentomicsAccumulators::new.
    let fragmentomics_config = rna::fragmentomics::FragmentomicsConfig {
        paired_end: params.paired,
        mapq_cut: params.mapq_cut,
    };

    // === Build Qualimap exon index (if enabled) ===
    let qualimap_index = if params.config.qualimap.enabled {
        Some(rna::qualimap::QualimapIndex::from_genes(genes))
    } else {
        None
    };

    // Per-gene merged-exon and intron geometry, indexed by GeneIdx. Built
    // once per run; consumed by the per-gene Tier-2 Parquet writer.
    let per_gene_shape_index = rna::per_gene::shape::GeneShapeIndex::build(genes);

    // Phase 4 SNP fingerprint: load and resolve the SNP-site panel against
    // the BAM header before count_reads() so the worker dispatch can use it.
    // The bundled smoke panel is loaded by default; users override via
    // `--snp-panel <TSV>`.
    let snp_site_index = {
        use rust_htslib::bam::Read;
        let reader = rust_htslib::bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM for SNP-panel resolve: {}", bam_path))?;
        let header = reader.header().clone();
        let user_panel = params.snp_panel.map(std::path::Path::new);
        rna::snp_fingerprint::SnpSiteIndex::build(&header, user_panel)?
    };

    // === dupRadar counting ===
    ui.blank();
    ui.section(&format!("Processing {}", bam_path));
    let pb = ui.progress_bar();
    let count_start = Instant::now();
    let mut count_result = rna::dupradar::counting::count_reads(
        bam_path,
        genes,
        params.stranded,
        params.paired,
        threads,
        params.chrom_mapping,
        params.chrom_prefix,
        params.fasta,
        params.skip_dup_check,
        params.biotype_attribute,
        if any_rseqc_enabled {
            Some(&rseqc_config)
        } else {
            None
        },
        if any_rseqc_enabled {
            Some(&rseqc_annotations)
        } else {
            None
        },
        qualimap_index.as_ref(),
        Some(&fragmentomics_config),
        params.fasta.map(std::path::Path::new),
        Some(&per_gene_shape_index),
        Some(&snp_site_index),
        Some(&pb),
    )?;
    let count_duration = count_start.elapsed();
    ui.finish_progress(&pb, count_result.stat_total_reads, count_duration);

    // Extract RSeQC accumulators from count_result.
    let rseqc_accums = count_result.rseqc.take();

    // Phase 2/3/4: finalize merged sample-level accumulators (fragmentomics
    // proper plus the Phase 4 chrom-metrics and cycle-quality riders).
    // BAM header refs are needed to resolve tid → contig name in the
    // chrom-metrics finalize, so we capture the result inline before the
    // unrelated `bam_header_refs` block below.
    let frag_outputs = count_result.fragmentomics.take().map(|frag| {
        let fragment_size = frag.fragment_size.map(|a| a.into_result());
        let periodicity = fragment_size
            .as_ref()
            .map(|fs| rna::fragmentomics::compute_periodicity(&fs.histogram));
        let end_motifs = frag.end_motifs.map(|a| a.into_result());
        let soft_clips = frag.soft_clips.map(|a| a.into_result());
        let cycle_quality = frag.cycle_quality.map(|a| a.into_result(params.paired));
        (
            fragment_size,
            end_motifs,
            soft_clips,
            periodicity,
            frag.chrom_metrics,
            cycle_quality,
        )
    });
    let (
        frag_size_result,
        end_motifs_result,
        soft_clips_result,
        periodicity_result,
        chrom_metrics_accum,
        cycle_quality_result,
    ) = match frag_outputs {
        Some((a, b, c, d, e, f)) => (a, b, c, d, e, f),
        None => (None, None, None, None, None, None),
    };

    // Apply the `--min-gene-reads` threshold and assemble per-gene Parquet rows.
    let per_gene_output = count_result
        .per_gene
        .take()
        .map(|pg| pg.finalize(params.min_gene_reads, genes, &count_result.gene_counts));

    // Summary stats for the box
    let total_mapped = count_result.stat_total_mapped;
    let total_dup = count_result.stat_total_dup;

    // Summary box
    ui.summary_box(
        &format!("{} — Counting Summary", sample_name),
        &[
            (
                "Total reads:",
                format_count(count_result.stat_total_reads),
                format!(
                    "({} fragments)",
                    format_count(count_result.stat_total_fragments)
                ),
            ),
            (
                "Assigned:",
                format_count(count_result.stat_assigned),
                format_pct(
                    count_result.stat_assigned,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "No features:",
                format_count(count_result.stat_no_features),
                format_pct(
                    count_result.stat_no_features,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "Ambiguous:",
                format_count(count_result.stat_ambiguous),
                format_pct(
                    count_result.stat_ambiguous,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "Duplicates:",
                format_count(total_dup),
                format_pct(total_dup, total_mapped),
            ),
            (
                "Multimappers:",
                format_count(count_result.fc_multimapping),
                format_pct(count_result.fc_multimapping, count_result.stat_total_reads),
            ),
        ],
    );

    // Output directories: nested subdirectories by default, flat if requested
    let fc_dir = if params.flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("featurecounts")
    };
    let dr_dir = if params.flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("dupradar")
    };

    // === featureCounts outputs ===
    ui.section("Writing outputs...");
    if config.any_featurecounts_output() {
        std::fs::create_dir_all(&fc_dir).with_context(|| {
            format!(
                "Failed to create featurecounts output directory: {}",
                fc_dir.display()
            )
        })?;

        if config.featurecounts.counts_file {
            let counts_path = fc_dir.join(format!("{}.featureCounts.tsv", sample_name));
            rna::featurecounts::output::write_counts_file(
                &counts_path,
                genes,
                &count_result,
                bam_path,
                params.command_line,
            )?;
            let p = counts_path.display().to_string();
            ui.output_item("featureCounts", &p);
            written_outputs.push(("featureCounts".into(), p));
        }

        if config.featurecounts.summary_file {
            let summary_path = fc_dir.join(format!("{}.featureCounts.tsv.summary", sample_name));
            rna::featurecounts::output::write_summary_file(&summary_path, &count_result, bam_path)?;
            let p = summary_path.display().to_string();
            ui.output_detail(&format!("Summary: {p}"));
            written_outputs.push(("featureCounts summary".into(), p));
        }

        // Biotype outputs (only if attribute was found in GTF)
        if params.biotype_in_gtf && config.any_biotype_output() {
            let biotype_counts =
                rna::featurecounts::output::aggregate_biotype_counts(&count_result);
            ui.detail(&format!(
                "Biotype counting: {} biotypes found",
                biotype_counts.len()
            ));

            if config.featurecounts.biotype_summary_file {
                let biotype_summary_path =
                    fc_dir.join(format!("{}.featureCounts.biotype.tsv.summary", sample_name));
                rna::featurecounts::output::write_biotype_summary_file(
                    &biotype_summary_path,
                    &count_result,
                    bam_path,
                )?;
                let p = biotype_summary_path.display().to_string();
                ui.output_detail(&format!("Biotype summary: {p}"));
                written_outputs.push(("biotype summary".into(), p));
            }

            if config.featurecounts.biotype_counts {
                let biotype_path = fc_dir.join(format!("{}.biotype_counts.tsv", sample_name));
                rna::featurecounts::output::write_biotype_counts(&biotype_path, &biotype_counts)?;
                let p = biotype_path.display().to_string();
                ui.output_detail(&format!("Biotype counts: {p}"));
                written_outputs.push(("biotype counts".into(), p));
            }

            if config.featurecounts.biotype_counts_mqc {
                let mqc_biotype_path =
                    fc_dir.join(format!("{}.biotype_counts_mqc.tsv", sample_name));
                rna::featurecounts::output::write_biotype_counts_mqc(
                    &mqc_biotype_path,
                    &biotype_counts,
                )?;
                let p = mqc_biotype_path.display().to_string();
                ui.output_detail(&format!("Biotype MultiQC: {p}"));
                written_outputs.push(("biotype MultiQC".into(), p));
            }

            if config.featurecounts.biotype_rrna_mqc {
                let mqc_rrna_path =
                    fc_dir.join(format!("{}.biotype_counts_rrna_mqc.tsv", sample_name));
                rna::featurecounts::output::write_biotype_rrna_mqc(
                    &mqc_rrna_path,
                    &biotype_counts,
                    count_result.fc_biotype_assigned,
                    &sample_name,
                )?;
                let p = mqc_rrna_path.display().to_string();
                ui.output_detail(&format!("rRNA MultiQC: {p}"));
                written_outputs.push(("rRNA MultiQC".into(), p));
            }
        }
    }

    // === dupRadar outputs ===
    let mut dupradar_fit: Option<rna::dupradar::fitting::FitResult> = None;
    let mut dupradar_genes_with_reads: u64 = 0;
    let mut dupradar_genes_with_dups: u64 = 0;

    if config.any_dupradar_output() {
        std::fs::create_dir_all(&dr_dir).with_context(|| {
            format!(
                "Failed to create dupradar output directory: {}",
                dr_dir.display()
            )
        })?;
        let dup_matrix = rna::dupradar::dupmatrix::DupMatrix::build(genes, &count_result);

        let stats = dup_matrix.get_stats();
        dupradar_genes_with_reads = stats.n_regions_covered as u64;
        dupradar_genes_with_dups = stats.n_regions_duplication as u64;

        // Write duplication matrix
        if config.dupradar.dup_matrix {
            let matrix_path = dr_dir.join(format!("{}_dupMatrix.txt", sample_name));
            dup_matrix.write_tsv(&matrix_path)?;
            let p = matrix_path.display().to_string();
            written_outputs.push(("dupRadar matrix".into(), p));
        }

        // Fit logistic regression model (needed for intercept/slope, density plot, and MultiQC)
        let need_fit = config.dupradar.intercept_slope
            || config.dupradar.density_scatter_plot
            || config.dupradar.multiqc_intercept
            || config.dupradar.multiqc_curve;

        if need_fit {
            let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
            let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

            match rna::dupradar::fitting::duprate_exp_fit(&rpk_values, &dup_rate_values) {
                Ok(fit) => {
                    ui.output_detail(&format!(
                        "Model fit: intercept={:.6}, slope={:.6}",
                        fit.intercept, fit.slope
                    ));
                    if config.dupradar.intercept_slope {
                        let fit_path = dr_dir.join(format!("{}_intercept_slope.txt", sample_name));
                        rna::dupradar::plots::write_intercept_slope(&fit, &sample_name, &fit_path)?;
                        let p = fit_path.display().to_string();
                        ui.output_detail(&format!("Fit results: {p}"));
                        written_outputs.push(("dupRadar fit".into(), p));
                    }
                    dupradar_fit = Some(fit);
                }
                Err(e) => {
                    ui.warn(&format!("Could not fit dupRadar model: {}", e));
                }
            }
        }
        let fit_ok = dupradar_fit.as_ref();

        // Generate plots (in parallel — all plots read shared immutable data)
        let any_plot = config.dupradar.density_scatter_plot
            || config.dupradar.boxplot
            || config.dupradar.expression_histogram;

        if any_plot {
            let rpkm_threshold = 0.5;
            let rpkm_threshold_rpk = fit_ok.and_then(|_fit| {
                let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
                let rpkm_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpkm).collect();
                rna::dupradar::fitting::compute_rpkm_threshold_rpk(
                    &rpk_values,
                    &rpkm_values,
                    rpkm_threshold,
                )
            });

            let density_path = dr_dir.join(format!("{}_duprateExpDens.png", sample_name));
            let boxplot_path = dr_dir.join(format!("{}_duprateExpBoxplot.png", sample_name));
            let histogram_path = dr_dir.join(format!("{}_expressionHist.png", sample_name));

            let sample_name_str = sample_name.as_str();
            std::thread::scope(|s| -> Result<()> {
                // Density scatter plot (only if fit succeeded and enabled)
                let density_handle = if config.dupradar.density_scatter_plot {
                    fit_ok.map(|fit| {
                        let dm_ref = &dup_matrix;
                        let thresh = rpkm_threshold_rpk;
                        let path = &density_path;
                        s.spawn(move || {
                            rna::dupradar::plots::density_scatter_plot(
                                dm_ref,
                                fit,
                                thresh,
                                rpkm_threshold,
                                sample_name_str,
                                path,
                            )
                        })
                    })
                } else {
                    None
                };

                // Boxplot
                let boxplot_handle = if config.dupradar.boxplot {
                    let dm_ref = &dup_matrix;
                    let path = &boxplot_path;
                    Some(s.spawn(move || {
                        rna::dupradar::plots::duprate_boxplot(dm_ref, sample_name_str, path)
                    }))
                } else {
                    None
                };

                // Histogram
                let histogram_handle = if config.dupradar.expression_histogram {
                    let dm_ref = &dup_matrix;
                    let path = &histogram_path;
                    Some(s.spawn(move || {
                        rna::dupradar::plots::expression_histogram(dm_ref, sample_name_str, path)
                    }))
                } else {
                    None
                };

                // Collect results
                if let Some(handle) = density_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("density scatter plot thread panicked"))??;
                }
                if let Some(handle) = boxplot_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("boxplot thread panicked"))??;
                }
                if let Some(handle) = histogram_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("histogram thread panicked"))??;
                }

                Ok(())
            })?;

            let plots_dir = dr_dir.display().to_string();
            written_outputs.push(("dupRadar plots".into(), plots_dir));
        }

        // Write MultiQC-compatible output files
        if let Some(fit) = fit_ok {
            if config.dupradar.multiqc_intercept {
                let mqc_intercept_path =
                    dr_dir.join(format!("{}_dup_intercept_mqc.txt", sample_name));
                rna::dupradar::plots::write_mqc_intercept(fit, &sample_name, &mqc_intercept_path)?;
                written_outputs.push((
                    "dupRadar MultiQC".into(),
                    mqc_intercept_path.display().to_string(),
                ));
            }

            if config.dupradar.multiqc_curve {
                let mqc_curve_path =
                    dr_dir.join(format!("{}_duprateExpDensCurve_mqc.txt", sample_name));
                rna::dupradar::plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
                written_outputs.push((
                    "dupRadar MultiQC curve".into(),
                    mqc_curve_path.display().to_string(),
                ));
            }
        }

        // Consolidated dupRadar output line
        ui.output_item("dupRadar", &format!("{}/*", dr_dir.display()));
        ui.output_detail(&format!(
            "{} genes, {} with reads",
            format_count(stats.n_regions as u64),
            format_count(stats.n_regions_covered as u64),
        ));
        if let Some(fit) = fit_ok {
            ui.output_detail(&format!(
                "Model fit: intercept={:.6}, slope={:.6}",
                fit.intercept, fit.slope,
            ));
        }
    }
    // === Qualimap RNA-Seq QC output ===
    if let (Some(ref qm_result), Some(ref qm_index)) = (&count_result.qualimap, &qualimap_index) {
        let qm_dir = if params.flat_output {
            outdir.to_path_buf()
        } else {
            outdir.join("qualimap")
        };

        rna::qualimap::output::write_qualimap_results(
            qm_result,
            qm_index,
            bam_path,
            params.gtf_path,
            params.stranded,
            &qm_dir,
            &sample_name,
        )?;
        let p = qm_dir.display().to_string();
        ui.output_item("Qualimap", &format!("{p}/*"));
        written_outputs.push(("Qualimap".into(), p));
    }

    // === RSeQC analyses (post-processing of single-pass accumulators) ===
    let rseqc_accums = rseqc_accums.unwrap_or_else(|| {
        ui.detail("No RSeQC tools enabled, skipping");
        RseqcAccumulators::empty()
    });
    // Extract BAM header info (reference names + lengths) for samtools-compatible outputs
    let bam_header_refs = {
        let reader = rust_htslib::bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM for header: {}", bam_path))?;
        let header = reader.header();
        (0..header.target_count())
            .map(|tid| {
                let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                let len = header.target_len(tid).unwrap_or(0);
                (name, len)
            })
            .collect::<Vec<(String, u64)>>()
    };
    // Finalize Phase 4 chrom_metrics now that `bam_header_refs` is in scope
    // (chrom_metrics resolves tid → contig name from this slice).
    let chrom_metrics_result = chrom_metrics_accum.map(|cm| cm.into_result(&bam_header_refs));

    let rseqc_outputs = write_rseqc_outputs(
        bam_path,
        &sample_name,
        params,
        rseqc_accums,
        &bam_header_refs,
    )?;
    written_outputs.extend(rseqc_outputs.written);

    let bam_duration = bam_start.elapsed();

    // === v1 envelope (canonical per-sample output) ===
    let (bam_md5, bam_size_bytes) = params
        .bam_hashes
        .get(bam_path)
        .cloned()
        .with_context(|| format!("internal: missing precomputed BAM hash for {bam_path}"))?;
    let bam_path_abs = absolutize(bam_path);
    let biotype_counts = rna::featurecounts::output::aggregate_biotype_counts(&count_result);
    let biotype_rrna_count = biotype_counts.get("rRNA").copied().unwrap_or(0);

    let qc_flags_vec = qc_flags::evaluate(&qc_flags::QcContext {
        paired_end: params.paired,
        strandedness: params.stranded,
        library_prep: params.library_prep,
        bam_stat: rseqc_outputs.bam_stat.as_ref(),
        count_result: Some(&count_result),
        preseq: rseqc_outputs.preseq.as_ref(),
        dupradar_fit: dupradar_fit.as_ref(),
        dupradar_genes_with_reads,
        featurecounts_biotype_rrna: biotype_rrna_count,
        fragment_size: frag_size_result.as_ref(),
        read_distribution: rseqc_outputs.read_distribution.as_ref(),
    });

    // Phase 4 gene-class fractions: bundled hemoglobin / RP / apolipoprotein
    // symbol lists, summed against `count_result.gene_counts.fc_reads` and
    // normalized to `fc_assigned`. Always emitted (only blank when the GTF
    // carried no `gene_name` attribute and no genes matched).
    let gene_class_index = rna::gene_class::GeneClassIndex::build(genes);
    let gene_class_fractions =
        gene_class_index.finalize(&count_result.gene_counts, count_result.fc_assigned);

    // Phase 4 marker-panel aggregation: bundled lm22 / tabula_sapiens_cfrna /
    // vorperian by default, plus optional user TSV via `--panels-tsv`.
    // `--panels` CSV picks among the bundled panels; default = all.
    let panel_selection: Option<Vec<&str>> = params.panels.map(|csv| {
        csv.split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect()
    });
    let panel_user_paths: Vec<&Path> = params
        .panels_tsv
        .iter()
        .map(|p| Path::new(p.as_str()))
        .collect();
    let panel_index =
        rna::panels::PanelIndex::build(genes, panel_selection.as_deref(), &panel_user_paths)?;
    let panels_result = if panel_index.is_empty() {
        None
    } else {
        Some(panel_index.aggregate(&count_result.gene_counts, count_result.fc_assigned))
    };

    // Phase 4 sex inference: pure finalize-time computation against
    // chrom-metrics + gene_counts. No CLI metadata input — emits a
    // `predicted_sex` label and supporting metrics; downstream consumers
    // compare against their own subject metadata.
    let sex_inference_result = rna::sex_infer::compute(
        chrom_metrics_result.as_ref(),
        genes,
        &count_result.gene_counts,
        count_result.fc_assigned,
    );

    // Phase 4 saturation curve: deterministic per-read hash bucketing,
    // computed at the requested fraction grid. Fragments accumulated by
    // the dispatcher; finalize is a pure traversal of the merged buckets.
    let saturation_fractions: Vec<f64> = match params.saturation_fractions {
        Some(csv) => csv
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse::<f64>().ok())
            .collect(),
        None => rna::saturation::DEFAULT_FRACTIONS.to_vec(),
    };
    let saturation_result = count_result
        .saturation
        .take()
        .map(|s| s.finalize(&saturation_fractions));

    // Phase 4 SNP fingerprint finalize: emit per-site allele counts +
    // fractions. Skipped when no sites resolved against the BAM header.
    let snp_fingerprint_result = match (count_result.snp_fingerprint.take(), &snp_site_index) {
        (Some(accum), idx) if !idx.is_empty() => Some(accum.finalize(idx)),
        _ => None,
    };

    // Phase 4 splice-site dinucleotides — finalize-time tally over the
    // unique junctions in junction_annotation. Requires --fasta. Skipped
    // silently when the user didn't supply one (no qc_flag added; the
    // existing `end_motifs_skipped_no_fasta` already signals that case).
    let splice_dinuc_result = match (rseqc_outputs.junction_annotation.as_ref(), params.fasta) {
        (Some(jr), Some(fasta_path)) => {
            match rna::splice_dinuc::compute(jr, std::path::Path::new(fasta_path)) {
                Ok(r) => Some(r),
                Err(e) => {
                    ui.detail(&format!("splice_site_dinucleotides skipped: {}", e));
                    None
                }
            }
        }
        _ => None,
    };

    // Write the per-gene Parquet sibling and assemble its envelope block.
    let per_gene_block = if let Some(output) = per_gene_output.as_ref() {
        let parquet_filename = format!("{sample_name}.per_gene.parquet");
        let parquet_path = outdir.join(&parquet_filename);
        let strandedness_str = params.stranded.to_string();
        let meta = rna::per_gene::PerGeneFileMeta {
            sample_id: sample_name.clone(),
            extractor_version: env!("CARGO_PKG_VERSION").to_string(),
            git_commit: env!("GIT_SHORT_HASH").to_string(),
            per_gene_schema_version: rna::per_gene::writer::PER_GENE_SCHEMA_VERSION.to_string(),
            bam_md5: bam_md5.clone(),
            gtf_md5: params.gtf_md5.to_string(),
            library_prep: params.library_prep.to_string(),
            min_gene_reads_threshold: output.min_gene_reads_threshold,
            n_genes_total: output.n_genes_total,
            n_genes_emitted: output.n_genes_emitted,
            paired_end: params.paired,
            strandedness: strandedness_str,
        };
        let parquet_md5 =
            rna::per_gene::writer::write_per_gene_parquet(&output.rows, &meta, &parquet_path)?;
        ui.output_item("per_gene_parquet", &parquet_path.display().to_string());
        Some(envelope::PerGeneBlock {
            parquet_path: parquet_filename,
            parquet_md5,
            n_genes_total: output.n_genes_total,
            n_genes_emitted: output.n_genes_emitted,
            min_gene_reads_threshold: output.min_gene_reads_threshold,
            schema_version_per_gene: rna::per_gene::writer::PER_GENE_SCHEMA_VERSION.to_string(),
            columns: rna::per_gene::column_names(),
        })
    } else {
        None
    };

    let envelope_value = envelope::build(envelope::BuildInputs {
        extractor_version: env!("CARGO_PKG_VERSION"),
        git_commit: env!("GIT_SHORT_HASH"),
        sample_id: &sample_name,
        bam_path: &bam_path_abs,
        bam_md5: &bam_md5,
        bam_size_bytes,
        gtf_path: params.gtf_path_abs,
        gtf_md5: params.gtf_md5,
        reference_fasta_md5: params.reference_fasta_md5,
        library_prep: params.library_prep,
        paired_end: params.paired,
        strandedness: params.stranded,
        filters: envelope::Filters::from_mapq(params.mapq_cut),
        runtime_seconds: bam_duration.as_secs_f64(),
        peak_rss_mb: runtime_stats::peak_rss_mb(),
        qc_flags: qc_flags_vec,
        bam_stat: rseqc_outputs.bam_stat.as_ref(),
        infer_experiment: rseqc_outputs.infer_experiment.as_ref(),
        read_distribution: rseqc_outputs.read_distribution.as_ref(),
        read_duplication: rseqc_outputs.read_duplication.as_ref(),
        junction_annotation: rseqc_outputs.junction_annotation.as_ref(),
        junction_saturation: rseqc_outputs.junction_saturation.as_ref(),
        inner_distance: rseqc_outputs.inner_distance.as_ref(),
        tin: rseqc_outputs.tin.as_ref(),
        preseq: rseqc_outputs.preseq.as_ref(),
        qualimap: count_result.qualimap.as_ref(),
        count_result: Some(&count_result),
        dupradar_fit: dupradar_fit.as_ref(),
        dupradar_genes_with_reads,
        dupradar_genes_with_duplication: dupradar_genes_with_dups,
        fragment_size: frag_size_result.as_ref(),
        end_motifs: end_motifs_result.as_ref(),
        soft_clips: soft_clips_result.as_ref(),
        periodicity: periodicity_result.as_ref(),
        fasta_provided: params.fasta.is_some(),
        per_gene: per_gene_block.as_ref(),
        splice_site_dinucleotides: splice_dinuc_result.as_ref(),
        per_chromosome: chrom_metrics_result.as_ref(),
        cycle_quality: cycle_quality_result.as_ref(),
        gene_class_fractions: Some(&gene_class_fractions),
        panels: panels_result.as_ref(),
        sex_inference: Some(&sex_inference_result),
        saturation: saturation_result.as_ref(),
        snp_fingerprint: snp_fingerprint_result.as_ref(),
    });
    let envelope_path = outdir.join(format!("{sample_name}.liquidqc.json"));
    envelope::write(&envelope_value, &envelope_path)?;
    ui.output_item("envelope", &envelope_path.display().to_string());

    ui.finish(bam_stem, bam_duration);

    Ok(BamResult {
        duration: bam_duration,
        infer_experiment: rseqc_outputs.infer_experiment,
    })
}

/// Convert a possibly-relative path string into an absolute string without
/// resolving symlinks (so `/private/var/...` doesn't surface on macOS).
fn absolutize(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string())
}

// ============================================================================
// RSeQC output writing (post-processing of single-pass accumulators)
// ============================================================================

/// Results returned from `write_rseqc_outputs`, bundling output file paths
/// with the per-tool result structs needed for the v1 envelope and downstream
/// checks (e.g. strandedness mismatch).
struct RseqcOutputs {
    /// Output files written (tool, path).
    written: Vec<(String, String)>,
    /// bam_stat result (also feeds samtools-equivalent outputs above).
    bam_stat: Option<rna::rseqc::bam_stat::BamStatResult>,
    /// infer_experiment result, if the tool was enabled and produced data.
    infer_experiment: Option<rna::rseqc::infer_experiment::InferExperimentResult>,
    read_distribution: Option<rna::rseqc::read_distribution::ReadDistributionResult>,
    read_duplication: Option<rna::rseqc::read_duplication::ReadDuplicationResult>,
    junction_annotation: Option<rna::rseqc::junction_annotation::JunctionResults>,
    junction_saturation: Option<rna::rseqc::junction_saturation::SaturationResult>,
    inner_distance: Option<rna::rseqc::inner_distance::InnerDistanceResult>,
    tin: Option<rna::rseqc::tin::TinResults>,
    preseq: Option<rna::preseq::PreseqResult>,
}

/// Write all RSeQC outputs from the single-pass accumulators.
///
/// Converts accumulated data to tool-specific result types and writes all
/// output files, plots, and summaries. Returns the written output paths and
/// the infer_experiment result (if available) for strandedness mismatch checking.
fn write_rseqc_outputs(
    bam_path: &str,
    sample_name: &str,
    params: &SharedParams,
    accums: RseqcAccumulators,
    bam_header_refs: &[(String, u64)],
) -> Result<RseqcOutputs> {
    let ui = params.ui;
    let outdir = params.outdir;
    let mut written: Vec<(String, String)> = Vec::new();
    let mut infer_experiment_result: Option<rna::rseqc::infer_experiment::InferExperimentResult> =
        None;
    let mut read_distribution_result: Option<
        rna::rseqc::read_distribution::ReadDistributionResult,
    > = None;
    let mut read_duplication_result: Option<rna::rseqc::read_duplication::ReadDuplicationResult> =
        None;
    let mut junction_annotation_result: Option<rna::rseqc::junction_annotation::JunctionResults> =
        None;
    let mut junction_saturation_result: Option<rna::rseqc::junction_saturation::SaturationResult> =
        None;
    let mut inner_distance_result: Option<rna::rseqc::inner_distance::InnerDistanceResult> = None;
    let mut tin_result: Option<rna::rseqc::tin::TinResults> = None;
    let mut preseq_result: Option<rna::preseq::PreseqResult> = None;

    // Build tool-specific output directories: nested subdirectories by default, flat if requested
    let flat = params.flat_output;
    let rseqc_bam_stat_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("bam_stat")
    };
    let rseqc_read_dup_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("read_duplication")
    };
    let rseqc_infer_exp_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("infer_experiment")
    };
    let rseqc_read_dist_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("read_distribution")
    };
    let rseqc_junc_annot_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("junction_annotation")
    };
    let rseqc_junc_sat_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("junction_saturation")
    };
    let rseqc_inner_dist_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("inner_distance")
    };

    // --- samtools-compatible outputs (flagstat, idxstats, stats) ---
    let samtools_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("samtools")
    };

    // Compute bam_stat result once — used by both samtools and RSeQC outputs.
    let bam_stat_result = accums.bam_stat.map(|accum| accum.into_result());

    if let Some(ref result) = bam_stat_result {
        let has_samtools = params.config.flagstat.enabled
            || params.config.idxstats.enabled
            || params.config.samtools_stats.enabled;
        if has_samtools {
            ui.output_group("samtools");
        }

        // --- flagstat ---
        if params.config.flagstat.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let flagstat_path = samtools_dir.join(format!("{}.flagstat", sample_name));
            rna::rseqc::flagstat::write_flagstat(result, &flagstat_path)?;
            let p = flagstat_path.display().to_string();
            ui.output_item("flagstat", &p);
            written.push(("flagstat".into(), p));
        }

        // --- idxstats ---
        if params.config.idxstats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let idxstats_path = samtools_dir.join(format!("{}.idxstats", sample_name));
            rna::rseqc::idxstats::write_idxstats(result, bam_header_refs, &idxstats_path)?;
            let p = idxstats_path.display().to_string();
            ui.output_item("idxstats", &p);
            written.push(("idxstats".into(), p));
        }

        // --- samtools stats SN ---
        if params.config.samtools_stats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let stats_path = samtools_dir.join(format!("{}.stats", sample_name));
            rna::rseqc::stats::write_stats(result, &stats_path)?;
            let p = stats_path.display().to_string();
            ui.output_item("stats", &p);
            written.push(("samtools stats".into(), p));
        }
    }

    // --- RSeQC outputs ---
    {
        let has_rseqc = (bam_stat_result.is_some() && params.config.bam_stat.enabled)
            || accums.read_dup.is_some()
            || accums.infer_exp.is_some()
            || accums.read_dist.is_some()
            || accums.junc_annot.is_some()
            || accums.junc_sat.is_some()
            || accums.inner_dist.is_some()
            || accums.tin.is_some();
        if has_rseqc {
            ui.output_group("RSeQC");
        }
    }

    // --- bam_stat ---
    if let Some(ref result) = bam_stat_result {
        if params.config.bam_stat.enabled {
            std::fs::create_dir_all(&rseqc_bam_stat_dir)?;
            let output_path = rseqc_bam_stat_dir.join(format!("{}.bam_stat.txt", sample_name));
            rna::rseqc::bam_stat::write_bam_stat(result, &output_path)?;
            let p = output_path.display().to_string();
            ui.output_item("bam_stat", &p);
            written.push(("bam_stat".into(), p));
        }
    }

    // --- read_duplication ---
    if let Some(accum) = accums.read_dup {
        std::fs::create_dir_all(&rseqc_read_dup_dir)?;
        let result = accum.into_result();
        rna::rseqc::read_duplication::write_read_duplication(
            &result,
            &rseqc_read_dup_dir,
            sample_name,
        )?;
        let plot_path = rseqc_read_dup_dir.join(format!("{}.DupRate_plot.png", sample_name));
        rna::rseqc::plots::read_duplication_plot(&result, sample_name, &plot_path)?;
        let p = rseqc_read_dup_dir.display().to_string();
        ui.output_item("read_duplication", &format!("{p}/{sample_name}.*"));
        written.push(("read_duplication".into(), p));
        read_duplication_result = Some(result);
    }

    // --- infer_experiment ---
    if let Some(accum) = accums.infer_exp {
        std::fs::create_dir_all(&rseqc_infer_exp_dir)?;
        let result = accum.into_result();
        let output_path = rseqc_infer_exp_dir.join(format!("{}.infer_experiment.txt", sample_name));
        rna::rseqc::infer_experiment::write_infer_experiment(&result, &output_path)?;
        let p = output_path.display().to_string();
        ui.output_item("infer_experiment", &p);
        ui.output_detail(&format!(
            "{} usable reads sampled",
            format_count(result.total_sampled),
        ));
        written.push(("infer_experiment".into(), p));
        infer_experiment_result = Some(result);
    }

    // --- read_distribution ---
    if let Some(accum) = accums.read_dist {
        std::fs::create_dir_all(&rseqc_read_dist_dir)?;
        let rd_regions = params
            .rd_regions
            .context("rd_regions must be Some when read_distribution accumulator exists")?;
        let result = accum.into_result(rd_regions);
        let output_path =
            rseqc_read_dist_dir.join(format!("{}.read_distribution.txt", sample_name));
        rna::rseqc::read_distribution::write_read_distribution(&result, &output_path)?;
        let p = output_path.display().to_string();
        ui.output_item("read_distribution", &p);
        ui.output_detail(&format!(
            "{} reads, {} tags, {} assigned",
            format_count(result.total_reads),
            format_count(result.total_tags),
            format_count(result.total_tags - result.unassigned_tags),
        ));
        written.push(("read_distribution".into(), p));
        read_distribution_result = Some(result);
    }

    // --- junction_annotation ---
    if let Some(accum) = accums.junc_annot {
        std::fs::create_dir_all(&rseqc_junc_annot_dir)?;
        let prefix = rseqc_junc_annot_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(bam_header_refs);

        let xls_path = rseqc_junc_annot_dir.join(format!("{}.junction.xls", sample_name));
        rna::rseqc::junction_annotation::write_junction_xls(&results, &xls_path)?;

        let bed_out_path = rseqc_junc_annot_dir.join(format!("{}.junction.bed", sample_name));
        rna::rseqc::junction_annotation::write_junction_bed(&results, &bed_out_path)?;

        let interact_path =
            rseqc_junc_annot_dir.join(format!("{}.junction.Interact.bed", sample_name));
        rna::rseqc::junction_annotation::write_junction_interact_bed(
            &results,
            bam_path,
            &interact_path,
        )?;

        let r_path = rseqc_junc_annot_dir.join(format!("{}.junction_plot.r", sample_name));
        rna::rseqc::junction_annotation::write_junction_plot_r(&results, &prefix, &r_path)?;

        rna::rseqc::plots::junction_annotation_plot(&results, &prefix, sample_name)?;

        let summary_path =
            rseqc_junc_annot_dir.join(format!("{}.junction_annotation.log", sample_name));
        rna::rseqc::junction_annotation::write_summary(&results, &summary_path, params.gtf_path)?;

        // Only print the detailed junction summary in verbose mode
        if ui.is_verbose() {
            rna::rseqc::junction_annotation::print_summary(&results);
        }

        let p = rseqc_junc_annot_dir.display().to_string();
        ui.output_item("junction_annotation", &format!("{p}/{sample_name}.*"));
        written.push(("junction_annotation".into(), p));
        junction_annotation_result = Some(results);
    }

    // --- junction_saturation ---
    if let Some(accum) = accums.junc_sat {
        std::fs::create_dir_all(&rseqc_junc_sat_dir)?;
        let prefix = rseqc_junc_sat_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let known_junctions = params
            .known_junctions
            .context("known_junctions must be Some when junction_saturation accumulator exists")?;
        let results = accum.into_result(
            known_junctions,
            params.junction_saturation_percentile_floor as u32,
            params.junction_saturation_percentile_ceiling as u32,
            params.junction_saturation_percentile_step as u32,
            params.junction_saturation_min_coverage as u32,
            params.junction_saturation_seed,
        );

        rna::rseqc::junction_saturation::write_r_script(&results, &prefix)?;

        let plot_path =
            rseqc_junc_sat_dir.join(format!("{}.junctionSaturation_plot.png", sample_name));
        rna::rseqc::plots::junction_saturation_plot(&results, sample_name, &plot_path)?;

        let summary_path = format!("{prefix}.junctionSaturation_summary.txt");
        rna::rseqc::junction_saturation::write_summary(&results, &summary_path)?;

        let p = rseqc_junc_sat_dir.display().to_string();
        ui.output_item("junction_saturation", &format!("{p}/{sample_name}.*"));
        written.push(("junction_saturation".into(), p));
        junction_saturation_result = Some(results);
    }

    // --- inner_distance ---
    if let Some(accum) = accums.inner_dist {
        std::fs::create_dir_all(&rseqc_inner_dist_dir)?;
        let prefix = rseqc_inner_dist_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            params.inner_distance_step,
        )?;

        let detail_path = format!("{prefix}.inner_distance.txt");
        rna::rseqc::inner_distance::write_detail_file(&results, &detail_path)?;

        let freq_path = format!("{prefix}.inner_distance_freq.txt");
        rna::rseqc::inner_distance::write_freq_file(&results, &freq_path)?;

        let r_path = format!("{prefix}.inner_distance_plot.r");
        rna::rseqc::inner_distance::write_r_script(
            &results,
            &prefix,
            &r_path,
            params.inner_distance_step,
        )?;

        let plot_path =
            rseqc_inner_dist_dir.join(format!("{}.inner_distance_plot.png", sample_name));
        rna::rseqc::plots::inner_distance_plot(
            &results,
            params.inner_distance_step,
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            sample_name,
            &plot_path,
        )?;

        let summary_path = format!("{prefix}.inner_distance_summary.txt");
        rna::rseqc::inner_distance::write_summary(&results, &summary_path)?;

        let mean_path = format!("{prefix}.inner_distance_mean.txt");
        rna::rseqc::inner_distance::write_mean_file(&results, sample_name, &mean_path)?;

        let p = rseqc_inner_dist_dir.display().to_string();
        ui.output_item("inner_distance", &format!("{p}/{sample_name}.*"));
        ui.output_detail(&format!(
            "{} read pairs processed",
            format_count(results.total_pairs),
        ));
        written.push(("inner_distance".into(), p));
        inner_distance_result = Some(results);
    }

    // --- TIN ---
    if let Some(accum) = accums.tin {
        let rseqc_tin_dir = if flat {
            outdir.to_path_buf()
        } else {
            outdir.join("rseqc").join("tin")
        };
        std::fs::create_dir_all(&rseqc_tin_dir)?;
        let prefix = rseqc_tin_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let tin_index = params
            .tin_index
            .as_ref()
            .context("TIN index must exist when TIN accumulator is present")?;
        let results = accum.into_result(tin_index);
        rna::rseqc::tin::write_tin(&results, Path::new(&format!("{prefix}.tin.xls")))?;
        rna::rseqc::tin::write_tin_summary(
            &results,
            bam_path,
            Path::new(&format!("{prefix}.summary.txt")),
        )?;
        let p = format!("{prefix}.tin.xls");
        ui.output_item("TIN", &p);
        ui.output_detail(&format!(
            "{} transcripts",
            format_count(results.len() as u64)
        ));
        written.push(("TIN".into(), p));
        tin_result = Some(results);
    }

    // --- preseq (library complexity) ---
    if let Some(mut accum) = accums.preseq {
        let preseq_dir = if flat {
            outdir.to_path_buf()
        } else {
            outdir.join("preseq")
        };
        std::fs::create_dir_all(&preseq_dir)?;
        let output_path = preseq_dir.join(format!("{}.lc_extrap.txt", sample_name));
        accum.finalize();
        let total_reads = accum.total_fragments;
        let n_distinct = accum.n_distinct();
        let histogram = accum.into_histogram();
        debug!(
            "preseq: {} histogram bins, {} total reads, {} distinct",
            histogram.len(),
            total_reads,
            n_distinct,
        );
        match rna::preseq::estimate_complexity(
            &histogram,
            total_reads,
            n_distinct,
            &params.config.preseq,
        ) {
            Ok(result) => {
                rna::preseq::write_output(
                    &result,
                    &output_path,
                    params.config.preseq.confidence_level,
                )?;
                let p = output_path.display().to_string();
                ui.output_item("preseq", &p);
                ui.output_detail(&format!("{} extrapolation points", result.curve.len(),));
                written.push(("preseq".into(), p));
                preseq_result = Some(result);
            }
            Err(e) => {
                ui.warn(&format!("preseq: skipped — {}", e));
            }
        }
    }

    Ok(RseqcOutputs {
        written,
        bam_stat: bam_stat_result,
        infer_experiment: infer_experiment_result,
        read_distribution: read_distribution_result,
        read_duplication: read_duplication_result,
        junction_annotation: junction_annotation_result,
        junction_saturation: junction_saturation_result,
        inner_distance: inner_distance_result,
        tin: tin_result,
        preseq: preseq_result,
    })
}
