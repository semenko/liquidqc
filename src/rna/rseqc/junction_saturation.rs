//! Junction saturation analysis for RNA-seq data.
//!
//! Reimplements RSeQC's `junction_saturation.py`: subsamples splice junction
//! observations at increasing percentages of total reads and reports how many
//! known / novel / total unique junctions are detected at each level.

use super::common::{self, KnownJunctionSet};
use anyhow::{Context, Result};
use log::{debug, info};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::io::Write;

// ============================= Saturation result =============================

/// Results from junction saturation analysis at each sampling percentage.
#[derive(Debug)]
pub struct SaturationResult {
    /// Sampling percentages (e.g., 5, 10, ..., 100).
    pub percentages: Vec<u32>,
    /// Number of known junctions detected at each percentage.
    pub known_counts: Vec<usize>,
    /// Number of novel junctions detected at each percentage.
    pub novel_counts: Vec<usize>,
    /// Number of all unique junctions detected at each percentage.
    pub all_counts: Vec<usize>,
}

// ============================= Main analysis =============================

/// Run junction saturation analysis on a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to the input BAM file.
/// * `known_junctions` - Pre-built set of known splice junctions (from BED or GTF).
/// * `min_intron` - Minimum intron size to consider (default: 50).
/// * `mapq_cut` - Minimum mapping quality (default: 30).
/// * `min_coverage` - Minimum supporting reads for known junctions (default: 1).
/// * `sample_start` - Starting percentage (default: 5).
/// * `sample_end` - Ending percentage (default: 100).
/// * `sample_step` - Step between percentages (default: 5).
///
/// # Returns
/// A `SaturationResult` with counts at each sampling level.
#[allow(clippy::too_many_arguments)]
pub fn junction_saturation(
    bam_path: &str,
    known_junctions: &KnownJunctionSet,
    min_intron: u64,
    mapq_cut: u8,
    min_coverage: u32,
    sample_start: u32,
    sample_end: u32,
    sample_step: u32,
    reference: Option<&str>,
) -> Result<SaturationResult> {
    // Build set of reference chromosomes (uppercased)
    let ref_chroms: std::collections::HashSet<String> = known_junctions
        .junctions
        .iter()
        .map(|k| k.split(':').next().unwrap().to_string())
        .collect();

    // Open BAM/CRAM file
    let mut bam = if let Some(ref_path) = reference {
        let mut reader =
            bam::Reader::from_path(bam_path).with_context(|| format!("opening BAM: {bam_path}"))?;
        reader.set_reference(ref_path)?;
        reader
    } else {
        bam::Reader::from_path(bam_path).with_context(|| format!("opening BAM: {bam_path}"))?
    };
    let header = bam.header().clone();

    // Map tid -> uppercase chromosome name
    let tid_to_chrom: HashMap<i32, String> = (0..header.target_count() as i32)
        .map(|tid| {
            let name = String::from_utf8_lossy(header.tid2name(tid as u32)).to_uppercase();
            (tid, name)
        })
        .collect();

    // Phase 1: Collect all splice junction observations
    info!("Collecting splice junction observations from {bam_path}...");
    let mut all_observations: Vec<String> = Vec::new();
    let mut total_reads: u64 = 0;
    let mut skipped_qcfail: u64 = 0;
    let mut skipped_dup: u64 = 0;
    let mut skipped_secondary: u64 = 0;
    let mut skipped_unmapped: u64 = 0;
    let mut skipped_mapq: u64 = 0;

    for result in bam.records() {
        let record = result.context("reading BAM record")?;
        total_reads += 1;

        // Flag filtering (same order as RSeQC)
        if record.flags() & 0x200 != 0 {
            skipped_qcfail += 1;
            continue;
        }
        if record.flags() & 0x400 != 0 {
            skipped_dup += 1;
            continue;
        }
        if record.flags() & 0x100 != 0 {
            skipped_secondary += 1;
            continue;
        }
        if record.is_unmapped() {
            skipped_unmapped += 1;
            continue;
        }
        if record.mapq() < mapq_cut {
            skipped_mapq += 1;
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let chrom = match tid_to_chrom.get(&tid) {
            Some(c) => c,
            None => continue,
        };

        // Skip chromosomes not in reference
        if !ref_chroms.contains(chrom) {
            continue;
        }

        let start = record.pos() as u64;
        let cigar = record.cigar();
        let introns = common::fetch_introns(start, cigar.as_ref());

        for (istart, iend) in introns {
            if iend - istart < min_intron {
                continue;
            }
            let key = format!("{chrom}:{istart}-{iend}");
            all_observations.push(key);
        }
    }

    info!(
        "Total reads: {total_reads}, QC fail: {skipped_qcfail}, Dup: {skipped_dup}, \
         Secondary: {skipped_secondary}, Unmapped: {skipped_unmapped}, Low MAPQ: {skipped_mapq}"
    );
    info!(
        "Total splice junction observations: {}",
        all_observations.len()
    );

    // Phase 2: Shuffle with a deterministic seed for reproducibility
    let mut rng = ChaCha8Rng::seed_from_u64(42);
    all_observations.shuffle(&mut rng);

    // Phase 3: Incremental sampling
    let sr_num = all_observations.len();
    let mut percentages: Vec<u32> = (sample_start..=sample_end)
        .step_by(sample_step as usize)
        .collect();
    if !percentages.contains(&100) {
        percentages.push(100);
    }

    let mut unique_junctions: HashMap<String, u32> = HashMap::new();
    let mut known_counts: Vec<usize> = Vec::new();
    let mut novel_counts: Vec<usize> = Vec::new();
    let mut all_counts: Vec<usize> = Vec::new();

    let mut prev_end: usize = 0;

    for &pct in &percentages {
        let index_end = (sr_num as f64 * (pct as f64 / 100.0)) as usize;

        // Add new observations from this increment
        for obs in &all_observations[prev_end..index_end] {
            *unique_junctions.entry(obs.clone()).or_insert(0) += 1;
        }
        prev_end = index_end;

        // Count junctions at this sampling level
        let all_count = unique_junctions.len();
        let mut known_count = 0usize;
        let mut novel_count = 0usize;

        for (junction, &count) in &unique_junctions {
            if known_junctions.junctions.contains(junction) {
                if count >= min_coverage {
                    known_count += 1;
                }
            } else {
                novel_count += 1;
            }
        }

        all_counts.push(all_count);
        known_counts.push(known_count);
        novel_counts.push(novel_count);

        debug!(
            "{}%: all={}, known={}, novel={}",
            pct, all_count, known_count, novel_count
        );
    }

    Ok(SaturationResult {
        percentages,
        known_counts,
        novel_counts,
        all_counts,
    })
}

// ============================= Output formatting =============================

/// Write the junction saturation R script.
///
/// Produces an R script matching RSeQC's `junctionSaturation_plot.r` format.
pub fn write_r_script(result: &SaturationResult, prefix: &str) -> Result<()> {
    let r_path = format!("{prefix}.junctionSaturation_plot.r");
    let mut f =
        std::fs::File::create(&r_path).with_context(|| format!("creating R script: {r_path}"))?;

    let pdf_path = format!("{prefix}.junctionSaturation_plot.pdf");

    // x values
    let x_str: Vec<String> = result.percentages.iter().map(|p| p.to_string()).collect();

    // y = known, z = all, w = novel
    let y_str: Vec<String> = result.known_counts.iter().map(|c| c.to_string()).collect();
    let z_str: Vec<String> = result.all_counts.iter().map(|c| c.to_string()).collect();
    let w_str: Vec<String> = result.novel_counts.iter().map(|c| c.to_string()).collect();

    // m = max of last values / 1000, n = min of first values / 1000
    let known_last = *result.known_counts.last().unwrap_or(&0) as i64;
    let all_last = *result.all_counts.last().unwrap_or(&0) as i64;
    let novel_last = *result.novel_counts.last().unwrap_or(&0) as i64;
    let known_first = *result.known_counts.first().unwrap_or(&0) as i64;
    let all_first = *result.all_counts.first().unwrap_or(&0) as i64;
    let novel_first = *result.novel_counts.first().unwrap_or(&0) as i64;

    let m = (known_last / 1000)
        .max(all_last / 1000)
        .max(novel_last / 1000);
    let n = (known_first / 1000)
        .min(all_first / 1000)
        .min(novel_first / 1000);

    writeln!(f, "pdf('{pdf_path}')")?;
    writeln!(f, "x=c({})", x_str.join(","))?;
    writeln!(f, "y=c({})", y_str.join(","))?;
    writeln!(f, "z=c({})", z_str.join(","))?;
    writeln!(f, "w=c({})", w_str.join(","))?;
    let all_last_k = all_last / 1000;
    let novel_last_k = novel_last / 1000;
    let all_first_k = all_first / 1000;
    let novel_first_k = novel_first / 1000;
    writeln!(f, "m=max({m},{all_last_k},{novel_last_k})")?;
    writeln!(f, "n=min({n},{all_first_k},{novel_first_k})")?;
    writeln!(
        f,
        "plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))"
    )?;
    writeln!(f, "points(x,y/1000,type='o',col='red')")?;
    writeln!(f, "points(x,w/1000,type='o',col='green')")?;
    writeln!(
        f,
        "legend(5,{}, legend=c(\"All junctions\",\"known junctions\", \"novel junctions\"),col=c(\"blue\",\"red\",\"green\"),lwd=1,pch=1)",
        all_last / 1000
    )?;
    writeln!(f, "dev.off()")?;

    info!("Wrote R script: {r_path}");
    Ok(())
}

/// Write a summary text file with junction saturation statistics.
pub fn write_summary(result: &SaturationResult, path: &str) -> Result<()> {
    let mut f = std::fs::File::create(path).with_context(|| format!("creating summary: {path}"))?;

    writeln!(
        f,
        "Percent\tAll_Junctions\tKnown_Junctions\tNovel_Junctions"
    )?;
    for i in 0..result.percentages.len() {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            result.percentages[i],
            result.all_counts[i],
            result.known_counts[i],
            result.novel_counts[i],
        )?;
    }

    info!("Wrote summary: {path}");
    Ok(())
}

// ============================= Tests =============================

// Note: CIGAR intron extraction tests and BED junction parsing tests have been
// moved to the common module (rseqc::common::tests). Tests here focus on
// saturation-specific logic.
