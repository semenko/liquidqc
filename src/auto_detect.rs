//! Pre-pass auto-detection for library properties that liquidqc would
//! otherwise force the user to declare on the CLI.
//!
//! Both helpers open the BAM independently of the main counting pass and
//! read at most a small bounded sample of records. The cost is sub-second
//! on indexed BAMs and serves to make `liquidqc rna sample.bam` Just Work
//! without `--paired` or `--stranded` for typical inputs.

use crate::cli::Strandedness;
use crate::genome::ChromStyle;
use crate::rna::rseqc::accumulators::InferExpAccum;
use crate::rna::rseqc::infer_experiment::{infer_strandedness, GeneModel, InferredStrandedness};
use anyhow::{Context, Result};
use rust_htslib::bam::{Read as BamRead, Reader};
use std::collections::HashMap;

/// Default number of primary records to sample for paired-end detection.
pub const PAIRED_SAMPLE_SIZE: usize = 1_000;

/// Default number of usable records to sample for strandedness inference.
///
/// Each record must overlap a gene model to count, so the actual reader
/// loop scans more BAM records than this — bounded by [`STRANDED_HARD_CAP`].
pub const STRANDED_SAMPLE_SIZE: u64 = 200_000;

/// Hard cap on records read during the strandedness pre-pass, regardless
/// of how many overlap the gene model. Prevents pathological BAMs (e.g.
/// reads that all sit in unannotated decoys) from running the pre-pass
/// indefinitely.
const STRANDED_HARD_CAP: u64 = 1_000_000;

/// How often to re-check the early-exit condition during the strandedness
/// pre-pass. Summing the accumulator's per-class counts on every record
/// would dominate the loop; checking every 1024 records is amortized free.
const STRANDED_EXIT_CHECK_EVERY: u64 = 1024;

/// Minimum primary-record sample size below which the absence of duplicate
/// flags is not sufficient evidence to auto-disable the dup-check. Below this
/// threshold the existing post-pass error handles the case.
pub const DUP_CHECK_MIN_SAMPLE: usize = 500;

/// Result of the first-pass BAM scan.
///
/// One BAM open covers paired-end detection and detection of whether
/// duplicate marks are present, since both questions are answered from
/// the same primary-record sample.
pub struct FirstPassSample {
    /// True when more than half of primary records carry `BAM_FPAIRED`.
    pub paired: bool,
    /// Total primary (non-secondary, non-supplementary, mapped) records sampled.
    pub n_primary_sampled: usize,
    /// Primary records carrying the `0x400` (PCR/optical duplicate) flag.
    pub n_duplicate_flagged: usize,
}

/// Run a single first-pass scan and return paired-end + duplicate-mark presence.
///
/// Stops after `sample_size` primary records, or end-of-file. See
/// [`FirstPassSample`] for the returned shape.
pub fn detect_first_pass(bam_path: &str, sample_size: usize) -> Result<FirstPassSample> {
    use rust_htslib::bam::record::Record;

    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("Cannot open BAM '{}' for first-pass auto-detect", bam_path))?;

    let mut record = Record::new();
    let mut n_paired: usize = 0;
    let mut n_dup: usize = 0;
    let mut n_total: usize = 0;
    while n_total < sample_size {
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            Some(Err(e)) => {
                return Err(e).with_context(|| {
                    format!("Failed reading BAM record during first-pass auto-detect: {bam_path}")
                });
            }
            None => break,
        }
        if record.is_secondary() || record.is_supplementary() || record.is_unmapped() {
            continue;
        }
        if record.is_paired() {
            n_paired += 1;
        }
        if record.is_duplicate() {
            n_dup += 1;
        }
        n_total += 1;
    }
    Ok(FirstPassSample {
        paired: n_total > 0 && n_paired * 2 > n_total,
        n_primary_sampled: n_total,
        n_duplicate_flagged: n_dup,
    })
}

/// Decision returned by [`detect_chr_bridge`] for the BAM↔GTF chromosome
/// naming comparison.
#[derive(Debug, PartialEq, Eq)]
pub enum ChrBridge {
    /// BAM and GTF naming agree, or the user already specified an explicit
    /// `chromosome_prefix`/`chromosome_mapping` and we should not override.
    /// Either way, no auto-bridge is needed.
    None,
    /// BAM is no-prefix, GTF is `chr`-prefixed → auto-apply `chr` as a
    /// prefix at runtime so BAM names match GTF names. Mirrors what the
    /// user would type as `chromosome_prefix: "chr"` in the YAML config.
    AddChrPrefix,
    /// BAM is `chr`-prefixed, GTF is no-prefix. Not auto-bridgeable with
    /// the existing prefix-only translation in counting; surface a friendly
    /// bail with a [`crate::config::RnaConfig::chromosome_mapping`] example.
    StripChrUnsupported,
}

/// Decide whether to auto-bridge a BAM/GTF chromosome naming mismatch.
///
/// User-supplied `chromosome_prefix` or `chromosome_mapping` always wins
/// (`user_explicit == true` short-circuits to [`ChrBridge::None`]) so this
/// auto-decision can never silently override an explicit config.
pub fn detect_chr_bridge(
    bam_style: ChromStyle,
    gtf_style: ChromStyle,
    user_explicit: bool,
) -> ChrBridge {
    if user_explicit {
        return ChrBridge::None;
    }
    match (bam_style, gtf_style) {
        (ChromStyle::NoChr, ChromStyle::Chr) => ChrBridge::AddChrPrefix,
        (ChromStyle::Chr, ChromStyle::NoChr) => ChrBridge::StripChrUnsupported,
        _ => ChrBridge::None,
    }
}

/// Run a pre-pass over a small sample of records and return the inferred
/// strandedness, or `None` when the sample is too sparse to reach a verdict.
///
/// `chrom_prefix` and `chrom_mapping` mirror the runtime translation done
/// in the main counting loop (so the gene model — keyed on GTF chrom names
/// — is queryable from BAM tids that may differ in naming style).
pub fn detect_strandedness(
    bam_path: &str,
    gene_model: &GeneModel,
    mapq_cut: u8,
    target_overlapping: u64,
    chrom_prefix: Option<&str>,
    chrom_mapping: &HashMap<String, String>,
) -> Result<Option<Strandedness>> {
    use rust_htslib::bam::record::Record;

    let mut reader = Reader::from_path(bam_path).with_context(|| {
        format!(
            "Cannot open BAM '{}' for strandedness auto-detect",
            bam_path
        )
    })?;
    let header = reader.header().clone();
    // Pre-cache tid → chrom name so the per-record loop avoids allocating
    // a fresh String on every iteration (up to 1M iterations). Translate
    // BAM names into GTF names so gene-model lookups succeed regardless of
    // chr-prefix style.
    let tid_to_bam_name: Vec<String> = (0..header.target_count())
        .map(|t| String::from_utf8_lossy(header.tid2name(t)).into_owned())
        .collect();
    let tid_to_name: Vec<String> =
        crate::genome::resolve_gtf_chrom_names(&tid_to_bam_name, chrom_prefix, chrom_mapping);

    let mut accum = InferExpAccum::default();
    let mut record = Record::new();
    let mut n_read: u64 = 0;
    while n_read < STRANDED_HARD_CAP {
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            Some(Err(e)) => {
                return Err(e).with_context(|| {
                    format!("Failed reading BAM record during strandedness auto-detect: {bam_path}")
                });
            }
            None => break,
        }
        n_read += 1;

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let chrom = match tid_to_name.get(tid as usize) {
            Some(s) => s.as_str(),
            None => continue,
        };
        accum.process_read(&record, chrom, gene_model, mapq_cut);

        if n_read.is_multiple_of(STRANDED_EXIT_CHECK_EVERY) {
            let p: u64 = accum.p_strandness.values().sum();
            let s: u64 = accum.s_strandness.values().sum();
            if p + s >= target_overlapping {
                break;
            }
        }
    }
    Ok(match infer_strandedness(&accum.into_result()) {
        InferredStrandedness::Forward => Some(Strandedness::Forward),
        InferredStrandedness::Reverse => Some(Strandedness::Reverse),
        InferredStrandedness::Unstranded => Some(Strandedness::Unstranded),
        InferredStrandedness::Undetermined => None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::Read;
    use std::path::Path;

    fn write_synthetic_bam(path: &Path, sam_lines: &[&str]) {
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1000_i64),
        );
        let header_view = bam::HeaderView::from_header(&header);
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            for line in sam_lines {
                let record = bam::Record::from_sam(&header_view, line.as_bytes()).unwrap();
                writer.write(&record).unwrap();
            }
        }
    }

    #[test]
    fn first_pass_paired_true_for_paired_records() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("paired.bam");
        write_synthetic_bam(
            &bam,
            &[
                "r1\t99\tchr1\t1\t60\t4M\t=\t10\t13\tACGT\tIIII",
                "r2\t99\tchr1\t2\t60\t4M\t=\t11\t13\tACGT\tIIII",
                "r3\t99\tchr1\t3\t60\t4M\t=\t12\t13\tACGT\tIIII",
            ],
        );
        assert!(
            detect_first_pass(bam.to_str().unwrap(), 100)
                .unwrap()
                .paired
        );
    }

    #[test]
    fn first_pass_paired_false_for_single_end_records() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("single.bam");
        write_synthetic_bam(
            &bam,
            &[
                "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII",
            ],
        );
        assert!(
            !detect_first_pass(bam.to_str().unwrap(), 100)
                .unwrap()
                .paired
        );
    }

    #[test]
    fn first_pass_paired_false_for_empty_bam() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("empty.bam");
        write_synthetic_bam(&bam, &[]);
        assert!(
            !detect_first_pass(bam.to_str().unwrap(), 100)
                .unwrap()
                .paired
        );
    }

    #[test]
    fn first_pass_skips_secondary_and_supplementary() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("sec.bam");
        write_synthetic_bam(
            &bam,
            &[
                "r1\t257\tchr1\t1\t60\t4M\t=\t10\t13\tACGT\tIIII",
                "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII",
            ],
        );
        let mut reader = bam::Reader::from_path(&bam).unwrap();
        let mut record = bam::Record::new();
        let mut count = 0;
        while let Some(Ok(())) = reader.read(&mut record) {
            count += 1;
        }
        assert_eq!(count, 2);
        assert!(
            !detect_first_pass(bam.to_str().unwrap(), 100)
                .unwrap()
                .paired
        );
    }

    #[test]
    fn first_pass_counts_duplicate_flagged_reads() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("dup.bam");
        // Flag 0x400 (1024) is the duplicate bit. Mix 2 dup-flagged + 1 plain.
        write_synthetic_bam(
            &bam,
            &[
                "r1\t1024\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "r2\t1024\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "r3\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII",
            ],
        );
        let fp = detect_first_pass(bam.to_str().unwrap(), 100).unwrap();
        assert_eq!(fp.n_primary_sampled, 3);
        assert_eq!(fp.n_duplicate_flagged, 2);
    }

    #[test]
    fn first_pass_zero_dup_when_unmarked() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("nodup.bam");
        write_synthetic_bam(
            &bam,
            &[
                "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII",
            ],
        );
        let fp = detect_first_pass(bam.to_str().unwrap(), 100).unwrap();
        assert_eq!(fp.n_primary_sampled, 2);
        assert_eq!(fp.n_duplicate_flagged, 0);
    }

    #[test]
    fn chr_bridge_no_op_when_user_explicit() {
        // User-supplied prefix or mapping always wins, even on a real mismatch.
        assert_eq!(
            detect_chr_bridge(ChromStyle::NoChr, ChromStyle::Chr, true),
            ChrBridge::None
        );
    }

    #[test]
    fn chr_bridge_adds_chr_for_nochr_bam_chr_gtf() {
        assert_eq!(
            detect_chr_bridge(ChromStyle::NoChr, ChromStyle::Chr, false),
            ChrBridge::AddChrPrefix
        );
    }

    #[test]
    fn chr_bridge_unsupported_for_chr_bam_nochr_gtf() {
        assert_eq!(
            detect_chr_bridge(ChromStyle::Chr, ChromStyle::NoChr, false),
            ChrBridge::StripChrUnsupported
        );
    }

    #[test]
    fn chr_bridge_no_op_when_styles_match() {
        assert_eq!(
            detect_chr_bridge(ChromStyle::Chr, ChromStyle::Chr, false),
            ChrBridge::None
        );
        assert_eq!(
            detect_chr_bridge(ChromStyle::NoChr, ChromStyle::NoChr, false),
            ChrBridge::None
        );
    }

    #[test]
    fn chr_bridge_no_op_when_unknown() {
        assert_eq!(
            detect_chr_bridge(ChromStyle::Unknown, ChromStyle::Chr, false),
            ChrBridge::None
        );
        assert_eq!(
            detect_chr_bridge(ChromStyle::Chr, ChromStyle::Unknown, false),
            ChrBridge::None
        );
    }
}
