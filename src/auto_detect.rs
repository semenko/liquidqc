//! Pre-pass auto-detection for library properties that liquidqc would
//! otherwise force the user to declare on the CLI.
//!
//! Both helpers open the BAM independently of the main counting pass and
//! read at most a small bounded sample of records. The cost is sub-second
//! on indexed BAMs and serves to make `liquidqc rna sample.bam` Just Work
//! without `--paired` or `--stranded` for typical inputs.

use crate::cli::Strandedness;
use crate::rna::rseqc::accumulators::InferExpAccum;
use crate::rna::rseqc::infer_experiment::{infer_strandedness, GeneModel, InferredStrandedness};
use anyhow::{Context, Result};
use rust_htslib::bam::{Read as BamRead, Reader};

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

/// Auto-detect paired-end based on a sample of primary records.
///
/// Returns `true` if more than half of the inspected primary records have
/// the `BAM_FPAIRED` flag set. Falls back to `false` (single-end) when
/// the BAM is empty or every primary record is filtered.
pub fn detect_paired(bam_path: &str, sample_size: usize) -> Result<bool> {
    use rust_htslib::bam::record::Record;

    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("Cannot open BAM '{}' for paired auto-detect", bam_path))?;

    let mut record = Record::new();
    let mut n_paired: usize = 0;
    let mut n_total: usize = 0;
    while n_total < sample_size {
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            Some(Err(e)) => {
                return Err(e).with_context(|| {
                    format!("Failed reading BAM record during paired auto-detect: {bam_path}")
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
        n_total += 1;
    }
    if n_total == 0 {
        return Ok(false);
    }
    Ok(n_paired * 2 > n_total)
}

/// Run a pre-pass over a small sample of records and return the inferred
/// strandedness, or `None` when the sample is too sparse to reach a verdict.
pub fn detect_strandedness(
    bam_path: &str,
    gene_model: &GeneModel,
    mapq_cut: u8,
    target_overlapping: u64,
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
    // a fresh String on every iteration (up to 1M iterations).
    let tid_to_name: Vec<String> = (0..header.target_count())
        .map(|t| String::from_utf8_lossy(header.tid2name(t)).into_owned())
        .collect();

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
    fn detect_paired_returns_true_for_paired_records() {
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
        assert!(detect_paired(bam.to_str().unwrap(), 100).unwrap());
    }

    #[test]
    fn detect_paired_returns_false_for_single_end_records() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("single.bam");
        write_synthetic_bam(
            &bam,
            &[
                "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "r2\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII",
            ],
        );
        assert!(!detect_paired(bam.to_str().unwrap(), 100).unwrap());
    }

    #[test]
    fn detect_paired_returns_false_for_empty_bam() {
        let tmp = tempfile::tempdir().unwrap();
        let bam = tmp.path().join("empty.bam");
        write_synthetic_bam(&bam, &[]);
        assert!(!detect_paired(bam.to_str().unwrap(), 100).unwrap());
    }

    #[test]
    fn detect_paired_skips_secondary_and_supplementary() {
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
        assert!(!detect_paired(bam.to_str().unwrap(), 100).unwrap());
    }
}
