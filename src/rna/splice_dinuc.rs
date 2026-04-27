//! Splice-site dinucleotide composition of detected junctions.
//!
//! Tally each unique junction's donor/acceptor 2-mer into `gt_ag` (canonical
//! U2, also `CT`/`AC` reverse-complemented), `gc_ag`, `at_ac`, or `other`.
//! Unique-junction tally — depth bias would otherwise dominate.

use anyhow::{Context, Result};
use rust_htslib::faidx;
use serde::Serialize;
use std::path::Path;

use crate::rna::rseqc::junction_annotation::{Junction, JunctionResults};

/// Canonical four-class result.
#[derive(Debug, Default, Clone, Copy, Serialize)]
pub struct SpliceSiteDinucleotidesResult {
    pub gt_ag: u64,
    pub gc_ag: u64,
    pub at_ac: u64,
    pub other: u64,
    /// Junctions skipped because the donor or acceptor 2-mer ran off the
    /// contig or the contig is missing from the FASTA.
    pub skipped_oob: u64,
    /// Total unique junctions considered (sum of the four classes plus
    /// `skipped_oob`).
    pub junctions_total: u64,
}

/// Compute the dinucleotide tally over the unique junctions in `results`.
///
/// One [`faidx::Reader`] is opened for the duration of the call. Errors only
/// when the FASTA index is unreadable; per-junction lookup failures are
/// counted in `skipped_oob` and never abort the run.
pub fn compute(
    results: &JunctionResults,
    fasta_path: &Path,
) -> Result<SpliceSiteDinucleotidesResult> {
    let reader = faidx::Reader::from_path(fasta_path).with_context(|| {
        format!(
            "Failed to open FASTA for splice-site dinucleotides: {}",
            fasta_path.display()
        )
    })?;

    let mut out = SpliceSiteDinucleotidesResult::default();
    for junction in results.junctions.keys() {
        out.junctions_total += 1;
        match classify(&reader, junction) {
            Some(SpliceClass::GtAg) => out.gt_ag += 1,
            Some(SpliceClass::GcAg) => out.gc_ag += 1,
            Some(SpliceClass::AtAc) => out.at_ac += 1,
            Some(SpliceClass::Other) => out.other += 1,
            None => out.skipped_oob += 1,
        }
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SpliceClass {
    GtAg,
    GcAg,
    AtAc,
    Other,
}

fn classify(reader: &faidx::Reader, junction: &Junction) -> Option<SpliceClass> {
    // RSeQC stores chrom uppercased internally; FASTA names follow the BAM
    // header. Try the uppercased chrom first, then the lowercased "chr…"
    // variant typical of UCSC-style references — both forms are common in
    // practice. (The aligner-generated BAM and the FASTA always agree, so
    // exactly one of the two will succeed.)
    let chrom = junction.chrom.as_str();
    let donor = fetch_2mer(reader, chrom, junction.intron_start)
        .or_else(|| fetch_2mer(reader, &chrom_lower(chrom), junction.intron_start))?;
    // `intron_end` is exclusive (first base after the intron); donor is at
    // `intron_start..intron_start+2`, acceptor is at `intron_end-2..intron_end`.
    let acceptor_start = junction.intron_end.checked_sub(2)?;
    let acceptor = fetch_2mer(reader, chrom, acceptor_start)
        .or_else(|| fetch_2mer(reader, &chrom_lower(chrom), acceptor_start))?;
    Some(classify_pair(donor, acceptor))
}

fn fetch_2mer(reader: &faidx::Reader, chrom: &str, pos: u64) -> Option<[u8; 2]> {
    if reader.fetch_seq_len(chrom) == 0 {
        return None;
    }
    let bases = reader
        .fetch_seq(chrom, pos as usize, (pos + 1) as usize)
        .ok()?;
    if bases.len() != 2 {
        return None;
    }
    Some([bases[0].to_ascii_uppercase(), bases[1].to_ascii_uppercase()])
}

fn chrom_lower(chrom: &str) -> String {
    // RSeQC junction tracking uppercases the chrom name; map back to the
    // standard lowercase "chr…" form for FASTA fallback.
    chrom.replace("CHR", "chr")
}

fn classify_pair(donor: [u8; 2], acceptor: [u8; 2]) -> SpliceClass {
    // Reference-orientation pairs for the three biological classes — both
    // strands listed since the GTF doesn't tell us the intron's strand.
    const GT: [u8; 2] = [b'G', b'T'];
    const AG: [u8; 2] = [b'A', b'G'];
    const CT: [u8; 2] = [b'C', b'T'];
    const AC: [u8; 2] = [b'A', b'C'];
    const GC: [u8; 2] = [b'G', b'C'];
    const AT: [u8; 2] = [b'A', b'T'];
    match (donor, acceptor) {
        // GT-AG canonical (+) strand, or CT-AC reverse-complemented from (-).
        (GT, AG) | (CT, AC) => SpliceClass::GtAg,
        // GC-AG alternative (+) strand, or CT-GC from (-).
        (GC, AG) | (CT, GC) => SpliceClass::GcAg,
        // AT-AC U12-type (+) strand, or GT-AT from (-).
        (AT, AC) | (GT, AT) => SpliceClass::AtAc,
        _ => SpliceClass::Other,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn d(b: &[u8; 2]) -> [u8; 2] {
        *b
    }

    #[test]
    fn gt_ag_canonical() {
        assert_eq!(classify_pair(d(b"GT"), d(b"AG")), SpliceClass::GtAg);
    }

    #[test]
    fn gt_ag_minus_strand_form() {
        assert_eq!(classify_pair(d(b"CT"), d(b"AC")), SpliceClass::GtAg);
    }

    #[test]
    fn gc_ag_classes_match() {
        assert_eq!(classify_pair(d(b"GC"), d(b"AG")), SpliceClass::GcAg);
        assert_eq!(classify_pair(d(b"CT"), d(b"GC")), SpliceClass::GcAg);
    }

    #[test]
    fn at_ac_classes_match() {
        assert_eq!(classify_pair(d(b"AT"), d(b"AC")), SpliceClass::AtAc);
        assert_eq!(classify_pair(d(b"GT"), d(b"AT")), SpliceClass::AtAc);
    }

    #[test]
    fn other_falls_through() {
        assert_eq!(classify_pair(d(b"NN"), d(b"AG")), SpliceClass::Other);
        assert_eq!(classify_pair(d(b"GG"), d(b"GG")), SpliceClass::Other);
    }

    #[test]
    fn chrom_lower_strips_uppercase_chr_prefix() {
        assert_eq!(chrom_lower("CHR1"), "chr1");
        assert_eq!(chrom_lower("CHRX"), "chrX");
        assert_eq!(chrom_lower("1"), "1"); // numeric chroms unchanged
    }
}
