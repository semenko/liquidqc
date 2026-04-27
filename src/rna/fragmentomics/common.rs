//! Shared helpers across the four Phase 2 fragmentomics accumulators:
//! per-record proper-pair / leftmost-mate filtering, 4-mer encoding /
//! decoding, IUPAC nibble decoding, and reverse-complement.
//!
//! Per-pair dedup convention (used by [`fragment_size`](super::fragment_size),
//! [`end_motifs`](super::end_motifs)): only the leftmost mate of a properly
//! paired alignment contributes one observation. Selected as the record
//! whose `insert_size()` is strictly positive (the SAM TLEN sign convention
//! makes the leftmost mate's TLEN positive). This is bitwise equivalent to
//! the `pos < mpos || (pos == mpos && is_first_in_template)` dedup used by
//! the inherited [`InnerDistAccum`](crate::rna::rseqc::accumulators::InnerDistAccum)
//! for properly paired records and avoids the chromosome-spanning logic
//! since cross-chromosome pairs do not have a defined fragment length.

use indexmap::IndexMap;
use rust_htslib::bam::record::{Cigar, Record};

use crate::rna::bam_flags::{
    BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREVERSE,
    BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
};

/// True iff `record` is the leftmost mate of a properly paired alignment
/// suitable for one fragmentomics observation per pair.
///
/// Filters: not unmapped, mate not unmapped, primary alignment, not
/// duplicate / QC-fail / secondary / supplementary, paired + proper pair,
/// same chromosome as mate, MAPQ ≥ `mapq_cut`, `insert_size() > 0`.
pub fn is_leftmost_proper_pair_mate(record: &Record, mapq_cut: u8) -> bool {
    let flags = record.flags();
    if flags & BAM_FUNMAP != 0
        || flags & BAM_FMUNMAP != 0
        || flags & BAM_FSECONDARY != 0
        || flags & BAM_FSUPPLEMENTARY != 0
        || flags & BAM_FQCFAIL != 0
        || flags & BAM_FDUP != 0
    {
        return false;
    }
    if flags & BAM_FPAIRED == 0 || flags & BAM_FPROPER_PAIR == 0 {
        return false;
    }
    if record.tid() != record.mtid() {
        return false;
    }
    if record.mapq() < mapq_cut {
        return false;
    }
    record.insert_size() > 0
}

/// True iff `record` is a primary, mapped alignment suitable for soft-clip
/// k-mer counting (does not require pairing).
pub fn is_primary_mapped(record: &Record, mapq_cut: u8) -> bool {
    let flags = record.flags();
    if flags & BAM_FUNMAP != 0
        || flags & BAM_FSECONDARY != 0
        || flags & BAM_FSUPPLEMENTARY != 0
        || flags & BAM_FQCFAIL != 0
        || flags & BAM_FDUP != 0
    {
        return false;
    }
    record.mapq() >= mapq_cut
}

/// True iff this primary record is on the reverse strand.
pub fn is_reverse_strand(record: &Record) -> bool {
    record.flags() & BAM_FREVERSE != 0
}

/// Length in nt of a leading CIGAR `S` op, or 0 if the alignment doesn't
/// start with one. SAM/BAM stores CIGAR in reference orientation, so this
/// is the soft clip at the leftmost reference position regardless of
/// read strand.
pub fn leading_softclip_len(cigar: &rust_htslib::bam::record::CigarStringView) -> usize {
    match cigar.iter().next() {
        Some(Cigar::SoftClip(n)) => *n as usize,
        _ => 0,
    }
}

/// Length in nt of a trailing CIGAR `S` op, or 0.
pub fn trailing_softclip_len(cigar: &rust_htslib::bam::record::CigarStringView) -> usize {
    match cigar.iter().next_back() {
        Some(Cigar::SoftClip(n)) => *n as usize,
        _ => 0,
    }
}

/// Strand-corrected soft-clip lengths at the read's 5' and 3' ends.
/// For forward-strand reads the leading CIGAR `S` is the 5' clip; for
/// reverse-strand reads the leading CIGAR `S` is the read's 3' clip.
pub fn softclip_lens_5p3p(record: &Record) -> (usize, usize) {
    let cigar = record.cigar();
    let leading = leading_softclip_len(&cigar);
    let trailing = trailing_softclip_len(&cigar);
    if is_reverse_strand(record) {
        (trailing, leading)
    } else {
        (leading, trailing)
    }
}

/// Strand-corrected presence of soft clips at the read's 5' and 3' ends.
pub fn softclip_status_5p3p(record: &Record) -> (bool, bool) {
    let (a, b) = softclip_lens_5p3p(record);
    (a > 0, b > 0)
}

// ---------------------------------------------------------------------------
// 4-mer encoding (A=0, C=1, G=2, T=3) packed into a u8 (2 bits per base).
// ---------------------------------------------------------------------------

/// Encode a 4-base ACGT slice into a u8 (`AAAA`=0, `TTTT`=255). Returns
/// `None` if the slice contains any non-ACGT base (case-insensitive).
pub fn encode_kmer4(bases: &[u8]) -> Option<u8> {
    debug_assert_eq!(bases.len(), 4);
    let mut k: u8 = 0;
    for b in bases {
        let two_bit = match *b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        k = (k << 2) | two_bit;
    }
    Some(k)
}

/// Decode a 4-mer back to a 4-byte ACGT string.
pub fn decode_kmer4(k: u8) -> [u8; 4] {
    let mut out = [0u8; 4];
    for (i, b) in out.iter_mut().enumerate() {
        let shift = (3 - i) * 2;
        let two_bit = (k >> shift) & 0b11;
        *b = match two_bit {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
    }
    out
}

/// Reverse-complement a 4-mer encoded by [`encode_kmer4`].
pub fn rc_kmer4(k: u8) -> u8 {
    let mut out: u8 = 0;
    for i in 0..4 {
        let shift = i * 2;
        let two_bit = (k >> shift) & 0b11;
        // Complement: A<->T (0<->3), C<->G (1<->2). XOR with 0b11 swaps both.
        let comp = two_bit ^ 0b11;
        out = (out << 2) | comp;
    }
    out
}

// ---------------------------------------------------------------------------
// IUPAC nibble decoding (BAM-encoded 4-bit base codes).
// htslib uses: =0 A1 C2 M3 G4 R5 S6 V7 T8 W9 Y10 H11 K12 D13 B14 N15.
// We only accept ACGT; other codes return None.
// ---------------------------------------------------------------------------

/// Decode a single 4-bit IUPAC code to an ACGT ASCII byte. Returns `None`
/// for ambiguous / non-ACGT codes (which fragmentomics callers drop).
pub fn iupac_nibble_to_acgt(code: u8) -> Option<u8> {
    match code {
        1 => Some(b'A'),
        2 => Some(b'C'),
        4 => Some(b'G'),
        8 => Some(b'T'),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// 4-mer histogram → IndexMap conversion (shared finalization helper).
// ---------------------------------------------------------------------------

/// Trait for the integer count types stored in 256-entry 4-mer histograms.
/// Lets the entropy / JSD / serialization helpers below accept either
/// `[u64; 256]` (Tier-1, sample-level) or `[u32; 256]` (Tier-2, per-gene).
pub trait KmerCount: Copy {
    fn as_u64(self) -> u64;
    fn as_f64(self) -> f64;
    fn is_zero(self) -> bool;
}
impl KmerCount for u64 {
    fn as_u64(self) -> u64 {
        self
    }
    fn as_f64(self) -> f64 {
        self as f64
    }
    fn is_zero(self) -> bool {
        self == 0
    }
}
impl KmerCount for u32 {
    fn as_u64(self) -> u64 {
        self as u64
    }
    fn as_f64(self) -> f64 {
        self as f64
    }
    fn is_zero(self) -> bool {
        self == 0
    }
}

/// Walk a 256-entry 4-mer count array and emit only the nonzero entries
/// as `(decoded_kmer_string, count)` pairs in 4-mer-index order. Used by
/// both `IndexMap<String, _>` and `Vec<(String, _)>` consumers.
pub fn kmer_array_iter_nonzero<T: KmerCount>(
    arr: &[T; 256],
) -> impl Iterator<Item = (String, T)> + '_ {
    arr.iter()
        .enumerate()
        .filter_map(|(k, c)| if c.is_zero() { None } else { Some((k, *c)) })
        .map(|(k, c)| {
            let bytes = decode_kmer4(k as u8);
            let key = String::from_utf8(bytes.to_vec()).expect("ACGT bytes are valid UTF-8");
            (key, c)
        })
}

/// Convert a 256-entry 4-mer count array into an `IndexMap<String, T>`
/// suitable for JSON serialization, dropping zero entries.
pub fn kmer_array_to_map<T: KmerCount>(arr: &[T; 256]) -> IndexMap<String, T> {
    kmer_array_iter_nonzero(arr).collect()
}

/// Normalize a 256-entry count array to a probability distribution.
/// Returns `None` when the array sums to zero.
pub fn normalize_kmer_counts<T: KmerCount>(arr: &[T; 256]) -> Option<[f64; 256]> {
    let total: u64 = arr.iter().map(|c| c.as_u64()).sum();
    if total == 0 {
        return None;
    }
    let total_f = total as f64;
    let mut out = [0f64; 256];
    for (i, c) in arr.iter().enumerate() {
        out[i] = c.as_f64() / total_f;
    }
    Some(out)
}

/// Shannon entropy (base 2) of a 4-mer distribution. Returns `None` when
/// the array is all zeros (no observations).
pub fn shannon_entropy_log2<T: KmerCount>(arr: &[T; 256]) -> Option<f64> {
    shannon_entropy_log2_from_probs(normalize_kmer_counts(arr)?)
}

/// Jensen-Shannon divergence (base 2) between two 4-mer distributions.
/// Returns `None` when either input has zero total count.
pub fn jsd_log2<T: KmerCount>(a: &[T; 256], b: &[T; 256]) -> Option<f64> {
    let pa = normalize_kmer_counts(a)?;
    let pb = normalize_kmer_counts(b)?;
    Some(jsd_log2_from_probs(&pa, &pb))
}

/// Shannon entropy of a pre-normalized probability distribution.
/// Lets callers compute entropy + JSD without re-normalizing the inputs
/// (saves two `O(256)` sums and divides per gene at finalize time).
pub fn shannon_entropy_log2_from_probs(probs: [f64; 256]) -> Option<f64> {
    let mut h = 0.0;
    for p in probs {
        if p > 0.0 {
            h -= p * p.log2();
        }
    }
    Some(h)
}

/// Jensen-Shannon divergence of two pre-normalized probability distributions.
pub fn jsd_log2_from_probs(pa: &[f64; 256], pb: &[f64; 256]) -> f64 {
    let mut jsd = 0.0;
    for i in 0..256 {
        let p = pa[i];
        let q = pb[i];
        let m = 0.5 * (p + q);
        if p > 0.0 && m > 0.0 {
            jsd += 0.5 * p * (p / m).log2();
        }
        if q > 0.0 && m > 0.0 {
            jsd += 0.5 * q * (q / m).log2();
        }
    }
    // Numerical floor to clamp tiny negative values from FP noise.
    jsd.max(0.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_encode_round_trip() {
        for k in 0u8..=255 {
            let bases = decode_kmer4(k);
            assert_eq!(encode_kmer4(&bases), Some(k));
        }
    }

    #[test]
    fn kmer_known_values() {
        assert_eq!(encode_kmer4(b"AAAA"), Some(0));
        assert_eq!(encode_kmer4(b"AAAC"), Some(1));
        assert_eq!(encode_kmer4(b"TTTT"), Some(0xFF));
        assert_eq!(encode_kmer4(b"ACGT"), Some(0b00_01_10_11));
    }

    #[test]
    fn kmer_rejects_non_acgt() {
        assert_eq!(encode_kmer4(b"ACGN"), None);
        assert_eq!(encode_kmer4(b"acgX"), None);
    }

    #[test]
    fn kmer_case_insensitive() {
        assert_eq!(encode_kmer4(b"acgt"), encode_kmer4(b"ACGT"));
    }

    #[test]
    fn reverse_complement_examples() {
        // "ACGT" -> "ACGT" (palindrome)
        assert_eq!(
            rc_kmer4(encode_kmer4(b"ACGT").unwrap()),
            encode_kmer4(b"ACGT").unwrap()
        );
        // "AAAA" -> "TTTT"
        assert_eq!(
            rc_kmer4(encode_kmer4(b"AAAA").unwrap()),
            encode_kmer4(b"TTTT").unwrap()
        );
        // "ATCG" -> "CGAT"
        assert_eq!(
            rc_kmer4(encode_kmer4(b"ATCG").unwrap()),
            encode_kmer4(b"CGAT").unwrap()
        );
        // Involution.
        for k in 0u8..=255 {
            assert_eq!(rc_kmer4(rc_kmer4(k)), k);
        }
    }

    use rust_htslib::bam::record::CigarString;

    fn build_record(cigar: CigarString, reverse: bool) -> Record {
        // Construct a minimal record. Tests only exercise flags + CIGAR.
        let mut record = Record::new();
        let total: u32 = cigar.iter().map(|c| c.len()).sum();
        let bases = vec![b'A'; total.max(1) as usize];
        let quals = vec![20u8; bases.len()];
        record.set(b"r", Some(&cigar), &bases, &quals);
        let mut flags: u16 = 0;
        if reverse {
            flags |= BAM_FREVERSE;
        }
        record.set_flags(flags);
        record
    }

    fn cigar(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    #[test]
    fn softclip_5p_3p_forward_strand() {
        // Leading soft clip → 5' on a forward-strand read.
        let r = build_record(cigar(vec![Cigar::SoftClip(3), Cigar::Match(10)]), false);
        assert_eq!(softclip_status_5p3p(&r), (true, false));

        // Trailing soft clip → 3' on a forward-strand read.
        let r = build_record(cigar(vec![Cigar::Match(10), Cigar::SoftClip(3)]), false);
        assert_eq!(softclip_status_5p3p(&r), (false, true));
    }

    #[test]
    fn softclip_5p_3p_reverse_strand_swaps() {
        // Leading soft clip on a reverse read = read 3' end.
        let r = build_record(cigar(vec![Cigar::SoftClip(3), Cigar::Match(10)]), true);
        assert_eq!(softclip_status_5p3p(&r), (false, true));

        // Trailing soft clip on a reverse read = read 5' end.
        let r = build_record(cigar(vec![Cigar::Match(10), Cigar::SoftClip(3)]), true);
        assert_eq!(softclip_status_5p3p(&r), (true, false));
    }

    #[test]
    fn no_softclip_returns_false_false() {
        let r = build_record(cigar(vec![Cigar::Match(13)]), false);
        assert_eq!(softclip_status_5p3p(&r), (false, false));
    }

    #[test]
    fn iupac_decoder_passes_acgt_only() {
        assert_eq!(iupac_nibble_to_acgt(1), Some(b'A'));
        assert_eq!(iupac_nibble_to_acgt(2), Some(b'C'));
        assert_eq!(iupac_nibble_to_acgt(4), Some(b'G'));
        assert_eq!(iupac_nibble_to_acgt(8), Some(b'T'));
        // N=15, M=3, etc.
        assert_eq!(iupac_nibble_to_acgt(15), None);
        assert_eq!(iupac_nibble_to_acgt(3), None);
        assert_eq!(iupac_nibble_to_acgt(0), None);
    }
}
