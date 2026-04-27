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
use rust_htslib::bam::record::Record;

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

/// Convert a 256-entry 4-mer count array into an `IndexMap<String, u64>`
/// suitable for JSON serialization, dropping zero entries. Used by both
/// the end-motif and soft-clip accumulators at finalization.
pub fn kmer_array_to_map(arr: &[u64; 256]) -> IndexMap<String, u64> {
    let mut m = IndexMap::with_capacity(arr.iter().filter(|c| **c > 0).count());
    for (k, count) in arr.iter().enumerate() {
        if *count == 0 {
            continue;
        }
        let bytes = decode_kmer4(k as u8);
        let key = String::from_utf8(bytes.to_vec()).expect("ACGT bytes are valid UTF-8");
        m.insert(key, *count);
    }
    m
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
