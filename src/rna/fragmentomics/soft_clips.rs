//! Phase 2 Tier 1: soft-clip 4-mers from read sequences.
//!
//! For each primary mapped read we look at the 5' and 3' ends *in sequencing
//! direction* (orientation-corrected against the read strand). A leading
//! CIGAR `S` on a forward-strand read is a 5'-clip; on a reverse-strand
//! read it's a 3'-clip. The 4-mer recorded for each clipped end is the
//! 4 read bases at the read tip in sequencing direction:
//!
//! - 5'-clip 4-mer: `read[0..4]` — the very first 4 bases sequenced. When
//!   the clip length is < 4, only the first `clip_len` bases are from the
//!   clip; the rest are aligned bases adjacent to the cut site.
//! - 3'-clip 4-mer: `read[N-4..N]` — the last 4 bases sequenced.
//!
//! The 4-mer is binned by clip length on that end: 1 nt, 2 nt, or ≥3 nt.
//! Per-end clip rates (`soft_clip_rate_5p`, `soft_clip_rate_3p`) are over
//! all primary mapped reads passing the MAPQ filter — including those with
//! no clip, which contribute to the denominator but not to any k-mer bin.
//!
//! Reads shorter than 4 bases of resolved sequence are skipped from k-mer
//! counting (they still count toward the denominator).

use indexmap::IndexMap;
use rust_htslib::bam::record::Record;
use serde::Serialize;

use super::common::{
    encode_kmer4, is_primary_mapped, is_reverse_strand, iupac_nibble_to_acgt, kmer_array_to_map,
    leading_softclip_len, rc_kmer4, trailing_softclip_len,
};

/// Number of clip-length bins: 1 nt, 2 nt, ≥3 nt.
pub const N_CLIP_LEN_BINS: usize = 3;

/// Index assignment for the three clip-length bins.
fn clip_len_bin(clip_len: usize) -> Option<usize> {
    match clip_len {
        0 => None,
        1 => Some(0),
        2 => Some(1),
        _ => Some(2),
    }
}

#[derive(Debug, Clone)]
pub struct SoftClipAccum {
    /// Per-clip-length bin: 256-entry 4-mer counts at the read's 5' end.
    pub(crate) kmer_5p: [[u64; 256]; N_CLIP_LEN_BINS],
    /// Same for the read's 3' end.
    pub(crate) kmer_3p: [[u64; 256]; N_CLIP_LEN_BINS],
    /// Reads with any clip on the 5' end (≥ 1 nt). Numerator for `soft_clip_rate_5p`.
    pub reads_with_5p_clip: u64,
    /// Reads with any clip on the 3' end (≥ 1 nt). Numerator for `soft_clip_rate_3p`.
    pub reads_with_3p_clip: u64,
    /// Total primary mapped reads observed (both denominators).
    pub primary_reads_observed: u64,
}

impl Default for SoftClipAccum {
    fn default() -> Self {
        Self::new()
    }
}

impl SoftClipAccum {
    pub fn new() -> Self {
        Self {
            kmer_5p: [[0u64; 256]; N_CLIP_LEN_BINS],
            kmer_3p: [[0u64; 256]; N_CLIP_LEN_BINS],
            reads_with_5p_clip: 0,
            reads_with_3p_clip: 0,
            primary_reads_observed: 0,
        }
    }

    pub fn process_read(&mut self, record: &Record, mapq_cut: u8) {
        if !is_primary_mapped(record, mapq_cut) {
            return;
        }
        self.primary_reads_observed += 1;

        let cigar = record.cigar();
        let leading = leading_softclip_len(&cigar);
        let trailing = trailing_softclip_len(&cigar);
        let reverse = is_reverse_strand(record);

        // Orientation-correct: which end of the *read* is the leading clip on?
        let (clip_5p_len, clip_3p_len) = if reverse {
            (trailing, leading)
        } else {
            (leading, trailing)
        };

        if clip_5p_len > 0 {
            self.reads_with_5p_clip += 1;
        }
        if clip_3p_len > 0 {
            self.reads_with_3p_clip += 1;
        }

        let seq_len = record.seq_len();
        if seq_len < 4 {
            return; // Cannot form a 4-mer.
        }

        // SEQ stores the forward-strand sequence; for reverse-strand reads
        // we reverse-complement BAM-stored bases to recover sequencing direction.
        // 5' end of read = SEQ[0..4] (forward) or revcomp(SEQ[N-4..N]) (reverse).
        // 3' end of read = SEQ[N-4..N] (forward) or revcomp(SEQ[0..4]) (reverse).
        let seq = record.seq();
        let n_minus_4 = seq_len - 4;
        let (start_5p, start_3p) = if reverse {
            (n_minus_4, 0)
        } else {
            (0, n_minus_4)
        };
        if let Some(bin) = clip_len_bin(clip_5p_len) {
            if let Some(kmer) = read_kmer_at(&seq, start_5p, reverse) {
                self.kmer_5p[bin][kmer as usize] += 1;
            }
        }
        if let Some(bin) = clip_len_bin(clip_3p_len) {
            if let Some(kmer) = read_kmer_at(&seq, start_3p, reverse) {
                self.kmer_3p[bin][kmer as usize] += 1;
            }
        }
    }

    pub fn merge(&mut self, other: SoftClipAccum) {
        for bin in 0..N_CLIP_LEN_BINS {
            for k in 0..256 {
                self.kmer_5p[bin][k] += other.kmer_5p[bin][k];
                self.kmer_3p[bin][k] += other.kmer_3p[bin][k];
            }
        }
        self.reads_with_5p_clip += other.reads_with_5p_clip;
        self.reads_with_3p_clip += other.reads_with_3p_clip;
        self.primary_reads_observed += other.primary_reads_observed;
    }

    pub fn into_result(self) -> SoftClipResult {
        let total = self.primary_reads_observed.max(1);
        let rate_5p = if self.primary_reads_observed == 0 {
            0.0
        } else {
            self.reads_with_5p_clip as f64 / total as f64
        };
        let rate_3p = if self.primary_reads_observed == 0 {
            0.0
        } else {
            self.reads_with_3p_clip as f64 / total as f64
        };

        SoftClipResult {
            soft_clip_rate_5p: rate_5p,
            soft_clip_rate_3p: rate_3p,
            kmer_5p_by_clip_len: kmer_arrays_to_block(&self.kmer_5p),
            kmer_3p_by_clip_len: kmer_arrays_to_block(&self.kmer_3p),
            primary_reads_observed: self.primary_reads_observed,
        }
    }
}

fn kmer_arrays_to_block(arrays: &[[u64; 256]; N_CLIP_LEN_BINS]) -> SoftClipKmerByLen {
    SoftClipKmerByLen {
        len_1: kmer_array_to_map(&arrays[0]),
        len_2: kmer_array_to_map(&arrays[1]),
        len_3plus: kmer_array_to_map(&arrays[2]),
    }
}

/// Encode a 4-mer from BAM `seq[start..start+4]`, reverse-complementing
/// when `rc` is true. Returns `None` on any non-ACGT base or if the
/// window runs off the end of `seq`.
fn read_kmer_at(seq: &rust_htslib::bam::record::Seq<'_>, start: usize, rc: bool) -> Option<u8> {
    if start + 4 > seq.len() {
        return None;
    }
    let mut bases = [0u8; 4];
    for (i, b) in bases.iter_mut().enumerate() {
        *b = iupac_nibble_to_acgt(seq.encoded_base(start + i))?;
    }
    let k = encode_kmer4(&bases)?;
    Some(if rc { rc_kmer4(k) } else { k })
}

// ---------------------------------------------------------------------------
// Serializable output
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize)]
pub struct SoftClipResult {
    pub soft_clip_rate_5p: f64,
    pub soft_clip_rate_3p: f64,
    pub kmer_5p_by_clip_len: SoftClipKmerByLen,
    pub kmer_3p_by_clip_len: SoftClipKmerByLen,
    pub primary_reads_observed: u64,
}

#[derive(Debug, Serialize)]
pub struct SoftClipKmerByLen {
    pub len_1: IndexMap<String, u64>,
    pub len_2: IndexMap<String, u64>,
    pub len_3plus: IndexMap<String, u64>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clip_len_bin_assignment() {
        assert_eq!(clip_len_bin(0), None);
        assert_eq!(clip_len_bin(1), Some(0));
        assert_eq!(clip_len_bin(2), Some(1));
        assert_eq!(clip_len_bin(3), Some(2));
        assert_eq!(clip_len_bin(4), Some(2));
        assert_eq!(clip_len_bin(99), Some(2));
    }

    #[test]
    fn merge_is_additive() {
        let mut a = SoftClipAccum::new();
        a.kmer_5p[0][42] = 3;
        a.reads_with_5p_clip = 5;
        a.primary_reads_observed = 10;
        let mut b = SoftClipAccum::new();
        b.kmer_5p[0][42] = 2;
        b.reads_with_5p_clip = 1;
        b.primary_reads_observed = 7;
        a.merge(b);
        assert_eq!(a.kmer_5p[0][42], 5);
        assert_eq!(a.reads_with_5p_clip, 6);
        assert_eq!(a.primary_reads_observed, 17);
    }

    #[test]
    fn rates_are_zero_when_no_reads_observed() {
        let r = SoftClipAccum::new().into_result();
        assert_eq!(r.soft_clip_rate_5p, 0.0);
        assert_eq!(r.soft_clip_rate_3p, 0.0);
        assert_eq!(r.primary_reads_observed, 0);
    }

    #[test]
    fn block_drops_zero_keys() {
        let mut a = SoftClipAccum::new();
        a.primary_reads_observed = 1;
        a.reads_with_5p_clip = 1;
        a.kmer_5p[2][0] = 7; // AAAA bin
        let r = a.into_result();
        let m = &r.kmer_5p_by_clip_len.len_3plus;
        assert_eq!(m.len(), 1);
        assert_eq!(m.get("AAAA"), Some(&7));
        assert!(r.kmer_5p_by_clip_len.len_1.is_empty());
        assert!(r.kmer_5p_by_clip_len.len_2.is_empty());
    }
}
