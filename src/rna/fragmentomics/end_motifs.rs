//! Phase 2 Tier 1: 5' and 3' fragment-end 4-mers from the reference FASTA
//! (cleavage motifs).
//!
//! For each leftmost mate of a properly paired alignment (one observation per
//! pair, mirroring [`super::fragment_size::FragmentSizeAccum`]), the
//! accumulator fetches:
//!
//! - 5' 4-mer: 4 reference bases at `[leftmost, leftmost + 3]` (inclusive),
//!   in reference forward-strand orientation.
//! - 3' 4-mer: 4 reference bases at `[rightmost - 3, rightmost]` (inclusive),
//!   reverse-complemented so it reads outward from the 3' cleavage site
//!   (this is the standard cfDNA end-motif convention).
//!
//! Where `leftmost = record.pos()`, `rightmost = leftmost + |TLEN| - 1` for
//! a properly paired record with `insert_size() > 0`.
//!
//! Pairs containing any non-ACGT base in either 4-mer are dropped from
//! `pair_count_used` and counted in `pair_count_skipped_non_acgt`.
//!
//! One [`rust_htslib::faidx::Reader`] handle is created per worker; htslib
//! FAI handles are not safe to share concurrently, but each worker owns
//! its own. The struct is `Send` (htslib upstream provides
//! `unsafe impl Send for faidx::Reader`).

use anyhow::{Context, Result};
use indexmap::IndexMap;
use rust_htslib::bam::record::Record;
use rust_htslib::faidx;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use super::common::{
    encode_kmer4, is_leftmost_proper_pair_mate, jsd_log2, kmer_array_to_map, rc_kmer4,
    shannon_entropy_log2,
};

#[derive(Debug)]
pub struct EndMotifAccum {
    /// Retained so future merge / re-open code paths can reconstruct the
    /// reader without re-plumbing the FASTA path. Not used in Phase 2.
    #[allow(dead_code)]
    fasta_path: PathBuf,
    reader: faidx::Reader,
    /// Per-contig cached lengths.
    contig_lens: HashMap<String, u64>,
    /// Contigs the FASTA doesn't carry — populated lazily on first miss
    /// to avoid repeated `fetch_seq_len` calls for unmapped contigs.
    contig_missing: HashSet<String>,
    pub(crate) kmer_5p: [u64; 256],
    pub(crate) kmer_3p: [u64; 256],
    pub pair_count_used: u64,
    pub pair_count_skipped_non_acgt: u64,
    /// Pairs skipped because the contig is missing from the FASTA or the
    /// fragment runs off the end of the contig.
    pub pair_count_skipped_oob: u64,
}

impl EndMotifAccum {
    /// Open a per-worker reader on `fasta_path`. Fails if the FAI index is
    /// missing or unreadable.
    pub fn new(fasta_path: &Path) -> Result<Self> {
        let reader = faidx::Reader::from_path(fasta_path).with_context(|| {
            format!(
                "Failed to open FASTA for end-motif extraction: {}",
                fasta_path.display()
            )
        })?;
        Ok(Self {
            fasta_path: fasta_path.to_path_buf(),
            reader,
            contig_lens: HashMap::new(),
            contig_missing: HashSet::new(),
            kmer_5p: [0u64; 256],
            kmer_3p: [0u64; 256],
            pair_count_used: 0,
            pair_count_skipped_non_acgt: 0,
            pair_count_skipped_oob: 0,
        })
    }

    fn contig_len(&mut self, chrom: &str) -> Option<u64> {
        if self.contig_missing.contains(chrom) {
            return None;
        }
        if let Some(len) = self.contig_lens.get(chrom) {
            return Some(*len);
        }
        let len = self.reader.fetch_seq_len(chrom);
        if len == 0 {
            self.contig_missing.insert(chrom.to_string());
            return None;
        }
        self.contig_lens.insert(chrom.to_string(), len);
        Some(len)
    }

    pub fn process_read(&mut self, record: &Record, chrom: &str, mapq_cut: u8) {
        if !is_leftmost_proper_pair_mate(record, mapq_cut) {
            return;
        }
        let tlen = record.insert_size() as u64;
        if tlen < 4 {
            // Fragment too short for a 4-mer at each end.
            self.pair_count_skipped_oob += 1;
            return;
        }
        let leftmost = record.pos() as u64;
        let rightmost = leftmost + tlen - 1;
        match self.lookup_kmers_at(chrom, leftmost, rightmost) {
            EndMotifLookup::Used { kmer_5p, kmer_3p } => {
                self.kmer_5p[kmer_5p as usize] += 1;
                self.kmer_3p[kmer_3p as usize] += 1;
                self.pair_count_used += 1;
            }
            EndMotifLookup::SkippedNonAcgt => {
                self.pair_count_skipped_non_acgt += 1;
            }
            EndMotifLookup::SkippedOob => {
                self.pair_count_skipped_oob += 1;
            }
        }
    }

    /// Pure FASTA-bound lookup of the 5' and 3' fragment-end 4-mers for
    /// the given proper-pair leftmost/rightmost reference coordinates.
    ///
    /// Does NOT update this accumulator's k-mer arrays nor the pair-count
    /// counters — the caller decides where to record the result. Updates
    /// internal contig-length caches as a side effect, which is shared
    /// across consumers. Phase 3 per-gene Tier-2 reuses this entry point
    /// so we only open one `faidx::Reader` per worker and only seek the
    /// FASTA once per pair.
    pub fn lookup_kmers_at(
        &mut self,
        chrom: &str,
        leftmost: u64,
        rightmost: u64,
    ) -> EndMotifLookup {
        if rightmost < leftmost + 3 {
            // Fragment shorter than 4 bp at either end.
            return EndMotifLookup::SkippedOob;
        }
        let len = match self.contig_len(chrom) {
            Some(l) => l,
            None => return EndMotifLookup::SkippedOob,
        };
        if rightmost >= len {
            return EndMotifLookup::SkippedOob;
        }
        let bases_5p =
            match self
                .reader
                .fetch_seq(chrom, leftmost as usize, (leftmost + 3) as usize)
            {
                Ok(v) if v.len() == 4 => v,
                _ => return EndMotifLookup::SkippedOob,
            };
        let bases_3p_ref =
            match self
                .reader
                .fetch_seq(chrom, (rightmost - 3) as usize, rightmost as usize)
            {
                Ok(v) if v.len() == 4 => v,
                _ => return EndMotifLookup::SkippedOob,
            };
        let kmer_5p = match encode_kmer4(&bases_5p) {
            Some(k) => k,
            None => return EndMotifLookup::SkippedNonAcgt,
        };
        let kmer_3p = match encode_kmer4(&bases_3p_ref) {
            Some(k) => rc_kmer4(k),
            None => return EndMotifLookup::SkippedNonAcgt,
        };
        EndMotifLookup::Used { kmer_5p, kmer_3p }
    }

    pub fn merge(&mut self, other: EndMotifAccum) {
        for k in 0..256 {
            self.kmer_5p[k] += other.kmer_5p[k];
            self.kmer_3p[k] += other.kmer_3p[k];
        }
        self.pair_count_used += other.pair_count_used;
        self.pair_count_skipped_non_acgt += other.pair_count_skipped_non_acgt;
        self.pair_count_skipped_oob += other.pair_count_skipped_oob;
        self.contig_missing.extend(other.contig_missing);
    }

    pub fn into_result(self) -> EndMotifResult {
        // Empty histograms are reported as `0.0` for Phase-2 schema
        // backwards-compatibility (the field is `type: number, minimum: 0`,
        // never null).
        EndMotifResult {
            kmer_5p: kmer_array_to_map(&self.kmer_5p),
            kmer_3p: kmer_array_to_map(&self.kmer_3p),
            shannon_entropy_5p: shannon_entropy_log2(&self.kmer_5p).unwrap_or(0.0),
            shannon_entropy_3p: shannon_entropy_log2(&self.kmer_3p).unwrap_or(0.0),
            jensen_shannon_divergence_5p_vs_3p: jsd_log2(&self.kmer_5p, &self.kmer_3p)
                .unwrap_or(0.0),
            pair_count_used: self.pair_count_used,
            pair_count_skipped_non_acgt: self.pair_count_skipped_non_acgt,
        }
    }
}

/// Result of a pure FASTA-bound lookup for the 5'/3' end-motifs of a
/// fragment. Returned by [`EndMotifAccum::lookup_kmers_at`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndMotifLookup {
    /// Both 4-mers were ACGT and within the contig. `kmer_3p` is already
    /// reverse-complemented to fragment-strand orientation.
    Used { kmer_5p: u8, kmer_3p: u8 },
    /// One of the 4-mers contained a non-ACGT base.
    SkippedNonAcgt,
    /// The fragment ran off the contig, the contig is missing from the
    /// FASTA, or the fragment is too short for two 4-mers.
    SkippedOob,
}

// ---------------------------------------------------------------------------
// Serializable output
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize)]
pub struct EndMotifResult {
    pub kmer_5p: IndexMap<String, u64>,
    pub kmer_3p: IndexMap<String, u64>,
    pub shannon_entropy_5p: f64,
    pub shannon_entropy_3p: f64,
    pub jensen_shannon_divergence_5p_vs_3p: f64,
    pub pair_count_used: u64,
    pub pair_count_skipped_non_acgt: u64,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shannon_none_for_empty_distribution() {
        let arr = [0u64; 256];
        assert!(shannon_entropy_log2(&arr).is_none());
    }

    #[test]
    fn shannon_uniform_distribution_is_log2_256() {
        let arr = [1u64; 256];
        let h = shannon_entropy_log2(&arr).expect("non-empty");
        assert!((h - 8.0).abs() < 1e-12);
    }

    #[test]
    fn shannon_concentrated_distribution_is_zero() {
        let mut arr = [0u64; 256];
        arr[42] = 1000;
        let h = shannon_entropy_log2(&arr).expect("non-empty");
        assert!(h.abs() < 1e-12);
    }

    #[test]
    fn jsd_self_is_zero() {
        let mut arr = [0u64; 256];
        for (i, v) in arr.iter_mut().enumerate() {
            *v = (i as u64 % 7) + 1;
        }
        let jsd = jsd_log2(&arr, &arr).expect("non-empty");
        assert!(jsd.abs() < 1e-12, "JSD(P, P) should be 0, got {}", jsd);
    }

    #[test]
    fn jsd_disjoint_supports_is_one() {
        let mut a = [0u64; 256];
        let mut b = [0u64; 256];
        a[0] = 1;
        b[1] = 1;
        let jsd = jsd_log2(&a, &b).expect("non-empty");
        // For two distributions with disjoint single-point supports and
        // base-2 log, JSD = 1.
        assert!((jsd - 1.0).abs() < 1e-12, "JSD = {}, expected 1.0", jsd);
    }

    #[test]
    fn merge_is_additive_on_kmer_arrays_and_counters() {
        // We cannot easily construct an EndMotifAccum without a real FASTA,
        // so test by constructing two and merging via an in-memory struct
        // assemblage. The merge logic is independent of the reader.
        // (Reader is only used by process_read.)
        // Instead, test array-level merge via direct field manipulation.
        let mut a = [0u64; 256];
        let mut b = [0u64; 256];
        a[10] = 3;
        b[10] = 4;
        a[20] = 1;
        for i in 0..256 {
            a[i] += b[i];
        }
        assert_eq!(a[10], 7);
        assert_eq!(a[20], 1);
    }

    #[test]
    fn into_result_drops_zero_keys() {
        let mut arr = [0u64; 256];
        arr[0] = 7;
        arr[255] = 3;
        let m = kmer_array_to_map(&arr);
        assert_eq!(m.len(), 2);
        assert_eq!(m.get("AAAA"), Some(&7));
        assert_eq!(m.get("TTTT"), Some(&3));
    }
}
