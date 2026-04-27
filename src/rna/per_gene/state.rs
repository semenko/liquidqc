//! Per-gene accumulator state.
//!
//! `PerGeneState` is allocated as `Vec<PerGeneState>` indexed by the global
//! `GeneIdx` from `dupradar::counting`. Each per-worker accumulator owns one
//! such vector with zeroed state for every gene in the GTF; merging across
//! workers is field-wise addition.
//!
//! The 4-mer histograms (`em_5p` / `em_3p`) are lazy: the 1024 B per array
//! is only allocated on the first end-motif observation. Genes that never
//! receive a paired+FASTA fragment keep `None` and contribute zero to peak
//! memory. On a 60k-gene GENCODE GTF with 4 workers, this drops peak per-
//! worker state from ~530 MB to ~30 MB until end-motif lookups start landing.

/// Number of distinct 4-mers (4^4).
pub const KMER4_CARD: usize = 256;

/// All-zero 4-mer histogram. Returned as the default view when a gene has
/// no end-motif observations yet — keeps the writer / entropy paths
/// agnostic to whether the lazy `Box<[u32; 256]>` has been allocated.
const ZERO_KMER_HIST: [u32; KMER4_CARD] = [0; KMER4_CARD];

/// Field-wise additive merge of two lazy 4-mer histograms. Cheap when
/// both are `None`; allocates the receiver only when the donor has data.
fn merge_kmer_lazy(dst: &mut Option<Box<[u32; KMER4_CARD]>>, src: Option<Box<[u32; KMER4_CARD]>>) {
    let Some(src) = src else { return };
    match dst {
        Some(d) => {
            for i in 0..KMER4_CARD {
                d[i] = d[i].saturating_add(src[i]);
            }
        }
        None => *dst = Some(src),
    }
}

/// Width of the eight-bin TLEN summary histogram (matches the sample-level
/// fragmentomics bin scheme: lt_50, 50_80, 80_120, 120_160, 160_200,
/// 200_300, 300_500, gt_500).
pub const TLEN_BIN_COUNT: usize = 8;

/// Number of decile slots for 5'→3' gene-body coverage.
pub const DECILE_COUNT: usize = 10;

/// Maximum `|TLEN|` tracked for the per-gene mean denominator before it
/// counts as overflow. Matches `crate::rna::fragmentomics::fragment_size::HIST_LEN`.
pub const TLEN_OVERFLOW_THRESHOLD: u32 = 1001;

/// Per-gene accumulator state. One instance per gene per worker.
///
/// All fields default to zero. `Default::default()` produces a fresh,
/// empty state. `Default` is implemented manually because `[u32; 256]`
/// exceeds the array length the derive macro auto-implements (32).
#[derive(Debug, Clone)]
pub struct PerGeneState {
    // ---------------------------------------------------------------
    // Fragment / read counts (independent of paired-end-ness)
    // ---------------------------------------------------------------
    /// Number of fragments unambiguously assigned to this gene by the
    /// counting dispatcher (`combined_genes.len() == 1` paths).
    pub fragment_count_assigned: u32,
    /// Primary mapped reads (single-end) or primary mapped mate1 reads
    /// (paired-end leftmost-mate convention) whose unambiguous gene
    /// assignment is this gene.
    pub primary_reads: u32,
    /// Primary reads with a soft clip at the 5' end of the fragment
    /// (strand-corrected by the dispatcher).
    pub primary_reads_with_5p_clip: u32,
    /// Primary reads with a soft clip at the 3' end of the fragment.
    pub primary_reads_with_3p_clip: u32,

    // ---------------------------------------------------------------
    // Fragment length (paired-end leftmost-mate proper pairs only)
    // ---------------------------------------------------------------
    /// 8-bin coarse histogram of `|TLEN|`, identical bin layout to the
    /// sample-level `FragmentSizeBins`.
    pub tlen_bins: [u32; TLEN_BIN_COUNT],
    /// Pairs with `|TLEN| >= TLEN_OVERFLOW_THRESHOLD` (1001).
    pub tlen_overflow: u32,
    /// Sum of in-range `|TLEN|` for mean computation.
    pub tlen_sum: u64,
    /// Count of in-range observations (excludes overflow). Used for mean.
    pub tlen_count: u32,
    /// Min/max in-range `|TLEN|`. Sentinel `u32::MAX` means unset for min;
    /// 0 for max means unset.
    pub tlen_min: u32,
    pub tlen_max: u32,

    // ---------------------------------------------------------------
    // End motifs (paired-end + FASTA provided only)
    // ---------------------------------------------------------------
    /// 5' fragment-end 4-mer counts (encoded ACGT base-4). Lazy-allocated
    /// — `None` until the first observation, which drops the bulk of
    /// per-gene memory for genes that never see a paired+FASTA fragment.
    pub em_5p: Option<Box<[u32; KMER4_CARD]>>,
    /// 3' fragment-end 4-mer counts (reverse-complemented to fragment
    /// strand). Same lazy-allocation contract as [`Self::em_5p`].
    pub em_3p: Option<Box<[u32; KMER4_CARD]>>,
    /// Pairs that contributed a 4-mer to both ends.
    pub em_pair_count_used: u32,
    /// Pairs dropped because either 4-mer contained a non-ACGT base.
    pub em_pair_count_skipped_non_acgt: u32,
    /// Pairs dropped because the fragment ran off the FASTA contig or the
    /// contig was missing.
    pub em_pair_count_skipped_oob: u32,

    // ---------------------------------------------------------------
    // 5'→3' gene-body decile coverage
    // ---------------------------------------------------------------
    /// One increment per primary aligned-block midpoint that maps to the
    /// gene's merged-exon mRNA coordinate space, binned into 10 deciles.
    /// `decile_coverage[0]` is the 5'-most decile (strand-corrected).
    pub decile_coverage: [u32; DECILE_COUNT],

    // ---------------------------------------------------------------
    // Intron-exon accounting
    // ---------------------------------------------------------------
    /// Aligned bp from primary reads that overlap this gene's merged exons.
    pub exon_aligned_bp: u64,
    /// Aligned bp from primary reads that overlap this gene's intronic
    /// regions (gene span minus merged-exon union).
    pub intron_aligned_bp: u64,
}

impl Default for PerGeneState {
    fn default() -> Self {
        Self {
            fragment_count_assigned: 0,
            primary_reads: 0,
            primary_reads_with_5p_clip: 0,
            primary_reads_with_3p_clip: 0,
            tlen_bins: [0; TLEN_BIN_COUNT],
            tlen_overflow: 0,
            tlen_sum: 0,
            tlen_count: 0,
            tlen_min: 0,
            tlen_max: 0,
            em_5p: None,
            em_3p: None,
            em_pair_count_used: 0,
            em_pair_count_skipped_non_acgt: 0,
            em_pair_count_skipped_oob: 0,
            decile_coverage: [0; DECILE_COUNT],
            exon_aligned_bp: 0,
            intron_aligned_bp: 0,
        }
    }
}

impl PerGeneState {
    /// Field-wise additive merge.
    pub fn merge(&mut self, other: PerGeneState) {
        self.fragment_count_assigned += other.fragment_count_assigned;
        self.primary_reads += other.primary_reads;
        self.primary_reads_with_5p_clip += other.primary_reads_with_5p_clip;
        self.primary_reads_with_3p_clip += other.primary_reads_with_3p_clip;

        for i in 0..TLEN_BIN_COUNT {
            self.tlen_bins[i] = self.tlen_bins[i].saturating_add(other.tlen_bins[i]);
        }
        self.tlen_overflow = self.tlen_overflow.saturating_add(other.tlen_overflow);
        self.tlen_sum = self.tlen_sum.saturating_add(other.tlen_sum);
        self.tlen_count = self.tlen_count.saturating_add(other.tlen_count);
        // tlen_min: take the smaller non-zero, treating 0 as unset on both sides.
        self.tlen_min = match (self.tlen_min, other.tlen_min) {
            (0, x) | (x, 0) => x,
            (a, b) => a.min(b),
        };
        self.tlen_max = self.tlen_max.max(other.tlen_max);

        merge_kmer_lazy(&mut self.em_5p, other.em_5p);
        merge_kmer_lazy(&mut self.em_3p, other.em_3p);
        self.em_pair_count_used += other.em_pair_count_used;
        self.em_pair_count_skipped_non_acgt += other.em_pair_count_skipped_non_acgt;
        self.em_pair_count_skipped_oob += other.em_pair_count_skipped_oob;

        for i in 0..DECILE_COUNT {
            self.decile_coverage[i] =
                self.decile_coverage[i].saturating_add(other.decile_coverage[i]);
        }
        self.exon_aligned_bp = self.exon_aligned_bp.saturating_add(other.exon_aligned_bp);
        self.intron_aligned_bp = self
            .intron_aligned_bp
            .saturating_add(other.intron_aligned_bp);
    }

    /// Borrow `em_5p` as a flat `&[u32; 256]`; returns the static zero
    /// histogram when no end-motif has been observed for this gene yet.
    pub fn em_5p_view(&self) -> &[u32; KMER4_CARD] {
        self.em_5p.as_deref().unwrap_or(&ZERO_KMER_HIST)
    }

    /// Borrow `em_3p` as a flat `&[u32; 256]`; same lazy-zero contract.
    pub fn em_3p_view(&self) -> &[u32; KMER4_CARD] {
        self.em_3p.as_deref().unwrap_or(&ZERO_KMER_HIST)
    }

    /// Mutable access to `em_5p`, allocating on first call.
    pub fn em_5p_or_alloc(&mut self) -> &mut [u32; KMER4_CARD] {
        self.em_5p
            .get_or_insert_with(|| Box::new([0; KMER4_CARD]))
            .as_mut()
    }

    /// Mutable access to `em_3p`, allocating on first call.
    pub fn em_3p_or_alloc(&mut self) -> &mut [u32; KMER4_CARD] {
        self.em_3p
            .get_or_insert_with(|| Box::new([0; KMER4_CARD]))
            .as_mut()
    }

    /// Map a 1-bp `|TLEN|` value into the 8-bin coarse histogram index.
    /// Returns `None` if the value is in overflow (≥ `TLEN_OVERFLOW_THRESHOLD`).
    pub fn tlen_bin_index(tlen_abs: u32) -> Option<usize> {
        // Half-open bins: [0,50), [50,80), [80,120), [120,160),
        // [160,200), [200,300), [300,500), [500, OVERFLOW).
        if tlen_abs >= TLEN_OVERFLOW_THRESHOLD {
            return None;
        }
        Some(if tlen_abs < 50 {
            0
        } else if tlen_abs < 80 {
            1
        } else if tlen_abs < 120 {
            2
        } else if tlen_abs < 160 {
            3
        } else if tlen_abs < 200 {
            4
        } else if tlen_abs < 300 {
            5
        } else if tlen_abs < 500 {
            6
        } else {
            7
        })
    }

    /// Record one `|TLEN|` observation.
    pub fn record_tlen(&mut self, tlen_abs: u32) {
        if let Some(bin) = Self::tlen_bin_index(tlen_abs) {
            self.tlen_bins[bin] = self.tlen_bins[bin].saturating_add(1);
            self.tlen_sum = self.tlen_sum.saturating_add(tlen_abs as u64);
            self.tlen_count = self.tlen_count.saturating_add(1);
            if self.tlen_min == 0 || tlen_abs < self.tlen_min {
                self.tlen_min = tlen_abs;
            }
            if tlen_abs > self.tlen_max {
                self.tlen_max = tlen_abs;
            }
        } else {
            self.tlen_overflow = self.tlen_overflow.saturating_add(1);
        }
    }

    /// Mean `|TLEN|` over in-range observations. `None` when `tlen_count == 0`.
    pub fn tlen_mean(&self) -> Option<f64> {
        if self.tlen_count == 0 {
            None
        } else {
            Some(self.tlen_sum as f64 / self.tlen_count as f64)
        }
    }

    /// Approximate median TLEN using bin midpoints and overflow lower bound.
    /// Returns `None` when no observations exist.
    pub fn tlen_median_approx(&self) -> Option<f64> {
        let total = self.tlen_count as u64 + self.tlen_overflow as u64;
        if total == 0 {
            return None;
        }
        // Cumulative walk over bins; treat overflow as a single bin at
        // the lower bound (TLEN_OVERFLOW_THRESHOLD).
        let half = total.div_ceil(2);
        let bin_midpoints: [f64; TLEN_BIN_COUNT] = [
            25.0,  // [0, 50)
            65.0,  // [50, 80)
            100.0, // [80, 120)
            140.0, // [120, 160)
            180.0, // [160, 200)
            250.0, // [200, 300)
            400.0, // [300, 500)
            750.0, // [500, OVERFLOW)
        ];
        let mut cum: u64 = 0;
        for (i, &count) in self.tlen_bins.iter().enumerate() {
            cum += count as u64;
            if cum >= half {
                return Some(bin_midpoints[i]);
            }
        }
        // Falls in overflow.
        Some(TLEN_OVERFLOW_THRESHOLD as f64)
    }

    /// Soft-clip rate at the 5' end relative to `primary_reads`.
    /// `None` when no primary reads were observed.
    pub fn soft_clip_rate_5p(&self) -> Option<f64> {
        if self.primary_reads == 0 {
            None
        } else {
            Some(self.primary_reads_with_5p_clip as f64 / self.primary_reads as f64)
        }
    }

    /// Soft-clip rate at the 3' end relative to `primary_reads`.
    pub fn soft_clip_rate_3p(&self) -> Option<f64> {
        if self.primary_reads == 0 {
            None
        } else {
            Some(self.primary_reads_with_3p_clip as f64 / self.primary_reads as f64)
        }
    }

    /// Intron-exon ratio (intron / (intron + exon)). `None` when neither
    /// region was hit by any primary read.
    pub fn intron_exon_ratio(&self) -> Option<f64> {
        let total = self.exon_aligned_bp + self.intron_aligned_bp;
        if total == 0 {
            None
        } else {
            Some(self.intron_aligned_bp as f64 / total as f64)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tlen_bin_boundaries() {
        assert_eq!(PerGeneState::tlen_bin_index(0), Some(0));
        assert_eq!(PerGeneState::tlen_bin_index(49), Some(0));
        assert_eq!(PerGeneState::tlen_bin_index(50), Some(1));
        assert_eq!(PerGeneState::tlen_bin_index(79), Some(1));
        assert_eq!(PerGeneState::tlen_bin_index(80), Some(2));
        assert_eq!(PerGeneState::tlen_bin_index(119), Some(2));
        assert_eq!(PerGeneState::tlen_bin_index(120), Some(3));
        assert_eq!(PerGeneState::tlen_bin_index(159), Some(3));
        assert_eq!(PerGeneState::tlen_bin_index(160), Some(4));
        assert_eq!(PerGeneState::tlen_bin_index(199), Some(4));
        assert_eq!(PerGeneState::tlen_bin_index(200), Some(5));
        assert_eq!(PerGeneState::tlen_bin_index(299), Some(5));
        assert_eq!(PerGeneState::tlen_bin_index(300), Some(6));
        assert_eq!(PerGeneState::tlen_bin_index(499), Some(6));
        assert_eq!(PerGeneState::tlen_bin_index(500), Some(7));
        assert_eq!(PerGeneState::tlen_bin_index(1000), Some(7));
        assert_eq!(PerGeneState::tlen_bin_index(1001), None);
    }

    #[test]
    fn record_tlen_updates_bins_sum_min_max() {
        let mut s = PerGeneState::default();
        s.record_tlen(75); // bin 1
        s.record_tlen(150); // bin 3
        s.record_tlen(2000); // overflow
        assert_eq!(s.tlen_bins[1], 1);
        assert_eq!(s.tlen_bins[3], 1);
        assert_eq!(s.tlen_overflow, 1);
        assert_eq!(s.tlen_count, 2);
        assert_eq!(s.tlen_sum, 75 + 150);
        assert_eq!(s.tlen_min, 75);
        assert_eq!(s.tlen_max, 150);
        assert!((s.tlen_mean().unwrap() - 112.5).abs() < 1e-12);
    }

    #[test]
    fn tlen_mean_none_when_empty() {
        let s = PerGeneState::default();
        assert!(s.tlen_mean().is_none());
        assert!(s.tlen_median_approx().is_none());
    }

    #[test]
    fn tlen_median_falls_in_overflow_when_warranted() {
        let mut s = PerGeneState::default();
        s.record_tlen(75); // bin 1 (1 in-range)
        s.record_tlen(2000); // overflow
        s.record_tlen(2000); // overflow
                             // total = 3, half = 2; cum after bin 1 = 1 < 2 → median in overflow.
        let med = s.tlen_median_approx().unwrap();
        assert!((med - TLEN_OVERFLOW_THRESHOLD as f64).abs() < 1e-12);
    }

    #[test]
    fn merge_is_additive() {
        let mut a = PerGeneState::default();
        a.tlen_bins[2] = 3;
        a.em_5p_or_alloc()[10] = 5;
        a.exon_aligned_bp = 100;
        a.tlen_min = 80;
        a.tlen_max = 120;

        let mut b = PerGeneState::default();
        b.tlen_bins[2] = 2;
        b.em_5p_or_alloc()[10] = 1;
        b.em_5p_or_alloc()[20] = 7;
        b.exon_aligned_bp = 50;
        b.intron_aligned_bp = 25;
        b.tlen_min = 50;
        b.tlen_max = 200;

        a.merge(b);
        assert_eq!(a.tlen_bins[2], 5);
        assert_eq!(a.em_5p_view()[10], 6);
        assert_eq!(a.em_5p_view()[20], 7);
        assert_eq!(a.exon_aligned_bp, 150);
        assert_eq!(a.intron_aligned_bp, 25);
        assert_eq!(a.tlen_min, 50);
        assert_eq!(a.tlen_max, 200);
    }

    #[test]
    fn merge_em_lazy_takes_donor_when_receiver_unset() {
        let a_before_total = {
            let mut a = PerGeneState::default();
            assert!(a.em_5p.is_none());
            let mut b = PerGeneState::default();
            b.em_5p_or_alloc()[42] = 9;
            a.merge(b);
            // Receiver picked up donor's box — no double allocation.
            assert!(a.em_5p.is_some());
            assert_eq!(a.em_5p_view()[42], 9);
            a.em_5p_view().iter().sum::<u32>()
        };
        assert_eq!(a_before_total, 9);
    }

    #[test]
    fn merge_em_lazy_no_alloc_when_both_unset() {
        let mut a = PerGeneState::default();
        let b = PerGeneState::default();
        a.merge(b);
        assert!(a.em_5p.is_none());
        assert!(a.em_3p.is_none());
    }

    #[test]
    fn merge_min_treats_zero_as_unset() {
        let mut a = PerGeneState::default(); // tlen_min = 0 (unset)
        let mut b = PerGeneState::default();
        b.tlen_min = 75;
        a.merge(b);
        assert_eq!(a.tlen_min, 75);
    }

    #[test]
    fn rates_are_none_without_primary_reads() {
        let s = PerGeneState::default();
        assert!(s.soft_clip_rate_5p().is_none());
        assert!(s.soft_clip_rate_3p().is_none());
        assert!(s.intron_exon_ratio().is_none());
    }

    #[test]
    fn rates_use_primary_reads_denominator() {
        let mut s = PerGeneState::default();
        s.primary_reads = 100;
        s.primary_reads_with_5p_clip = 25;
        s.primary_reads_with_3p_clip = 5;
        s.exon_aligned_bp = 800;
        s.intron_aligned_bp = 200;
        assert!((s.soft_clip_rate_5p().unwrap() - 0.25).abs() < 1e-12);
        assert!((s.soft_clip_rate_3p().unwrap() - 0.05).abs() < 1e-12);
        assert!((s.intron_exon_ratio().unwrap() - 0.2).abs() < 1e-12);
    }
}
