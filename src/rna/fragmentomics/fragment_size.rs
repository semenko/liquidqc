//! Phase 2 Tier 1: fragment-length distribution from `|TLEN|` of proper pairs.
//!
//! Single observation per pair: only the leftmost mate (first encountered in
//! coordinate-sorted order, i.e. `record.pos() <= record.mpos()`) contributes
//! `|record.insert_size()|` to the histogram. The same dedup mirror as
//! [`super::end_motifs::EndMotifAccum`] so size and end-motif counts agree on
//! `total_pairs_observed`.

use indexmap::IndexMap;
use rust_htslib::bam::record::Record;
use serde::Serialize;

use super::common::is_leftmost_proper_pair_mate;

/// Maximum fragment length tracked in 1 nt bins. Anything ≥ this falls into
/// the overflow counter and the `gt_500` named bin (the largest summary bin).
/// 1001 (i.e., 0..=1000) gives enough headroom for the 200-300 / 300-500 /
/// >500 summary bins while staying in single-page memory (8 KB per worker).
pub const HIST_LEN: usize = 1001;

#[derive(Debug, Default, Clone)]
pub struct FragmentSizeAccum {
    /// 1 nt bins for `|TLEN|` in [0, HIST_LEN).
    pub hist: Vec<u64>,
    /// Pairs observed with `|TLEN|` >= HIST_LEN.
    pub overflow: u64,
    /// Phase 4: leftmost proper-pair mates whose `|TLEN|` is shorter than
    /// `record.seq_len()` — i.e. the read literally extends past the
    /// opposite end of the fragment into adapter sequence. Counted on the
    /// same dispatch filter as the histogram so the denominator equals the
    /// histogram total + overflow. Earlier code used `|TLEN| < 2 * qlen`,
    /// which fires whenever the two reads merely overlap (the norm for
    /// short cfRNA fragments) and reported ~75% on clean libraries; the
    /// `|TLEN| < qlen` definition matches the actual adapter-contamination
    /// regime where each read's 3' end runs into adapter bases.
    pub adapter_readthrough_pairs: u64,
}

impl FragmentSizeAccum {
    pub fn new() -> Self {
        Self {
            hist: vec![0; HIST_LEN],
            overflow: 0,
            adapter_readthrough_pairs: 0,
        }
    }

    pub fn process_read(&mut self, record: &Record, mapq_cut: u8) {
        if !is_leftmost_proper_pair_mate(record, mapq_cut) {
            return;
        }
        let isize_abs = record.insert_size().unsigned_abs() as usize;
        if isize_abs < HIST_LEN {
            self.hist[isize_abs] += 1;
        } else {
            self.overflow += 1;
        }
        // Phase 4 adapter readthrough: pair where |TLEN| < leftmost_qlen,
        // i.e. the read length exceeds the fragment so each read's 3' end
        // extends past the opposite end of the fragment into adapter.
        let qlen = record.seq_len() as u64;
        if qlen > 0 && (isize_abs as u64) < qlen {
            self.adapter_readthrough_pairs += 1;
        }
    }

    pub fn merge(&mut self, other: FragmentSizeAccum) {
        if self.hist.len() != other.hist.len() {
            // Both should always be HIST_LEN; defensive in case the constant
            // ever changes mid-run.
            self.hist.resize(other.hist.len().max(self.hist.len()), 0);
        }
        for (i, v) in other.hist.iter().enumerate() {
            self.hist[i] += *v;
        }
        self.overflow += other.overflow;
        self.adapter_readthrough_pairs += other.adapter_readthrough_pairs;
    }

    pub fn into_result(self) -> FragmentSizeResult {
        let total: u64 = self.hist.iter().sum::<u64>() + self.overflow;
        let bins = compute_bins(&self.hist, self.overflow);

        let frac = |count: u64| -> f64 {
            if total == 0 {
                0.0
            } else {
                count as f64 / total as f64
            }
        };

        let count_lt_80: u64 = self.hist.iter().take(80).sum();
        // Half-open `[300, ∞)` to match the bin convention (`b300_500` = [300, 500),
        // `gt_500` = [500, ∞)); summing those two bins reproduces this count.
        let count_ge_300: u64 = self.hist.iter().skip(300).sum::<u64>() + self.overflow;

        let (mean, median) = mean_and_median(&self.hist, self.overflow, total);
        let adapter_readthrough_rate = frac(self.adapter_readthrough_pairs);

        FragmentSizeResult {
            bins,
            frac_lt_80: frac(count_lt_80),
            frac_ge_300: frac(count_ge_300),
            mean,
            median,
            total_pairs_observed: total,
            adapter_readthrough_rate,
            adapter_readthrough_pairs: self.adapter_readthrough_pairs,
            histogram: self.hist,
            overflow: self.overflow,
        }
    }
}

fn compute_bins(hist: &[u64], overflow: u64) -> FragmentSizeBins {
    // Half-open semantics matching the brief / plan:
    //   [0, 50), [50, 80), [80, 120), [120, 160), [160, 200),
    //   [200, 300), [300, 500), [500, ∞)
    let sum_range = |lo: usize, hi: usize| -> u64 {
        let hi = hi.min(hist.len());
        let lo = lo.min(hi);
        hist[lo..hi].iter().sum()
    };
    FragmentSizeBins {
        lt_50: sum_range(0, 50),
        b50_80: sum_range(50, 80),
        b80_120: sum_range(80, 120),
        b120_160: sum_range(120, 160),
        b160_200: sum_range(160, 200),
        b200_300: sum_range(200, 300),
        b300_500: sum_range(300, 500),
        gt_500: sum_range(500, hist.len()) + overflow,
    }
}

fn mean_and_median(hist: &[u64], overflow: u64, total: u64) -> (f64, u64) {
    if total == 0 {
        return (0.0, 0);
    }
    // Mean: overflow contributes its bin-floor (HIST_LEN as a lower bound) to
    // avoid misleading inflation; honest under-estimate is preferable to an
    // arbitrary midpoint when we don't know the upper tail.
    let in_range_sum: u128 = hist
        .iter()
        .enumerate()
        .map(|(i, c)| (i as u128) * (*c as u128))
        .sum();
    let total_sum = in_range_sum + (HIST_LEN as u128) * (overflow as u128);
    let mean = total_sum as f64 / total as f64;

    // Median: walk cumulative counts.
    let half = total.div_ceil(2);
    let mut cum: u64 = 0;
    let mut median: u64 = 0;
    for (i, c) in hist.iter().enumerate() {
        cum += *c;
        if cum >= half {
            median = i as u64;
            return (mean, median);
        }
    }
    // Median falls in overflow.
    if cum < half {
        median = HIST_LEN as u64;
    }
    (mean, median)
}

/// Serializable result.
#[derive(Debug, Serialize)]
pub struct FragmentSizeResult {
    pub bins: FragmentSizeBins,
    pub frac_lt_80: f64,
    pub frac_ge_300: f64,
    pub mean: f64,
    pub median: u64,
    pub total_pairs_observed: u64,
    /// Phase 4: pair-level adapter-readthrough rate. Defined as the fraction
    /// of leftmost proper-pair mates with `|TLEN| < record.seq_len()`.
    /// Surfaced into the envelope as a top-level scalar (see
    /// `adapter_readthrough_rate` in the schema).
    #[serde(skip)]
    pub adapter_readthrough_rate: f64,
    /// Raw numerator behind `adapter_readthrough_rate` (debug visibility).
    #[serde(skip)]
    #[allow(dead_code)]
    pub adapter_readthrough_pairs: u64,
    /// Raw histogram (1 nt bins, 0..HIST_LEN). Not serialized into the
    /// envelope JSON — used by the periodicity FFT downstream.
    #[serde(skip)]
    pub histogram: Vec<u64>,
    /// Pairs with |TLEN| >= HIST_LEN. Not serialized.
    /// Retained on the result so downstream consumers can detect tail
    /// truncation without re-reading the histogram.
    #[serde(skip)]
    #[allow(dead_code)]
    pub overflow: u64,
}

/// Eight-bin fragment-length histogram. Field names use leading-letter
/// prefixes (`b50_80` etc.) so they are valid Rust identifiers, then
/// rename via serde to the JSON keys (`50_80`, etc.) the schema expects.
#[derive(Debug, Serialize)]
pub struct FragmentSizeBins {
    pub lt_50: u64,
    #[serde(rename = "50_80")]
    pub b50_80: u64,
    #[serde(rename = "80_120")]
    pub b80_120: u64,
    #[serde(rename = "120_160")]
    pub b120_160: u64,
    #[serde(rename = "160_200")]
    pub b160_200: u64,
    #[serde(rename = "200_300")]
    pub b200_300: u64,
    #[serde(rename = "300_500")]
    pub b300_500: u64,
    pub gt_500: u64,
}

impl FragmentSizeBins {
    #[allow(dead_code)]
    pub fn total(&self) -> u64 {
        self.lt_50
            + self.b50_80
            + self.b80_120
            + self.b120_160
            + self.b160_200
            + self.b200_300
            + self.b300_500
            + self.gt_500
    }

    /// Convert to an IndexMap with the schema's exact JSON keys, preserving
    /// declaration order. Used by tests that index by string key.
    #[allow(dead_code)]
    pub fn as_indexmap(&self) -> IndexMap<&'static str, u64> {
        let mut m = IndexMap::with_capacity(8);
        m.insert("lt_50", self.lt_50);
        m.insert("50_80", self.b50_80);
        m.insert("80_120", self.b80_120);
        m.insert("120_160", self.b120_160);
        m.insert("160_200", self.b160_200);
        m.insert("200_300", self.b200_300);
        m.insert("300_500", self.b300_500);
        m.insert("gt_500", self.gt_500);
        m
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn hist_with(values: &[(usize, u64)]) -> Vec<u64> {
        let mut h = vec![0u64; HIST_LEN];
        for (i, c) in values {
            h[*i] = *c;
        }
        h
    }

    #[test]
    fn bin_boundaries_are_left_closed_right_open() {
        // One observation at each boundary: 49 → lt_50, 50 → 50_80,
        // 79 → 50_80, 80 → 80_120, 119 → 80_120, 120 → 120_160, etc.
        let hist = hist_with(&[
            (49, 1),
            (50, 1),
            (79, 1),
            (80, 1),
            (119, 1),
            (120, 1),
            (159, 1),
            (160, 1),
            (199, 1),
            (200, 1),
            (299, 1),
            (300, 1),
            (499, 1),
            (500, 1),
        ]);
        let bins = compute_bins(&hist, 0);
        assert_eq!(bins.lt_50, 1);
        assert_eq!(bins.b50_80, 2);
        assert_eq!(bins.b80_120, 2);
        assert_eq!(bins.b120_160, 2);
        assert_eq!(bins.b160_200, 2);
        assert_eq!(bins.b200_300, 2);
        assert_eq!(bins.b300_500, 2);
        assert_eq!(bins.gt_500, 1);
    }

    #[test]
    fn overflow_lands_in_gt_500() {
        let hist = vec![0u64; HIST_LEN];
        let bins = compute_bins(&hist, 7);
        assert_eq!(bins.gt_500, 7);
    }

    #[test]
    fn median_on_simple_histogram() {
        let mut h = vec![0u64; HIST_LEN];
        // Three observations at 100, 200, 300 — median is 200.
        h[100] = 1;
        h[200] = 1;
        h[300] = 1;
        let (_, median) = mean_and_median(&h, 0, 3);
        assert_eq!(median, 200);
    }

    #[test]
    fn median_falls_into_overflow_bin_when_warranted() {
        let mut h = vec![0u64; HIST_LEN];
        // Two in-range, three in overflow. Median index = 3 → overflow.
        h[10] = 1;
        h[20] = 1;
        let (_, median) = mean_and_median(&h, 3, 5);
        assert_eq!(median, HIST_LEN as u64);
    }

    #[test]
    fn fractions_use_total_including_overflow() {
        let mut accum = FragmentSizeAccum::new();
        accum.hist[40] = 1; // <80
        accum.hist[50] = 1; // <80
        accum.hist[400] = 1; // >300
        accum.overflow = 1; // >300
        let result = accum.into_result();
        assert_eq!(result.total_pairs_observed, 4);
        assert!((result.frac_lt_80 - 0.5).abs() < 1e-12);
        assert!((result.frac_ge_300 - 0.5).abs() < 1e-12);
    }

    #[test]
    fn empty_accumulator_returns_zeroed_result() {
        let accum = FragmentSizeAccum::new();
        let r = accum.into_result();
        assert_eq!(r.total_pairs_observed, 0);
        assert_eq!(r.median, 0);
        assert_eq!(r.mean, 0.0);
        assert_eq!(r.frac_lt_80, 0.0);
        assert_eq!(r.frac_ge_300, 0.0);
        assert_eq!(r.bins.total(), 0);
    }

    #[test]
    fn merge_is_additive() {
        let mut a = FragmentSizeAccum::new();
        a.hist[100] = 3;
        a.overflow = 1;
        a.adapter_readthrough_pairs = 2;
        let mut b = FragmentSizeAccum::new();
        b.hist[100] = 2;
        b.overflow = 4;
        b.adapter_readthrough_pairs = 5;
        a.merge(b);
        assert_eq!(a.hist[100], 5);
        assert_eq!(a.overflow, 5);
        assert_eq!(a.adapter_readthrough_pairs, 7);
    }

    #[test]
    fn adapter_readthrough_rate_is_zero_when_no_pairs() {
        let r = FragmentSizeAccum::new().into_result();
        assert_eq!(r.adapter_readthrough_pairs, 0);
        assert_eq!(r.adapter_readthrough_rate, 0.0);
    }

    #[test]
    fn adapter_readthrough_rate_is_numerator_over_total() {
        // 2 pairs in-range + 1 in overflow = 3 total. 1 of them flagged
        // as readthrough. Rate = 1/3.
        let mut accum = FragmentSizeAccum::new();
        accum.hist[100] = 2;
        accum.overflow = 1;
        accum.adapter_readthrough_pairs = 1;
        let r = accum.into_result();
        assert!((r.adapter_readthrough_rate - 1.0 / 3.0).abs() < 1e-12);
    }
}
