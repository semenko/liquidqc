//! Phase 2 Tier 1: FFT-based periodicity of the fragment-length histogram.
//!
//! Reads the 1 nt fragment-size histogram produced by
//! [`super::fragment_size::FragmentSizeAccum`], slices to the 50–500 nt
//! range (length 451), zero-pads to 512 (next power of 2 — better for
//! `realfft`'s mixed-radix planner), and computes the real FFT.
//!
//! The power spectrum is normalized to total power so each reported value
//! is a fraction in [0, 1]:
//!
//! - `power_helical_10_to_11_nt`: mean fractional power across periods 10
//!   and 11 nt (DNA helical pitch — a periodicity signature of nucleosome-
//!   protected fragmentation).
//! - `power_nucleosomal_145_to_170_nt`: mean fractional power across
//!   periods 145–170 nt (the dominant fragment-length mode for
//!   nucleosome-bound cfDNA / cfRNA cargo).
//! - `dominant_peak_period_nt`: period in nt where the spectrum peaks
//!   (excluding the DC bin). `None` if total_power == 0.
//! - `dominant_peak_power`: fractional power at that peak.
//! - `total_power`: sum of squared magnitudes (un-normalized) — multiply
//!   any of the above fractions by this to recover raw spectral weight.

use std::ops::Range;

use realfft::num_complex::Complex;
use realfft::RealFftPlanner;
use serde::Serialize;

/// Window over `tlen_hist` to FFT (50..=500 nt = length 451).
const HIST_RANGE: Range<usize> = 50..501;

/// Length to zero-pad to before the FFT. Next power of 2 above 451.
const PADDED_LEN: usize = 512;

/// Period bands of interest.
const HELICAL_PERIODS: Range<u32> = 10..12; // 10, 11
const NUCLEOSOMAL_PERIODS: Range<u32> = 145..171; // 145..=170

#[derive(Debug, Serialize)]
pub struct PeriodicityResult {
    pub power_helical_10_to_11_nt: f64,
    pub power_nucleosomal_145_to_170_nt: f64,
    pub dominant_peak_period_nt: Option<u32>,
    pub dominant_peak_power: f64,
    pub total_power: f64,
}

/// Compute the periodicity result from a 1 nt fragment-size histogram.
///
/// `tlen_hist[i]` is the count of pairs with `|TLEN| == i`. `tlen_hist`
/// must be at least `HIST_RANGE.end` long; shorter input is right-padded
/// with zeros. Returns an all-zero result when the windowed histogram is
/// uniformly zero.
pub fn compute(tlen_hist: &[u64]) -> PeriodicityResult {
    // Slice + zero-pad in one buffer. The non-zero region is signal[0..win_len].
    let win_len = HIST_RANGE.end - HIST_RANGE.start;
    let mut signal = vec![0.0f64; PADDED_LEN];
    for (out_i, hist_i) in HIST_RANGE.enumerate() {
        if hist_i < tlen_hist.len() {
            signal[out_i] = tlen_hist[hist_i] as f64;
        }
    }

    let signal_sum: f64 = signal[..win_len].iter().sum();
    if signal_sum == 0.0 {
        return PeriodicityResult {
            power_helical_10_to_11_nt: 0.0,
            power_nucleosomal_145_to_170_nt: 0.0,
            dominant_peak_period_nt: None,
            dominant_peak_power: 0.0,
            total_power: 0.0,
        };
    }

    // Subtract the in-window mean from the in-window samples only (leaving
    // zero-padding at zero). This removes the DC component without injecting
    // step discontinuities at the window/padding boundary, so a uniformly
    // flat windowed signal really does map to zero total power.
    let mean = signal_sum / (win_len as f64);
    for v in signal[..win_len].iter_mut() {
        *v -= mean;
    }

    let mut planner = RealFftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(PADDED_LEN);
    let mut spectrum = vec![Complex::<f64>::new(0.0, 0.0); PADDED_LEN / 2 + 1];
    fft.process(&mut signal, &mut spectrum)
        .expect("realfft size matches plan");

    // Power = |X|^2. Skip bin 0 (DC, already removed) when computing total.
    let powers: Vec<f64> = spectrum.iter().map(|c| c.norm_sqr()).collect();
    let total_power: f64 = powers.iter().skip(1).sum();
    if total_power == 0.0 {
        return PeriodicityResult {
            power_helical_10_to_11_nt: 0.0,
            power_nucleosomal_145_to_170_nt: 0.0,
            dominant_peak_period_nt: None,
            dominant_peak_power: 0.0,
            total_power: 0.0,
        };
    }

    // Period (in nt) for FFT bin k is PADDED_LEN / k. Iterate periods of
    // interest and accumulate fractional power.
    let band_power = |range: Range<u32>| -> f64 {
        let mut acc = 0.0;
        let mut n = 0u32;
        for period in range {
            // bin = PADDED_LEN / period (rounded). Skip if bin is 0 (DC) or
            // out of range.
            let bin = PADDED_LEN / period as usize;
            if bin == 0 || bin >= powers.len() {
                continue;
            }
            acc += powers[bin] / total_power;
            n += 1;
        }
        if n == 0 {
            0.0
        } else {
            acc / n as f64
        }
    };

    let power_helical = band_power(HELICAL_PERIODS);
    let power_nucleosomal = band_power(NUCLEOSOMAL_PERIODS);

    // Find dominant peak period. Skip bin 0 (DC, already removed) and bin 1
    // (period == PADDED_LEN, which is longer than the data window of length
    // `win_len`; not an observable periodicity, just residual low-frequency
    // content from any unimodal fragment-length hill). Cap from below so the
    // search range corresponds to periods <= `win_len`, the longest cycle
    // that can complete within the window.
    let min_bin = PADDED_LEN.div_ceil(win_len).max(2);
    let (peak_bin, peak_power) = powers
        .iter()
        .enumerate()
        .skip(min_bin)
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(k, p)| (k, *p))
        .unwrap_or((0, 0.0));
    let peak_period = if peak_bin == 0 {
        None
    } else {
        Some((PADDED_LEN / peak_bin) as u32)
    };
    let peak_frac = if total_power == 0.0 {
        0.0
    } else {
        peak_power / total_power
    };

    PeriodicityResult {
        power_helical_10_to_11_nt: power_helical,
        power_nucleosomal_145_to_170_nt: power_nucleosomal,
        dominant_peak_period_nt: peak_period,
        dominant_peak_power: peak_frac,
        total_power,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn empty_histogram_returns_zeros() {
        let r = compute(&[]);
        assert_eq!(r.total_power, 0.0);
        assert_eq!(r.dominant_peak_period_nt, None);
        assert_eq!(r.power_helical_10_to_11_nt, 0.0);
    }

    #[test]
    fn zero_histogram_returns_zeros() {
        let r = compute(&vec![0u64; 1000]);
        assert_eq!(r.total_power, 0.0);
        assert_eq!(r.dominant_peak_period_nt, None);
    }

    #[test]
    fn flat_histogram_has_no_dominant_periodicity() {
        // Constant signal: after mean subtraction it's all zeros, so total_power == 0.
        let r = compute(&vec![100u64; 1000]);
        assert_eq!(r.total_power, 0.0);
    }

    #[test]
    fn injected_period_167_recovers_near_167() {
        // Build a histogram with a strong sinusoid at period 167 nt.
        let period: f64 = 167.0;
        let mut hist = vec![0u64; 600];
        for (i, slot) in hist.iter_mut().enumerate().take(501).skip(50) {
            let phase = 2.0 * PI * (i as f64) / period;
            let amp = 1000.0 + 800.0 * phase.cos();
            *slot = amp.max(0.0) as u64;
        }
        let r = compute(&hist);
        assert!(r.total_power > 0.0);
        let peak = r.dominant_peak_period_nt.expect("non-empty signal");
        // Bin granularity at period 167 with PADDED_LEN=512 is rough — we
        // expect the bin nearest 167 to dominate, and PADDED_LEN/bin lands
        // close to the injected period. Assert within 25% to keep the test
        // robust to FFT-bin granularity at long periods.
        let err = (peak as f64 - period).abs() / period;
        assert!(
            err < 0.25,
            "peak {} too far from injected period {}",
            peak,
            period
        );
        // Nucleosomal band fractional power should dominate the helical
        // band for this input.
        assert!(r.power_nucleosomal_145_to_170_nt > r.power_helical_10_to_11_nt);
    }

    #[test]
    fn injected_period_10_recovers_near_10() {
        // Strong sinusoid at period 10 nt — DNA helical signature.
        let period: f64 = 10.0;
        let mut hist = vec![0u64; 600];
        for (i, slot) in hist.iter_mut().enumerate().take(501).skip(50) {
            let phase = 2.0 * PI * (i as f64) / period;
            let amp = 100.0 + 80.0 * phase.cos();
            *slot = amp.max(0.0) as u64;
        }
        let r = compute(&hist);
        assert!(r.total_power > 0.0);
        let peak = r.dominant_peak_period_nt.expect("non-empty signal");
        // Period 10 with PADDED_LEN=512 lands at bin 51 → period 10.0, very
        // close to integer.
        assert!((9..=11).contains(&peak), "peak {} not near 10", peak);
        // Helical band fraction should dominate the nucleosomal band here.
        assert!(r.power_helical_10_to_11_nt > r.power_nucleosomal_145_to_170_nt);
    }
}
