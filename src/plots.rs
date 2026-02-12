//! Plot generation for dupRust.
//!
//! Generates three types of plots matching dupRadar's output:
//! 1. Density scatter plot of duplication rate vs expression (primary plot)
//! 2. Boxplot of duplication rate by expression quantile bins
//! 3. Histogram of expression (RPK) distribution
//!
//! All plots are output as PNG files using the `plotters` crate.

use crate::dupmatrix::DupMatrix;
use crate::fitting::FitResult;
use anyhow::Result;
use plotters::prelude::*;

/// Color palette for density scatter plot (cyan -> blue -> green -> yellow -> red).
const DENSITY_COLORS: [(u8, u8, u8); 5] = [
    (0, 255, 255),   // cyan
    (0, 0, 255),     // blue
    (0, 255, 0),     // green
    (255, 255, 0),   // yellow
    (255, 0, 0),     // red
];

/// Interpolate between density colors based on a value 0.0 to 1.0.
fn density_color(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let n = DENSITY_COLORS.len() - 1;
    let idx = t * n as f64;
    let i = (idx.floor() as usize).min(n - 1);
    let frac = if i == n - 1 && idx >= n as f64 {
        1.0
    } else {
        idx - i as f64
    };

    let r = DENSITY_COLORS[i].0 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].0 as f64 * frac;
    let g = DENSITY_COLORS[i].1 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].1 as f64 * frac;
    let b = DENSITY_COLORS[i].2 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].2 as f64 * frac;

    RGBColor(r as u8, g as u8, b as u8)
}

/// Estimate 2D kernel density for scatter points using a simple grid-based approach.
///
/// Similar to R's `densCols()` with `nbin=500`.
fn estimate_density(x: &[f64], y: &[f64], nbins: usize) -> Vec<f64> {
    if x.is_empty() {
        return vec![];
    }

    let x_min = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let x_range = x_max - x_min;
    let y_range = y_max - y_min;

    if x_range == 0.0 || y_range == 0.0 {
        return vec![1.0; x.len()];
    }

    let x_step = x_range / nbins as f64;
    let y_step = y_range / nbins as f64;

    // Build 2D histogram
    let mut grid = vec![vec![0u32; nbins + 1]; nbins + 1];

    let to_bin = |val: f64, min: f64, step: f64| -> usize {
        ((val - min) / step).floor() as usize
    };

    for i in 0..x.len() {
        let bx = to_bin(x[i], x_min, x_step).min(nbins);
        let by = to_bin(y[i], y_min, y_step).min(nbins);
        grid[bx][by] += 1;
    }

    // Apply simple Gaussian smoothing (kernel radius = 2 bins)
    let radius = 2i32;
    let mut smoothed = vec![vec![0.0f64; nbins + 1]; nbins + 1];
    for bx in 0..=nbins {
        for by in 0..=nbins {
            if grid[bx][by] == 0 {
                continue;
            }
            let count = grid[bx][by] as f64;
            for dx in -radius..=radius {
                for dy in -radius..=radius {
                    let nx = bx as i32 + dx;
                    let ny = by as i32 + dy;
                    if nx >= 0 && nx <= nbins as i32 && ny >= 0 && ny <= nbins as i32 {
                        let weight =
                            (-((dx * dx + dy * dy) as f64) / (2.0 * 1.5 * 1.5)).exp();
                        smoothed[nx as usize][ny as usize] += count * weight;
                    }
                }
            }
        }
    }

    // Look up density for each point
    let mut densities = Vec::with_capacity(x.len());
    for i in 0..x.len() {
        let bx = to_bin(x[i], x_min, x_step).min(nbins);
        let by = to_bin(y[i], y_min, y_step).min(nbins);
        densities.push(smoothed[bx][by]);
    }

    // Normalize to [0, 1]
    let max_d = densities.iter().cloned().fold(0.0f64, f64::max);
    if max_d > 0.0 {
        for d in densities.iter_mut() {
            *d /= max_d;
        }
    }

    densities
}

/// Generate the density scatter plot of duplication rate vs expression.
///
/// This is the primary dupRadar plot (`duprateExpDensPlot`):
/// - X-axis: log10(RPK) with labels as 10^n
/// - Y-axis: duplication rate as percentage (0-100%)
/// - Points colored by local density (cyan -> red)
/// - Fitted logistic curve overlay
/// - Vertical threshold lines for 1 read/kbp and RPKM threshold
pub fn density_scatter_plot(
    dm: &DupMatrix,
    fit: &FitResult,
    rpkm_threshold: Option<f64>,
    output_path: &std::path::Path,
) -> Result<()> {
    let (width, height) = (800, 600);

    // Collect valid data points
    let mut x_data: Vec<f64> = Vec::new();
    let mut y_data: Vec<f64> = Vec::new();

    for row in &dm.rows {
        if row.rpk > 0.0 && row.dup_rate.is_finite() {
            x_data.push(row.rpk.log10());
            y_data.push(row.dup_rate * 100.0);
        }
    }

    if x_data.is_empty() {
        anyhow::bail!("No valid data points for density scatter plot");
    }

    // Compute density colors
    let densities = estimate_density(&x_data, &y_data, 500);

    // Determine axis ranges
    let x_min = x_data.iter().cloned().fold(f64::INFINITY, f64::min).floor();
    let x_max = x_data.iter().cloned().fold(f64::NEG_INFINITY, f64::max).ceil();
    let y_min = 0.0f64;
    let y_max = 100.0f64;

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("dupRadar - duplication rate", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc("Reads/kbp (RPK)")
        .y_desc("Duplication rate (%)")
        .x_label_formatter(&|v| {
            if *v == v.round() {
                format!("10^{}", *v as i32)
            } else {
                String::new()
            }
        })
        .draw()?;

    // Sort points by density (draw low density first, high density on top)
    let mut indexed: Vec<(usize, f64)> = densities.iter().enumerate().map(|(i, d)| (i, *d)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Draw scatter points
    for (i, _dens) in &indexed {
        let color = density_color(densities[*i]);
        chart.draw_series(std::iter::once(Circle::new(
            (x_data[*i], y_data[*i]),
            1,
            color.filled(),
        )))?;
    }

    // Draw fit curve
    let curve_points: Vec<(f64, f64)> = (0..200)
        .map(|i| {
            let x = x_min + (x_max - x_min) * i as f64 / 199.0;
            let y = fit.predict(x) * 100.0;
            (x, y)
        })
        .collect();

    chart.draw_series(LineSeries::new(
        curve_points,
        BLACK.stroke_width(2),
    ))?;

    // Draw "1 read/kbp" threshold (RPK = 1000 → log10 = 3)
    let rpk_1000 = 3.0f64;
    if rpk_1000 >= x_min && rpk_1000 <= x_max {
        chart.draw_series(LineSeries::new(
            vec![(rpk_1000, y_min), (rpk_1000, y_max)],
            RED.stroke_width(1),
        ))?;
    }

    // Draw RPKM threshold line
    if let Some(rpkm_thresh) = rpkm_threshold {
        if rpkm_thresh > 0.0 {
            let log_rpk = rpkm_thresh.log10();
            if log_rpk >= x_min && log_rpk <= x_max {
                chart.draw_series(LineSeries::new(
                    vec![(log_rpk, y_min), (log_rpk, y_max)],
                    GREEN.stroke_width(1),
                ))?;
            }
        }
    }

    // Add legend text
    let legend_text = format!(
        "int={:.4e}, slp={:.4e}",
        fit.intercept, fit.slope
    );
    // Draw legend as annotation
    chart.draw_series(std::iter::once(Text::new(
        legend_text,
        (x_min + 0.2, 95.0),
        ("sans-serif", 12).into_font(),
    )))?;

    root.present()?;
    Ok(())
}

/// Generate the duplication rate boxplot by expression quantile bins.
///
/// Matches dupRadar's `duprateExpBoxplot()`:
/// - Genes divided into 20 bins (5% quantile steps of RPK)
/// - Each box shows distribution of dupRate within that bin
pub fn duprate_boxplot(
    dm: &DupMatrix,
    output_path: &std::path::Path,
) -> Result<()> {
    let (width, height) = (900, 600);
    let step_size = 0.05;
    let n_bins = (1.0 / step_size) as usize;

    // Collect RPK and dupRate for genes with reads
    let mut gene_data: Vec<(f64, f64)> = dm
        .rows
        .iter()
        .filter(|r| r.all_counts > 0 && r.dup_rate.is_finite())
        .map(|r| (r.rpk, r.dup_rate))
        .collect();

    if gene_data.is_empty() {
        anyhow::bail!("No valid data for boxplot");
    }

    gene_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Compute RPK quantiles
    let rpk_values: Vec<f64> = gene_data.iter().map(|(r, _)| *r).collect();

    let mut bin_data: Vec<(String, Vec<f64>)> = Vec::new();

    for bin_idx in 0..n_bins {
        let p_low = bin_idx as f64 * step_size;
        let p_high = (bin_idx + 1) as f64 * step_size;

        let q_low = quantile(&rpk_values, p_low);
        let q_high = quantile(&rpk_values, p_high);

        let dup_rates: Vec<f64> = gene_data
            .iter()
            .filter(|(r, _)| {
                if bin_idx == 0 {
                    *r <= q_high
                } else {
                    *r > q_low && *r <= q_high
                }
            })
            .map(|(_, d)| *d)
            .collect();

        let mean_rpk: f64 = if dup_rates.is_empty() {
            0.0
        } else {
            gene_data
                .iter()
                .filter(|(r, _)| {
                    if bin_idx == 0 {
                        *r <= q_high
                    } else {
                        *r > q_low && *r <= q_high
                    }
                })
                .map(|(r, _)| *r)
                .sum::<f64>()
                / dup_rates.len() as f64
        };

        let label = format!(
            "{:.0}-{:.0}%\n{:.0}",
            p_low * 100.0,
            p_high * 100.0,
            mean_rpk
        );
        bin_data.push((label, dup_rates));
    }

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "dupRadar - duplication rate per expression bin",
            ("sans-serif", 18),
        )
        .margin(10)
        .x_label_area_size(80)
        .y_label_area_size(50)
        .build_cartesian_2d(0..n_bins, 0.0f64..1.0)?;

    chart
        .configure_mesh()
        .y_desc("Duplication rate")
        .x_desc("Expression quantile / mean RPK")
        .x_label_formatter(&|v| {
            let idx = *v;
            if idx < bin_data.len() {
                bin_data[idx].0.clone()
            } else {
                String::new()
            }
        })
        .x_labels(n_bins)
        .x_label_style(("sans-serif", 8).into_font().transform(FontTransform::Rotate270))
        .draw()?;

    // Draw boxplots
    for (idx, (_label, values)) in bin_data.iter().enumerate() {
        if values.is_empty() {
            continue;
        }

        let mut sorted = values.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let q1 = quantile(&sorted, 0.25);
        let median = quantile(&sorted, 0.5);
        let q3 = quantile(&sorted, 0.75);
        let iqr = q3 - q1;
        let whisker_low = sorted
            .iter()
            .find(|&&v| v >= q1 - 1.5 * iqr)
            .copied()
            .unwrap_or(q1);
        let whisker_high = sorted
            .iter()
            .rev()
            .find(|&&v| v <= q3 + 1.5 * iqr)
            .copied()
            .unwrap_or(q3);

        // Draw box

        // Convert to integer coordinates for the chart
        // We need to draw rectangles manually
        chart.draw_series(std::iter::once(Rectangle::new(
            [(idx, q1), (idx + 1, q3)],
            ShapeStyle {
                color: RGBAColor(70, 130, 180, 0.7),
                filled: true,
                stroke_width: 1,
            },
        )))?;

        // Box border
        chart.draw_series(std::iter::once(Rectangle::new(
            [(idx, q1), (idx + 1, q3)],
            BLACK.stroke_width(1),
        )))?;

        // Median line
        chart.draw_series(LineSeries::new(
            vec![(idx, median), (idx + 1, median)],
            BLACK.stroke_width(2),
        ))?;

        // Whiskers (vertical lines from center)
        // Lower whisker
        chart.draw_series(LineSeries::new(
            vec![(idx, whisker_low), (idx, q1)],
            BLACK.stroke_width(1),
        ))?;
        chart.draw_series(LineSeries::new(
            vec![(idx + 1, whisker_low), (idx + 1, q1)],
            BLACK.stroke_width(1),
        ))?;

        // Upper whisker  
        chart.draw_series(LineSeries::new(
            vec![(idx, q3), (idx, whisker_high)],
            BLACK.stroke_width(1),
        ))?;
        chart.draw_series(LineSeries::new(
            vec![(idx + 1, q3), (idx + 1, whisker_high)],
            BLACK.stroke_width(1),
        ))?;

        // Outliers
        for &v in &sorted {
            if v < whisker_low || v > whisker_high {
                chart.draw_series(std::iter::once(Circle::new(
                    (idx, v),
                    2,
                    BLACK.filled(),
                )))?;
            }
        }
    }

    root.present()?;
    Ok(())
}

/// Generate expression histogram.
///
/// Matches dupRadar's `expressionHist()`:
/// - Histogram of log10(RPK) with 100 bins
/// - Vertical red line at log10(1000) = 3 (1 read/kbp threshold)
pub fn expression_histogram(
    dm: &DupMatrix,
    output_path: &std::path::Path,
) -> Result<()> {
    let (width, height) = (800, 600);

    // Collect log10(RPK) for genes with reads
    let log_rpk: Vec<f64> = dm
        .rows
        .iter()
        .filter(|r| r.rpk > 0.0)
        .map(|r| r.rpk.log10())
        .collect();

    if log_rpk.is_empty() {
        anyhow::bail!("No valid data for expression histogram");
    }

    let x_min = log_rpk.iter().cloned().fold(f64::INFINITY, f64::min).floor();
    let x_max = log_rpk.iter().cloned().fold(f64::NEG_INFINITY, f64::max).ceil();

    // Build histogram
    let n_bins = 100;
    let bin_width = (x_max - x_min) / n_bins as f64;
    let mut hist = vec![0u32; n_bins];

    for &v in &log_rpk {
        let bin = ((v - x_min) / bin_width).floor() as usize;
        let bin = bin.min(n_bins - 1);
        hist[bin] += 1;
    }

    let y_max = *hist.iter().max().unwrap_or(&1) as f64;

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("dupRadar - expression histogram", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, 0.0..y_max * 1.05)?;

    chart
        .configure_mesh()
        .x_desc("Reads/kbp (RPK)")
        .y_desc("Frequency")
        .x_label_formatter(&|v| {
            if *v == v.round() {
                format!("10^{}", *v as i32)
            } else {
                String::new()
            }
        })
        .draw()?;

    // Draw histogram bars
    chart.draw_series(
        hist.iter().enumerate().map(|(i, &count)| {
            let x0 = x_min + i as f64 * bin_width;
            let x1 = x0 + bin_width;
            Rectangle::new(
                [(x0, 0.0), (x1, count as f64)],
                ShapeStyle {
                    color: RGBAColor(70, 130, 180, 0.8),
                    filled: true,
                    stroke_width: 1,
                },
            )
        }),
    )?;

    // Draw "1 read/kbp" threshold line at log10(1000) = 3
    let threshold = 3.0f64;
    if threshold >= x_min && threshold <= x_max {
        chart.draw_series(LineSeries::new(
            vec![(threshold, 0.0), (threshold, y_max * 1.05)],
            RED.stroke_width(2),
        ))?;
    }

    root.present()?;
    Ok(())
}

/// Write the intercept and slope values to a file (matching R's format: label\tvalue per line).
pub fn write_intercept_slope(fit: &FitResult, path: &std::path::Path) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "intercept\t{}", fit.intercept)?;
    writeln!(f, "slope\t{}", fit.slope)?;
    Ok(())
}

/// Write MultiQC-compatible intercept file.
pub fn write_mqc_intercept(
    fit: &FitResult,
    sample_name: &str,
    path: &std::path::Path,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "# id: dupRadar")?;
    writeln!(f, "# plot_type: 'generalstats'")?;
    writeln!(f, "# pconfig:")?;
    writeln!(f, "#     dupRadar_intercept:")?;
    writeln!(f, "#         title: 'dupRadar int'")?;
    writeln!(f, "#         namespace: 'dupRadar'")?;
    writeln!(f, "#         description: 'dupRadar duplication rate at low read counts'")?;
    writeln!(f, "#         max: 100")?;
    writeln!(f, "#         min: 0")?;
    writeln!(f, "#         format: '{{:.2f}}'")?;
    writeln!(f, "Sample\tdupRadar_intercept")?;
    writeln!(f, "{}\t{}", sample_name, fit.intercept)?;
    Ok(())
}

/// Write MultiQC-compatible duplication rate curve data.
pub fn write_mqc_curve(
    fit: &FitResult,
    dm: &DupMatrix,
    path: &std::path::Path,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;

    writeln!(f, "# id: 'dupradar'")?;
    writeln!(f, "# section_name: 'dupRadar'")?;
    writeln!(f, "# description: 'Duplication rate vs expression'")?;
    writeln!(f, "# plot_type: 'linegraph'")?;
    writeln!(f, "# pconfig:")?;
    writeln!(f, "#     title: 'dupRadar General Linear Model'")?;
    writeln!(f, "#     xlab: 'Expression (reads/kbp)'")?;
    writeln!(f, "#     ylab: 'Duplication rate (%)'")?;
    writeln!(f, "#     ymin: 0")?;
    writeln!(f, "#     ymax: 100")?;
    writeln!(f, "#     xlog: True")?;

    // Generate smooth curve from the logistic fit across the RPK range
    // Use 100 evenly-spaced points in log10 space for a smooth curve
    let rpk_values: Vec<f64> = dm
        .rows
        .iter()
        .filter(|r| r.rpk > 0.0)
        .map(|r| r.rpk)
        .collect();

    if rpk_values.is_empty() {
        return Ok(());
    }

    let min_rpk = rpk_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max_rpk = rpk_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    let log_min = min_rpk.log10();
    let log_max = max_rpk.log10();
    let n_points = 100;

    // Write in two-column format: RPK (reads/kbp) and dup rate (%)
    writeln!(f, "RPK\tDuplication Rate (%)")?;

    for i in 0..=n_points {
        let log_rpk = log_min + (log_max - log_min) * (i as f64) / (n_points as f64);
        let rpk = 10.0_f64.powf(log_rpk);
        let dup_pct = fit.predict_rpk(rpk) * 100.0;
        writeln!(f, "{}\t{}", rpk, dup_pct)?;
    }

    Ok(())
}

/// Compute a quantile from a sorted array.
fn quantile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if p <= 0.0 {
        return sorted[0];
    }
    if p >= 1.0 {
        return *sorted.last().unwrap();
    }

    let n = sorted.len();
    let idx = p * (n - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx.fract();

    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_color_endpoints() {
        let c0 = density_color(0.0);
        assert_eq!((c0.0, c0.1, c0.2), (0, 255, 255)); // cyan

        let c1 = density_color(1.0);
        assert_eq!((c1.0, c1.1, c1.2), (255, 0, 0)); // red
    }

    #[test]
    fn test_quantile() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((quantile(&data, 0.0) - 1.0).abs() < 1e-10);
        assert!((quantile(&data, 0.5) - 3.0).abs() < 1e-10);
        assert!((quantile(&data, 1.0) - 5.0).abs() < 1e-10);
        assert!((quantile(&data, 0.25) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_density() {
        // Simple test: clustered points should have higher density
        let x = vec![0.0, 0.1, 0.0, 0.1, 5.0];
        let y = vec![0.0, 0.1, 0.1, 0.0, 5.0];
        let d = estimate_density(&x, &y, 100);
        assert_eq!(d.len(), 5);
        // The cluster of 4 points should have higher density than the outlier
        assert!(d[0] > d[4]);
    }
}
