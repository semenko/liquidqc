//! 5'→3' decile binning along the merged-exon mRNA coordinate space.
//!
//! Per primary read assigned to a gene, each of its aligned blocks
//! contributes one increment to the decile that contains the block's
//! exonic-midpoint (in plus-strand mRNA coordinates, then strand-flipped
//! for `'-'` genes). Blocks lying entirely in introns/intergenic space
//! are skipped (they never overlap the gene's exons). Blocks that span
//! an exon boundary or partial-exon are mapped using the midpoint of
//! the exonic intersection with the block.

use super::shape::GeneShape;
use super::state::DECILE_COUNT;

/// Add decile-coverage increments for a read's aligned blocks against a
/// gene's merged-exon geometry.
///
/// `aligned_blocks` are 0-based half-open `[start, end)` reference
/// coords (as produced by `cigar_to_aligned_blocks` in counting.rs).
/// The exon coords on `shape` are 1-based inclusive (GTF convention),
/// so we convert when intersecting.
///
/// One block may straddle multiple exons of a spliced read; each
/// non-empty exonic-overlap segment contributes exactly one increment
/// using its segment midpoint. This makes the per-block contribution
/// stable under spliced-vs-unspliced reads of the same length.
pub fn add_block_deciles(
    shape: &GeneShape,
    aligned_blocks: &[(u64, u64)],
    deciles: &mut [u32; DECILE_COUNT],
) {
    if shape.total_exonic == 0 {
        return;
    }

    for &(b_start_0h, b_end_0h) in aligned_blocks {
        // Convert 0-based half-open block to 1-based inclusive.
        if b_end_0h <= b_start_0h {
            continue;
        }
        let block_lo = b_start_0h + 1;
        let block_hi = b_end_0h; // 1-based inclusive

        // Intersect with each merged exon.
        for (i, &(ex_lo, ex_hi)) in shape.exons.iter().enumerate() {
            if block_hi < ex_lo {
                break; // exons are sorted; no further overlap possible
            }
            if block_lo > ex_hi {
                continue;
            }
            let lo = block_lo.max(ex_lo);
            let hi = block_hi.min(ex_hi);
            if lo > hi {
                continue;
            }
            let mid_genome = (lo + hi) / 2;
            let mrna_pos = shape.cum_len[i] + (mid_genome - ex_lo);
            // Plus-strand fractional position (0..=1).
            let frac_plus = (mrna_pos as f64 + 0.5) / shape.total_exonic as f64;
            // Strand-correct: '-' strand's 5' end is at the largest mRNA
            // offset on the plus-strand coordinate.
            let frac = if shape.strand == '-' {
                1.0 - frac_plus
            } else {
                frac_plus
            };
            // Map to decile [0, 9].
            let mut idx = (frac * DECILE_COUNT as f64).floor() as i64;
            if idx < 0 {
                idx = 0;
            }
            if idx >= DECILE_COUNT as i64 {
                idx = (DECILE_COUNT - 1) as i64;
            }
            deciles[idx as usize] = deciles[idx as usize].saturating_add(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::per_gene::shape::GeneShape;

    fn empty_deciles() -> [u32; DECILE_COUNT] {
        [0; DECILE_COUNT]
    }

    #[test]
    fn plus_strand_5p_block_lands_in_first_decile() {
        // Single exon 100..=199 (100 bp). A read at the very 5' end (100..=109).
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(99, 109)], &mut d); // 0-based half-open == 1-based 100..=109
        assert_eq!(d[0], 1);
        for i in 1..DECILE_COUNT {
            assert_eq!(d[i], 0, "decile {} should be 0", i);
        }
    }

    #[test]
    fn minus_strand_reverses_decile_order() {
        // Same exon 100..=199, but strand '-'. A read at genome 5' end
        // (positions 100..=109, mRNA frac ~0.05) is at the 3' end of the
        // mRNA → decile 9.
        let shape = GeneShape::from_exons(vec![(100, 199)], '-');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(99, 109)], &mut d);
        assert_eq!(d[9], 1);
        for i in 0..(DECILE_COUNT - 1) {
            assert_eq!(d[i], 0, "decile {} should be 0", i);
        }
    }

    #[test]
    fn plus_strand_3p_end_lands_in_last_decile() {
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(189, 199)], &mut d); // 1-based 190..=199
        assert_eq!(d[9], 1);
    }

    #[test]
    fn intronic_block_is_skipped() {
        // Two exons: 100..=149 and 200..=249 (intron 150..=199).
        let shape = GeneShape::from_exons(vec![(100, 149), (200, 249)], '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(159, 179)], &mut d); // intronic 160..=179
        assert_eq!(d.iter().sum::<u32>(), 0);
    }

    #[test]
    fn block_straddling_two_exons_increments_both() {
        // Two exons: 100..=149 (mRNA 0..=49), 200..=249 (mRNA 50..=99).
        // A spliced block list [99..=149, 199..=249] simulates a read
        // spanning both exons via an N-skip (each block already exonic).
        let shape = GeneShape::from_exons(vec![(100, 149), (200, 249)], '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(99, 149), (199, 249)], &mut d);
        assert_eq!(d.iter().sum::<u32>(), 2);
        // Block 1 midpoint genome ~125 → mRNA ~25 → frac 0.25 → decile 2.
        // Block 2 midpoint genome ~225 → mRNA ~75 → frac 0.75 → decile 7.
        assert_eq!(d[2], 1);
        assert_eq!(d[7], 1);
    }

    #[test]
    fn empty_aligned_blocks_no_op() {
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[], &mut d);
        assert_eq!(d.iter().sum::<u32>(), 0);
    }

    #[test]
    fn empty_exons_no_op() {
        let shape = GeneShape::from_exons(Vec::<(u64, u64)>::new(), '+');
        let mut d = empty_deciles();
        add_block_deciles(&shape, &[(99, 109)], &mut d);
        assert_eq!(d.iter().sum::<u32>(), 0);
    }
}
