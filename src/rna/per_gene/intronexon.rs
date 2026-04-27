//! Per-gene exon vs intron aligned-bp accounting.
//!
//! Walks a read's aligned blocks once, distributing aligned bases across:
//! - bases overlapping any merged exon → `exon_aligned_bp`
//! - bases overlapping any intron of the gene → `intron_aligned_bp`
//!
//! Bases of the read that fall outside the gene's span are silently
//! dropped. The exon and intron interval sets are disjoint by
//! construction (`GeneShape::from_exons`), so a base is counted at
//! most once.

use super::shape::GeneShape;

/// Accumulator update: returns `(exon_bp_delta, intron_bp_delta)` for
/// the given aligned-block list against the gene's geometry.
///
/// `aligned_blocks` are 0-based half-open `[start, end)` reference
/// coordinates (matches `cigar_to_aligned_blocks`).
pub fn count_exon_intron_bp(shape: &GeneShape, aligned_blocks: &[(u64, u64)]) -> (u64, u64) {
    let mut exon_bp: u64 = 0;
    let mut intron_bp: u64 = 0;

    for &(b_start_0h, b_end_0h) in aligned_blocks {
        if b_end_0h <= b_start_0h {
            continue;
        }
        let block_lo = b_start_0h + 1; // → 1-based inclusive
        let block_hi = b_end_0h;

        // Sum exon overlap.
        for &(ex_lo, ex_hi) in &shape.exons {
            if block_hi < ex_lo {
                break;
            }
            if block_lo > ex_hi {
                continue;
            }
            let lo = block_lo.max(ex_lo);
            let hi = block_hi.min(ex_hi);
            if hi >= lo {
                exon_bp += hi - lo + 1;
            }
        }

        // Sum intron overlap.
        for &(in_lo, in_hi) in &shape.introns {
            if block_hi < in_lo {
                break;
            }
            if block_lo > in_hi {
                continue;
            }
            let lo = block_lo.max(in_lo);
            let hi = block_hi.min(in_hi);
            if hi >= lo {
                intron_bp += hi - lo + 1;
            }
        }
    }

    (exon_bp, intron_bp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::per_gene::shape::GeneShape;

    #[test]
    fn pure_exonic_read() {
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let (ex, intr) = count_exon_intron_bp(&shape, &[(99, 199)]); // 1-based 100..=199
        assert_eq!(ex, 100);
        assert_eq!(intr, 0);
    }

    #[test]
    fn read_spanning_one_intron() {
        // Two merged exons, one intron 150..=199.
        let shape = GeneShape::from_exons(vec![(100, 149), (200, 249)], '+');
        // Single contiguous block 130..=219 (1-based) — i.e., 129..219 0h.
        // Exonic overlap: 130..=149 (20 bp) + 200..=219 (20 bp) = 40 bp.
        // Intronic overlap: 150..=199 (50 bp).
        let (ex, intr) = count_exon_intron_bp(&shape, &[(129, 219)]);
        assert_eq!(ex, 40);
        assert_eq!(intr, 50);
    }

    #[test]
    fn intronic_only_read() {
        let shape = GeneShape::from_exons(vec![(100, 149), (200, 249)], '+');
        let (ex, intr) = count_exon_intron_bp(&shape, &[(159, 179)]); // 1-based 160..=179
        assert_eq!(ex, 0);
        assert_eq!(intr, 20);
    }

    #[test]
    fn intergenic_read_outside_gene_is_dropped() {
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let (ex, intr) = count_exon_intron_bp(&shape, &[(0, 50)]); // 1-based 1..=50
        assert_eq!(ex, 0);
        assert_eq!(intr, 0);
    }

    #[test]
    fn spliced_block_list_sums_correctly() {
        // Spliced read modeled as two M-blocks separated by N-skip.
        let shape = GeneShape::from_exons(vec![(100, 149), (200, 249)], '+');
        // Block 1 entirely in exon 1 (1-based 110..=139 → 30 bp).
        // Block 2 entirely in exon 2 (1-based 210..=239 → 30 bp).
        let (ex, intr) = count_exon_intron_bp(&shape, &[(109, 139), (209, 239)]);
        assert_eq!(ex, 60);
        assert_eq!(intr, 0);
    }

    #[test]
    fn empty_input_is_zero() {
        let shape = GeneShape::from_exons(vec![(100, 199)], '+');
        let (ex, intr) = count_exon_intron_bp(&shape, &[]);
        assert_eq!(ex, 0);
        assert_eq!(intr, 0);
    }
}
