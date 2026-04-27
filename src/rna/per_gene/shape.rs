//! Per-gene merged-exon and intron geometry, precomputed once per run.
//!
//! The per-gene accumulator maps each aligned read block onto its assigned
//! gene's mRNA-exonic coordinate system (for 5'→3' decile binning) and
//! accounts exon-vs-intron aligned bp. `Gene.exons` is un-merged and
//! possibly overlapping, so we materialize a companion structure with
//! merged exons and the intron complement to keep the per-record hot path
//! free of repeated merge work.

use indexmap::IndexMap;

use crate::gtf::Gene;

/// Geometry for one gene.
#[derive(Debug, Clone)]
pub struct GeneShape {
    /// Sorted, merged exon intervals (1-based inclusive, genome coords).
    pub exons: Vec<(u64, u64)>,
    /// `cum_len[i]` = total exonic bases preceding `exons[i]` on the mRNA,
    /// i.e., `sum_{j < i} (exons[j].1 - exons[j].0 + 1)`. Length matches
    /// `exons.len()`.
    pub cum_len: Vec<u64>,
    /// Total exonic length (matches `Gene.effective_length`).
    pub total_exonic: u64,
    /// Sorted intron intervals filling the gaps between merged exons within
    /// the gene span (1-based inclusive). Length is `exons.len() - 1` for
    /// genes with multiple merged exons.
    pub introns: Vec<(u64, u64)>,
    /// Sorted, merged CDS bounding-box intervals across all transcripts of
    /// the gene (1-based inclusive). Each transcript with `cds_start /
    /// cds_end` Some contributes one bbox; bboxes are then sorted and
    /// merged. Empty for non-coding genes. Used by Phase 4 Tier 2 CDS-vs-
    /// UTR coverage; intersected with `exons` to derive CDS-aligned bp.
    pub cds_bbox: Vec<(u64, u64)>,
    /// Strand: '+', '-', or '.' (treated as '+' for decile orientation).
    pub strand: char,
}

impl GeneShape {
    /// Build the merged-exon + intron geometry from raw GTF exon list and
    /// strand. Produces an empty `cds_bbox` — Phase-4 callers that need CDS
    /// classification should use [`GeneShape::from_gene`] instead.
    #[allow(dead_code)] // Test-only convenience constructor.
    pub fn from_exons<I>(exon_iter: I, strand: char) -> Self
    where
        I: IntoIterator<Item = (u64, u64)>,
    {
        Self::from_exons_and_cds(exon_iter, std::iter::empty(), strand)
    }

    /// Build the merged-exon + intron + CDS-bbox geometry directly from the
    /// parsed [`crate::gtf::Gene`]. CDS bbox per transcript is `(cds_start,
    /// cds_end)`; transcripts without CDS are skipped.
    pub fn from_gene(gene: &crate::gtf::Gene) -> Self {
        let exon_iter = gene.exons.iter().map(|e| (e.start, e.end));
        let cds_iter = gene
            .transcripts
            .iter()
            .filter_map(|t| match (t.cds_start, t.cds_end) {
                (Some(s), Some(e)) if e >= s => Some((s, e)),
                _ => None,
            });
        Self::from_exons_and_cds(exon_iter, cds_iter, gene.strand)
    }

    /// Internal builder: merges both exon and CDS-bbox interval streams.
    pub fn from_exons_and_cds<E, C>(exon_iter: E, cds_iter: C, strand: char) -> Self
    where
        E: IntoIterator<Item = (u64, u64)>,
        C: IntoIterator<Item = (u64, u64)>,
    {
        let merged_exons = sort_and_merge(exon_iter);
        let cds_bbox = sort_and_merge(cds_iter);

        // Cumulative exonic length and total.
        let mut cum_len = Vec::with_capacity(merged_exons.len());
        let mut total: u64 = 0;
        for &(s, e) in merged_exons.iter() {
            cum_len.push(total);
            total += e - s + 1;
        }

        // Introns = gaps strictly between consecutive merged exons.
        let mut introns = Vec::with_capacity(merged_exons.len().saturating_sub(1));
        for w in merged_exons.windows(2) {
            let prev_end = w[0].1;
            let next_start = w[1].0;
            if next_start > prev_end + 1 {
                introns.push((prev_end + 1, next_start - 1));
            }
        }

        Self {
            exons: merged_exons,
            cum_len,
            total_exonic: total,
            introns,
            cds_bbox,
            strand,
        }
    }
}

fn sort_and_merge<I>(iter: I) -> Vec<(u64, u64)>
where
    I: IntoIterator<Item = (u64, u64)>,
{
    let mut intervals: Vec<(u64, u64)> = iter.into_iter().collect();
    intervals.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(intervals.len());
    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 + 1 {
                if end > last.1 {
                    last.1 = end;
                }
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

/// One `GeneShape` per gene, indexed by `GeneIdx` (insertion order in the
/// counting interner). Built once per run from the parsed GTF.
#[derive(Debug)]
pub struct GeneShapeIndex {
    pub shapes: Vec<GeneShape>,
}

impl GeneShapeIndex {
    /// Build from the parsed GTF gene map. Iteration order matches
    /// `IndexMap` insertion order, which is exactly the `GeneIdInterner`
    /// order built by `dupradar::counting::GeneIdInterner::from_genes`.
    pub fn build(genes: &IndexMap<String, Gene>) -> Self {
        let shapes: Vec<GeneShape> = genes.values().map(GeneShape::from_gene).collect();
        Self { shapes }
    }

    pub fn get(&self, idx: usize) -> Option<&GeneShape> {
        self.shapes.get(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_overlapping_and_adjacent_exons() {
        // Two overlapping exons (10-30, 25-50) merge into one (10-50).
        // Adjacent exons (52-60, 61-70) merge into one (52-70).
        // Strict gap (10-50, 52-70) becomes an intron at 51-51.
        let s = GeneShape::from_exons(vec![(25, 50), (10, 30), (61, 70), (52, 60)], '+');
        assert_eq!(s.exons, vec![(10, 50), (52, 70)]);
        assert_eq!(s.introns, vec![(51, 51)]);
        assert_eq!(s.cum_len, vec![0, 41]); // first exon length = 41 bp
        assert_eq!(s.total_exonic, 41 + 19);
    }

    #[test]
    fn no_introns_for_single_exon_gene() {
        let s = GeneShape::from_exons(vec![(100, 200)], '-');
        assert_eq!(s.exons, vec![(100, 200)]);
        assert!(s.introns.is_empty());
        assert_eq!(s.total_exonic, 101);
        assert_eq!(s.strand, '-');
    }

    #[test]
    fn introns_only_between_strict_gaps() {
        // Adjacent exons (50-60, 61-70) merge — no intron between them.
        let s = GeneShape::from_exons(vec![(50, 60), (61, 70)], '+');
        assert_eq!(s.exons, vec![(50, 70)]);
        assert!(s.introns.is_empty());
    }
}
