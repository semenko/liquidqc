//! Per-gene fragmentomics accumulator (Tier 2).
//!
//! Hooks into the same per-chromosome BAM dispatcher used by
//! [`crate::rna::dupradar::counting`] and [`crate::rna::fragmentomics`]:
//! each worker owns a [`PerGeneAccumulator`] alongside its
//! [`crate::rna::rseqc::accumulators::RseqcAccumulators`] and
//! [`crate::rna::fragmentomics::FragmentomicsAccumulators`]; per-record
//! dispatch happens at the unambiguous gene-assignment branches inside
//! `process_counting_record` (single-gene single-end), the paired-end
//! same-chromosome merge, and the cross-chromosome mate-pair merge.
//!
//! Output: a sparse Parquet sibling next to the per-sample JSON envelope,
//! one row per gene that meets `--min-gene-reads`. The sample identifier
//! lives in Parquet file-level key/value metadata (not as a column), so
//! cross-sample concatenation is a pure stack with sample_id pulled from
//! file metadata.
//!
//! Submodule map:
//! - [`shape`] — per-gene merged exons / introns / cumulative exonic length
//!   (built once from the parsed GTF, shared across workers).
//! - [`state`] — `PerGeneState`: gene-indexed accumulator state (TLEN bins,
//!   end-motif 4-mer histograms, decile coverage, intron/exon bp, soft-clip
//!   counters).
//! - [`coverage`] — 5'→3' decile binning along merged-exon mRNA coords.
//! - [`intronexon`] — exon vs intron aligned-bp accounting.
//! - [`writer`] — Arrow schema + Parquet streaming writer.
//!
//! Per-gene state is `Vec<PerGeneState>` indexed by `GeneIdx` (matches the
//! `GeneIdInterner` in `dupradar::counting`). Allocation is per-worker;
//! merge is additive over fixed-length arrays. The per-worker FASTA
//! reader for end-motif lookups is borrowed from `EndMotifAccum` in the
//! parent fragmentomics module — we do not open a second `faidx::Reader`
//! handle per worker.

pub mod coverage;
pub mod intronexon;
pub mod shape;
pub mod state;
pub mod writer;

pub use shape::GeneShape;
pub use state::PerGeneState;
pub use writer::{column_names, PerGeneFileMeta, PerGeneOutput};

/// Per-worker container; one per chromosome batch. Field-wise additive
/// merge across workers; finalize filters by `--min-gene-reads`.
#[derive(Debug)]
pub struct PerGeneAccumulator {
    pub(crate) states: Vec<PerGeneState>,
}

impl PerGeneAccumulator {
    /// Construct a per-worker accumulator with zeroed state for `num_genes`
    /// genes.
    pub fn new(num_genes: usize) -> Self {
        let mut states = Vec::with_capacity(num_genes);
        states.resize_with(num_genes, PerGeneState::default);
        Self { states }
    }

    /// Merge another worker's state additively. Both accumulators must have
    /// the same `num_genes` (they were built from the same `GeneIdInterner`).
    pub fn merge(&mut self, other: PerGeneAccumulator) {
        debug_assert_eq!(self.states.len(), other.states.len());
        for (a, b) in self.states.iter_mut().zip(other.states) {
            a.merge(b);
        }
    }

    /// Increment the unambiguous-fragment-assignment counter for `gene_idx`.
    /// One call per fragment (single-end read or paired-end pair) that
    /// the counting dispatcher resolved to exactly one gene.
    pub fn contribute_fragment(&mut self, gene_idx: u32) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            state.fragment_count_assigned = state.fragment_count_assigned.saturating_add(1);
        }
    }

    /// Per primary mate: increments `primary_reads`, soft-clip flags, and
    /// the merged-exon decile + intron/exon aligned-bp accounting against
    /// the gene's geometry. Call once per primary mate (single-end: once
    /// per read; paired-end: once per mate of the pair, gated by primary
    /// mapping).
    pub fn contribute_read(
        &mut self,
        gene_idx: u32,
        gene_shape: &GeneShape,
        aligned_blocks: &[(u64, u64)],
        softclip_5p: bool,
        softclip_3p: bool,
    ) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            state.primary_reads = state.primary_reads.saturating_add(1);
            if softclip_5p {
                state.primary_reads_with_5p_clip =
                    state.primary_reads_with_5p_clip.saturating_add(1);
            }
            if softclip_3p {
                state.primary_reads_with_3p_clip =
                    state.primary_reads_with_3p_clip.saturating_add(1);
            }
            coverage::add_block_deciles(gene_shape, aligned_blocks, &mut state.decile_coverage);
            let (ex, intr, cds) = intronexon::count_exon_intron_bp(gene_shape, aligned_blocks);
            state.exon_aligned_bp = state.exon_aligned_bp.saturating_add(ex);
            state.intron_aligned_bp = state.intron_aligned_bp.saturating_add(intr);
            state.cds_aligned_bp = state.cds_aligned_bp.saturating_add(cds);
        }
    }

    /// Record a `|TLEN|` observation for the leftmost-proper-pair mate of
    /// a fragment unambiguously assigned to `gene_idx`. One call per
    /// fragment.
    pub fn contribute_tlen(&mut self, gene_idx: u32, tlen_abs: u32) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            state.record_tlen(tlen_abs);
        }
    }

    /// Record a successful end-motif lookup. `em_5p` / `em_3p` are the
    /// encoded 4-mers (the 3' value is already RC'd to fragment-strand
    /// orientation, matching the sample-level convention).
    pub fn contribute_end_motifs(&mut self, gene_idx: u32, em_5p: u8, em_3p: u8) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            let bin_5p = state.em_5p_or_alloc();
            bin_5p[em_5p as usize] = bin_5p[em_5p as usize].saturating_add(1);
            let bin_3p = state.em_3p_or_alloc();
            bin_3p[em_3p as usize] = bin_3p[em_3p as usize].saturating_add(1);
            state.em_pair_count_used = state.em_pair_count_used.saturating_add(1);
        }
    }

    /// Record an end-motif skip due to non-ACGT bases at one end.
    pub fn contribute_em_skipped_non_acgt(&mut self, gene_idx: u32) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            state.em_pair_count_skipped_non_acgt =
                state.em_pair_count_skipped_non_acgt.saturating_add(1);
        }
    }

    /// Record an end-motif skip due to fragment running off the contig
    /// (or contig missing from FASTA).
    pub fn contribute_em_skipped_oob(&mut self, gene_idx: u32) {
        if let Some(state) = self.states.get_mut(gene_idx as usize) {
            state.em_pair_count_skipped_oob = state.em_pair_count_skipped_oob.saturating_add(1);
        }
    }

    /// Dispatch per-gene contributions for a paired-end fragment that the
    /// counting dispatcher just resolved to exactly one gene. Caller
    /// supplies optional leftmost-mate metadata (TLEN/end-motif positions),
    /// per-mate primary-mapping flags + aligned blocks + soft-clip flags,
    /// and an optional borrowed FASTA reader for end-motif lookups. This
    /// helper exists so the same logic can run from
    /// `process_counting_record` (parallel & sequential) and from
    /// `ChromResult::merge` cross-chromosome reconciliation; the only
    /// load-bearing difference between those sites is whether a leftmost
    /// mate exists (cross-chrom never has one).
    #[allow(clippy::too_many_arguments)]
    pub fn dispatch_paired_fragment(
        &mut self,
        gene_idx: u32,
        gene_shape: &GeneShape,
        mate_a: PairedMateContribution<'_>,
        mate_b: PairedMateContribution<'_>,
        leftmost: Option<LeftmostContribution<'_>>,
        em_lookup: Option<&mut crate::rna::fragmentomics::EndMotifAccum>,
    ) {
        self.contribute_fragment(gene_idx);
        if mate_a.is_primary_mapped {
            self.contribute_read(
                gene_idx,
                gene_shape,
                mate_a.aligned_blocks,
                mate_a.softclip_5p,
                mate_a.softclip_3p,
            );
        }
        if mate_b.is_primary_mapped {
            self.contribute_read(
                gene_idx,
                gene_shape,
                mate_b.aligned_blocks,
                mate_b.softclip_5p,
                mate_b.softclip_3p,
            );
        }
        if let Some(lm) = leftmost {
            self.contribute_tlen(gene_idx, lm.tlen_abs);
            if let Some(em_accum) = em_lookup {
                use crate::rna::fragmentomics::EndMotifLookup;
                match em_accum.lookup_kmers_at(lm.fasta_chrom, lm.leftmost_pos, lm.rightmost_pos) {
                    EndMotifLookup::Used { kmer_5p, kmer_3p } => {
                        self.contribute_end_motifs(gene_idx, kmer_5p, kmer_3p);
                    }
                    EndMotifLookup::SkippedNonAcgt => {
                        self.contribute_em_skipped_non_acgt(gene_idx);
                    }
                    EndMotifLookup::SkippedOob => {
                        self.contribute_em_skipped_oob(gene_idx);
                    }
                }
            }
        }
    }

    /// Finalize per-worker state into Parquet rows after applying the
    /// `--min-gene-reads` filter. `genes` is the GTF-derived gene map
    /// (must be in `GeneIdInterner` insertion order — i.e., the same
    /// `IndexMap` passed to `count_reads`); `gene_counts` is the
    /// merged featureCounts result keyed by `gene_id`.
    pub fn finalize(
        self,
        min_gene_reads: u32,
        genes: &indexmap::IndexMap<String, crate::gtf::Gene>,
        gene_counts: &indexmap::IndexMap<String, crate::rna::dupradar::counting::GeneCounts>,
    ) -> PerGeneOutput {
        writer::build_rows(&self.states, min_gene_reads, genes, gene_counts)
    }
}

/// Per-mate inputs for [`PerGeneAccumulator::dispatch_paired_fragment`].
#[derive(Debug, Clone, Copy)]
pub struct PairedMateContribution<'a> {
    pub is_primary_mapped: bool,
    pub aligned_blocks: &'a [(u64, u64)],
    pub softclip_5p: bool,
    pub softclip_3p: bool,
}

/// Leftmost-proper-pair inputs for the same dispatch helper.
#[derive(Debug, Clone, Copy)]
pub struct LeftmostContribution<'a> {
    pub tlen_abs: u32,
    pub leftmost_pos: u64,
    pub rightmost_pos: u64,
    /// Raw BAM contig name — must match FASTA contig naming for end-motif lookups.
    pub fasta_chrom: &'a str,
}
