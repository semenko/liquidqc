//! Gene-detection saturation curve via deterministic crc32 qname bucketing.
//!
//! Each assigned read hashes to one of [`BUCKETS`] buckets per gene; at
//! finalize, fraction `f` includes buckets `[0, ceil(BUCKETS * f))`. Two
//! metrics per fraction: `genes_detected` (≥1 read) and
//! `genes_detected_min10` (≥10 reads, the standard "robustly detected"
//! threshold). Per-gene buckets are lazily boxed so cold genes cost zero.

use serde::Serialize;

/// Number of subsampling buckets. Fractions are computed at multiples of
/// `1.0 / BUCKETS`; with 100 buckets we get 1% resolution, which is more
/// than the canonical fraction grid needs.
pub const BUCKETS: usize = 100;

/// Bucket-array allocation type. Lazy per-gene `Box` so genes that are
/// never assigned a read keep zero memory.
type GeneBuckets = Box<[u32; BUCKETS]>;

#[derive(Debug, Clone)]
pub struct SaturationAccum {
    /// `Vec<Option<Box<…>>>` indexed by `GeneIdx`. Constructed at the same
    /// size as the gene set; per-gene boxes allocate on first read.
    per_gene_buckets: Vec<Option<GeneBuckets>>,
}

impl SaturationAccum {
    pub fn new(num_genes: usize) -> Self {
        Self {
            per_gene_buckets: vec![None; num_genes],
        }
    }

    pub fn observe(&mut self, gene_idx: u32, qname: &[u8]) {
        let idx = gene_idx as usize;
        if idx >= self.per_gene_buckets.len() {
            return;
        }
        let bucket = (crc32fast::hash(qname) as usize) % BUCKETS;
        let cell = self.per_gene_buckets[idx].get_or_insert_with(|| Box::new([0u32; BUCKETS]));
        cell[bucket] = cell[bucket].saturating_add(1);
    }

    pub fn merge(&mut self, other: SaturationAccum) {
        if self.per_gene_buckets.len() < other.per_gene_buckets.len() {
            self.per_gene_buckets
                .resize(other.per_gene_buckets.len(), None);
        }
        for (i, ob) in other.per_gene_buckets.into_iter().enumerate() {
            let Some(ob) = ob else {
                continue;
            };
            match self.per_gene_buckets[i] {
                Some(ref mut sb) => {
                    for k in 0..BUCKETS {
                        sb[k] = sb[k].saturating_add(ob[k]);
                    }
                }
                None => self.per_gene_buckets[i] = Some(ob),
            }
        }
    }

    /// Compute the saturation curve at the given fraction grid. Each
    /// fraction is rounded up to the nearest bucket cutoff so the
    /// last entry (`fraction = 1.0`) covers all buckets.
    pub fn finalize(self, fractions: &[f64]) -> SaturationResult {
        let mut total_assigned: u64 = 0;
        for cell in self.per_gene_buckets.iter().flatten() {
            for &c in cell.iter() {
                total_assigned = total_assigned.saturating_add(c as u64);
            }
        }
        let mut points = Vec::with_capacity(fractions.len());
        for &f in fractions {
            let f = f.clamp(0.0, 1.0);
            // Bucket cutoff: include buckets [0, c). f=1.0 → c=BUCKETS.
            let c = (f * BUCKETS as f64).ceil() as usize;
            let cutoff = c.min(BUCKETS);
            let mut genes_detected: u64 = 0;
            let mut genes_detected_min10: u64 = 0;
            let mut reads_in_subsample: u64 = 0;
            for cell in self.per_gene_buckets.iter().flatten() {
                let mut sum: u64 = 0;
                for k in 0..cutoff {
                    sum = sum.saturating_add(cell[k] as u64);
                }
                if sum > 0 {
                    genes_detected += 1;
                    reads_in_subsample = reads_in_subsample.saturating_add(sum);
                    if sum >= 10 {
                        genes_detected_min10 += 1;
                    }
                }
            }
            points.push(SaturationPoint {
                fraction: f,
                bucket_cutoff: cutoff as u32,
                reads_in_subsample,
                genes_detected,
                genes_detected_min10,
            });
        }
        SaturationResult {
            total_buckets: BUCKETS as u32,
            total_reads_assigned: total_assigned,
            points,
        }
    }
}

/// Default fraction grid. Matches the brief: 5%, 10%, 25%, 50%, 75%, 100%.
pub const DEFAULT_FRACTIONS: &[f64] = &[0.05, 0.10, 0.25, 0.50, 0.75, 1.00];

#[derive(Debug, Clone, Serialize)]
pub struct SaturationPoint {
    pub fraction: f64,
    pub bucket_cutoff: u32,
    pub reads_in_subsample: u64,
    pub genes_detected: u64,
    pub genes_detected_min10: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct SaturationResult {
    pub total_buckets: u32,
    pub total_reads_assigned: u64,
    pub points: Vec<SaturationPoint>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn observe_lazy_allocates_only_hit_genes() {
        let mut a = SaturationAccum::new(3);
        a.observe(1, b"read_a");
        assert!(a.per_gene_buckets[0].is_none());
        assert!(a.per_gene_buckets[1].is_some());
        assert!(a.per_gene_buckets[2].is_none());
    }

    #[test]
    fn same_qname_lands_in_same_bucket() {
        let mut a = SaturationAccum::new(2);
        a.observe(0, b"some_read_name");
        a.observe(0, b"some_read_name");
        let cell = a.per_gene_buckets[0].as_ref().unwrap();
        let nonzero: Vec<usize> = cell
            .iter()
            .enumerate()
            .filter_map(|(i, c)| if *c > 0 { Some(i) } else { None })
            .collect();
        assert_eq!(nonzero.len(), 1);
        assert_eq!(cell[nonzero[0]], 2);
    }

    #[test]
    fn finalize_genes_detected_is_monotone() {
        let mut a = SaturationAccum::new(50);
        // Each gene gets one read, one read per gene with a unique name so
        // buckets are spread.
        for g in 0..50u32 {
            for k in 0..20u32 {
                a.observe(g, format!("read_{g}_{k}").as_bytes());
            }
        }
        let r = a.finalize(DEFAULT_FRACTIONS);
        assert_eq!(r.points.len(), DEFAULT_FRACTIONS.len());
        // Monotone non-decreasing genes_detected as fraction grows.
        let mut prev = 0u64;
        for p in &r.points {
            assert!(p.genes_detected >= prev);
            prev = p.genes_detected;
        }
        // At fraction = 1.0, all 50 genes detected.
        assert_eq!(r.points.last().unwrap().genes_detected, 50);
    }

    #[test]
    fn finalize_genes_detected_min10_threshold() {
        let mut a = SaturationAccum::new(2);
        // Gene 0: 100 reads — should pass min10 at f=1.0.
        for k in 0..100u32 {
            a.observe(0, format!("g0_{k}").as_bytes());
        }
        // Gene 1: 5 reads — never passes min10.
        for k in 0..5u32 {
            a.observe(1, format!("g1_{k}").as_bytes());
        }
        let r = a.finalize(&[1.0]);
        assert_eq!(r.points[0].genes_detected, 2);
        assert_eq!(r.points[0].genes_detected_min10, 1);
    }

    #[test]
    fn merge_is_additive() {
        let mut a = SaturationAccum::new(3);
        let mut b = SaturationAccum::new(3);
        a.observe(0, b"r1");
        b.observe(0, b"r1");
        a.merge(b);
        let cell = a.per_gene_buckets[0].as_ref().unwrap();
        let total: u32 = cell.iter().sum();
        assert_eq!(total, 2);
    }

    #[test]
    fn fraction_zero_yields_zero_buckets() {
        let mut a = SaturationAccum::new(1);
        a.observe(0, b"r1");
        let r = a.finalize(&[0.0]);
        assert_eq!(r.points[0].bucket_cutoff, 0);
        assert_eq!(r.points[0].genes_detected, 0);
    }

    #[test]
    fn fraction_one_includes_all_buckets() {
        let mut a = SaturationAccum::new(1);
        a.observe(0, b"r1");
        let r = a.finalize(&[1.0]);
        assert_eq!(r.points[0].bucket_cutoff, BUCKETS as u32);
        assert_eq!(r.points[0].genes_detected, 1);
    }
}
