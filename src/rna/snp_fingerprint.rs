//! Genotype fingerprint at common SNP sites.
//!
//! The bundled `data/snps/common_snps_grch38_smoke.tsv` (12-site smoke panel,
//! NOT cohort-grade) is the default; override via `--snp-panel <TSV>`.
//! VCF/BCF input is reserved for a follow-up.

use anyhow::{bail, Context, Result};
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::HeaderView;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

use crate::rna::fragmentomics::common::is_primary_mapped;

const BUNDLED_SMOKE_PANEL: &str = include_str!("../../data/snps/common_snps_grch38_smoke.tsv");

/// One SNP site loaded from the panel TSV.
#[derive(Debug, Clone)]
pub struct SnpSite {
    pub chrom: String,
    /// 0-based reference position (TSV format is 1-based; converted on load).
    pub pos: u64,
    /// Reference base (uppercase ACGT).
    pub ref_base: u8,
    /// Alternate base (uppercase ACGT). Only single-ALT sites are loaded;
    /// multi-allelic rows are dropped at load with a warning.
    pub alt_base: u8,
    /// Optional dbSNP rsID for the result row.
    pub rsid: Option<String>,
}

/// Per-tid sorted SNP list, built once at sample start. Shared across
/// workers via `&` borrow.
#[derive(Debug, Default, Clone)]
pub struct SnpSiteIndex {
    /// Outer: indexed by BAM tid. Inner: SNP positions on that contig,
    /// sorted by `pos`. Sites whose contig is missing from the BAM header
    /// are dropped at build time (no entry produced).
    pub by_tid: Vec<Vec<SnpSite>>,
    /// Number of TSV rows actually loaded into the index (after skipping
    /// missing contigs / multi-allelic rows / non-ACGT alleles).
    pub sites_total: u64,
    /// Number of TSV rows whose contig is missing from the BAM header.
    pub sites_missing_contig: u64,
}

impl SnpSiteIndex {
    /// Build the index using the bundled smoke panel as the default. When
    /// `user_panel` is `Some(path)`, load that file instead — TSV format
    /// only in v1 (.vcf / .vcf.gz / .bcf return an error).
    pub fn build(header: &HeaderView, user_panel: Option<&Path>) -> Result<Self> {
        let sites = match user_panel {
            None => parse_tsv(BUNDLED_SMOKE_PANEL, "bundled:common_snps_grch38_smoke")?,
            Some(p) => {
                let ext = p
                    .extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                match ext.as_str() {
                    "tsv" | "txt" | "" => {
                        let text = std::fs::read_to_string(p).with_context(|| {
                            format!("failed to read SNP panel TSV {}", p.display())
                        })?;
                        parse_tsv(&text, &format!("file:{}", p.display()))?
                    }
                    "vcf" | "gz" | "bcf" => bail!(
                        "VCF/BCF SNP panels are reserved for a follow-up release; \
                         convert {} to TSV (chrom<TAB>pos<TAB>ref<TAB>alt[<TAB>rsid]) \
                         and pass that via `--snp-panel` for v1.",
                        p.display()
                    ),
                    _ => bail!(
                        "unrecognized SNP panel extension `{ext}` for {}; expected .tsv",
                        p.display()
                    ),
                }
            }
        };

        // Resolve contig name → tid against the BAM header.
        let mut name_to_tid: HashMap<String, i32> = HashMap::new();
        for tid in 0..header.target_count() {
            let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            name_to_tid.insert(name, tid as i32);
        }

        let mut by_tid: Vec<Vec<SnpSite>> =
            (0..header.target_count()).map(|_| Vec::new()).collect();
        let mut missing: u64 = 0;
        let mut loaded: u64 = 0;
        for site in sites {
            match name_to_tid.get(site.chrom.as_str()) {
                Some(&tid) => {
                    by_tid[tid as usize].push(site);
                    loaded += 1;
                }
                None => missing += 1,
            }
        }
        for v in by_tid.iter_mut() {
            v.sort_by_key(|s| s.pos);
        }
        Ok(Self {
            by_tid,
            sites_total: loaded,
            sites_missing_contig: missing,
        })
    }

    /// True when no sites resolved against the BAM header — fingerprint is
    /// silently skipped at finalize.
    pub fn is_empty(&self) -> bool {
        self.sites_total == 0
    }
}

fn parse_tsv(text: &str, source_label: &str) -> Result<Vec<SnpSite>> {
    let mut out = Vec::new();
    let mut header_seen = false;
    for (i, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if !(4..=5).contains(&fields.len()) {
            bail!(
                "{source_label}: line {} has {} tab-separated fields, expected 4 or 5 \
                 (chrom<TAB>pos<TAB>ref<TAB>alt[<TAB>rsid])",
                i + 1,
                fields.len()
            );
        }
        if !header_seen
            && fields[0].eq_ignore_ascii_case("chrom")
            && fields[1].eq_ignore_ascii_case("pos")
        {
            header_seen = true;
            continue;
        }
        header_seen = true;
        let chrom = fields[0].to_string();
        let pos_1b: u64 = fields[1].parse().with_context(|| {
            format!(
                "{source_label}: line {}: pos `{}` is not an integer",
                i + 1,
                fields[1]
            )
        })?;
        if pos_1b == 0 {
            bail!("{source_label}: line {}: pos must be 1-based (>= 1)", i + 1);
        }
        let ref_base = parse_single_base(fields[2]).with_context(|| {
            format!(
                "{source_label}: line {}: REF must be a single ACGT base",
                i + 1
            )
        })?;
        let alt_base = parse_single_base(fields[3]).with_context(|| {
            format!(
                "{source_label}: line {}: ALT must be a single ACGT base",
                i + 1
            )
        })?;
        let rsid = fields.get(4).map(|s| s.to_string());
        out.push(SnpSite {
            chrom,
            pos: pos_1b - 1,
            ref_base,
            alt_base,
            rsid,
        });
    }
    Ok(out)
}

fn parse_single_base(s: &str) -> Result<u8> {
    let bytes = s.as_bytes();
    if bytes.len() != 1 {
        bail!("expected single base, got `{}`", s);
    }
    match bytes[0].to_ascii_uppercase() {
        b @ (b'A' | b'C' | b'G' | b'T') => Ok(b),
        other => bail!("non-ACGT base `{}`", other as char),
    }
}

/// Per-worker pileup accumulator. Counts ACGT base calls at each SNP site
/// covered by the worker's BAM partition.
#[derive(Debug, Default, Clone)]
pub struct SnpFingerprintAccum {
    /// `(tid, pos)` → `[A, C, G, T]` counts. Lazy: only sites that get at
    /// least one base contribute an entry.
    counts: HashMap<(i32, u64), [u32; 4]>,
}

impl SnpFingerprintAccum {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn process_read(&mut self, record: &Record, index: &SnpSiteIndex, mapq_cut: u8) {
        if !is_primary_mapped(record, mapq_cut) {
            return;
        }
        let tid = record.tid();
        if tid < 0 {
            return;
        }
        let sites = match index.by_tid.get(tid as usize) {
            Some(v) if !v.is_empty() => v,
            _ => return,
        };
        // Bound the read's reference span so we can binary-search the
        // SNP list to the relevant slice.
        let cigar = record.cigar();
        let ref_start = record.pos() as u64;
        let ref_end = cigar.end_pos() as u64; // half-open
        if ref_end <= ref_start {
            return;
        }
        // First and last candidate SNP indices.
        let first = sites.partition_point(|s| s.pos < ref_start);
        let last = sites.partition_point(|s| s.pos < ref_end);
        if first >= last {
            return;
        }
        let seq = record.seq();
        for site in &sites[first..last] {
            let Some(base) = read_base_at_ref_pos(&cigar, ref_start, &seq, site.pos) else {
                continue;
            };
            let Some(idx) = base_index(base) else {
                continue;
            };
            let cell = self.counts.entry((tid, site.pos)).or_insert([0u32; 4]);
            cell[idx] = cell[idx].saturating_add(1);
        }
    }

    pub fn merge(&mut self, other: SnpFingerprintAccum) {
        for (k, ov) in other.counts {
            let cell = self.counts.entry(k).or_insert([0u32; 4]);
            for i in 0..4 {
                cell[i] = cell[i].saturating_add(ov[i]);
            }
        }
    }

    pub fn finalize(self, index: &SnpSiteIndex) -> SnpFingerprintResult {
        let mut entries = Vec::new();
        let mut sites_called: u64 = 0;
        let mut sites_with_any_depth: u64 = 0;
        for (tid, sites) in index.by_tid.iter().enumerate() {
            for site in sites {
                let counts = self
                    .counts
                    .get(&(tid as i32, site.pos))
                    .copied()
                    .unwrap_or([0u32; 4]);
                let depth: u64 = counts.iter().map(|c| *c as u64).sum::<u64>();
                let ref_idx = base_index(site.ref_base);
                let alt_idx = base_index(site.alt_base);
                let ref_count = ref_idx.map(|i| counts[i]).unwrap_or(0) as u64;
                let alt_count = alt_idx.map(|i| counts[i]).unwrap_or(0) as u64;
                let other_count = depth.saturating_sub(ref_count).saturating_sub(alt_count);
                let alt_fraction = if depth > 0 {
                    alt_count as f64 / depth as f64
                } else {
                    0.0
                };
                if depth > 0 {
                    sites_with_any_depth += 1;
                    if depth >= 4 {
                        sites_called += 1;
                    }
                }
                entries.push(SnpFingerprintSite {
                    chrom: site.chrom.clone(),
                    pos_1based: site.pos + 1,
                    ref_base: (site.ref_base as char).to_string(),
                    alt_base: (site.alt_base as char).to_string(),
                    rsid: site.rsid.clone(),
                    depth,
                    ref_count,
                    alt_count,
                    other_count,
                    alt_fraction,
                });
            }
        }
        SnpFingerprintResult {
            sites_total: index.sites_total,
            sites_missing_contig: index.sites_missing_contig,
            sites_with_any_depth,
            sites_called,
            min_depth_for_called: 4,
            sites: entries,
        }
    }
}

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Walk the CIGAR to locate the read base at a target reference position.
/// `cigar` is the raw CIGAR view; `ref_start` is the record's mapping
/// position (`record.pos() as u64`); `seq` is `record.seq()`. Returns
/// `None` if the position falls in a deletion / splice gap / soft-clip.
fn read_base_at_ref_pos(
    cigar: &rust_htslib::bam::record::CigarStringView,
    ref_start: u64,
    seq: &rust_htslib::bam::record::Seq,
    target_ref_pos: u64,
) -> Option<u8> {
    let mut ref_pos = ref_start;
    let mut read_pos: usize = 0;
    for op in cigar.iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = *n as u64;
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + n {
                    let offset = (target_ref_pos - ref_pos) as usize;
                    return Some(seq[read_pos + offset]);
                }
                ref_pos += n;
                read_pos += n as usize;
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                read_pos += *n as usize;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                ref_pos += *n as u64;
                if target_ref_pos < ref_pos {
                    return None;
                }
            }
            // HardClip / Pad consume neither.
            _ => {}
        }
    }
    None
}

#[derive(Debug, Clone, Serialize)]
pub struct SnpFingerprintSite {
    pub chrom: String,
    pub pos_1based: u64,
    #[serde(rename = "ref")]
    pub ref_base: String,
    #[serde(rename = "alt")]
    pub alt_base: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rsid: Option<String>,
    pub depth: u64,
    pub ref_count: u64,
    pub alt_count: u64,
    pub other_count: u64,
    pub alt_fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct SnpFingerprintResult {
    /// Number of SNP sites whose contig resolved against the BAM header.
    pub sites_total: u64,
    /// SNP TSV rows whose contig was missing from the BAM header.
    pub sites_missing_contig: u64,
    /// Sites with `depth > 0`.
    pub sites_with_any_depth: u64,
    /// Sites with `depth >= min_depth_for_called`.
    pub sites_called: u64,
    /// Depth threshold for `sites_called`.
    pub min_depth_for_called: u32,
    /// Per-site rows, in the same order as the bundled / user-supplied panel
    /// (per-tid order from the index, then by pos within tid).
    pub sites: Vec<SnpFingerprintSite>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tsv_accepts_header_and_optional_rsid() {
        let text = "chrom\tpos\tref\talt\trsid\nchr1\t100\tA\tG\trs1\nchr1\t200\tC\tT\n";
        let sites = parse_tsv(text, "test").unwrap();
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].pos, 99);
        assert_eq!(sites[0].ref_base, b'A');
        assert_eq!(sites[0].alt_base, b'G');
        assert_eq!(sites[0].rsid.as_deref(), Some("rs1"));
        assert!(sites[1].rsid.is_none());
    }

    #[test]
    fn parse_tsv_rejects_zero_pos() {
        let text = "chr1\t0\tA\tG\n";
        assert!(parse_tsv(text, "test").is_err());
    }

    #[test]
    fn parse_tsv_rejects_non_acgt() {
        let text = "chr1\t100\tA\tN\n";
        assert!(parse_tsv(text, "test").is_err());
    }

    #[test]
    fn bundled_smoke_panel_parses_cleanly() {
        let sites = parse_tsv(BUNDLED_SMOKE_PANEL, "bundled:smoke").unwrap();
        assert!(!sites.is_empty());
        for s in &sites {
            assert!(matches!(s.ref_base, b'A' | b'C' | b'G' | b'T'));
            assert!(matches!(s.alt_base, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn merge_is_additive() {
        let mut a = SnpFingerprintAccum::new();
        let mut b = SnpFingerprintAccum::new();
        a.counts.insert((0, 100), [1, 0, 0, 0]);
        b.counts.insert((0, 100), [2, 0, 0, 0]);
        b.counts.insert((0, 200), [0, 5, 0, 0]);
        a.merge(b);
        assert_eq!(a.counts[&(0, 100)][0], 3);
        assert_eq!(a.counts[&(0, 200)][1], 5);
    }
}
