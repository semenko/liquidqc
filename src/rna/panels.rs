//! Marker-panel read aggregation.
//!
//! 5-column TSVs (`gene_id<TAB>gene_symbol<TAB>panel<TAB>cell_type<TAB>weight`)
//! resolve each gene to one or more `(panel, cell_type)` cells; finalize sums
//! `fc_reads` per cell. Three bundled panels (`lm22`, `tabula_sapiens_cfrna`,
//! `vorperian`) are partial curations with citations in the file headers and
//! `NOTICE`. `--panels` selects bundled subsets; `--panels-tsv` adds user TSVs.
//! Bundled rows match `gene_symbol` against GTF `gene_name`; user rows match
//! by `gene_id` first then `gene_symbol`.

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use crate::gtf::Gene;
use crate::rna::dupradar::counting::GeneCounts;

/// Bundled lm22 (immune cell) panel TSV.
const LM22_TSV: &str = include_str!("../../data/panels/lm22.tsv");
/// Bundled Tabula Sapiens cfRNA-relevant panel TSV.
const TABULA_SAPIENS_CFRNA_TSV: &str = include_str!("../../data/panels/tabula_sapiens_cfrna.tsv");
/// Bundled Vorperian cfRNA cell-of-origin panel TSV.
const VORPERIAN_TSV: &str = include_str!("../../data/panels/vorperian.tsv");

/// Identifier for each bundled panel — value of `--panels` CSV entries.
pub const BUNDLED_PANEL_NAMES: &[&str] = &["lm22", "tabula_sapiens_cfrna", "vorperian"];

#[derive(Debug, Clone)]
struct PanelEntry {
    /// Optional Ensembl/GENCODE/etc. gene_id; matched first when present.
    gene_id: Option<String>,
    /// Required gene_symbol (used as a fallback to gene_id, and as the row's
    /// human-readable identifier in error messages).
    gene_symbol: String,
    /// Panel identifier (e.g. `"lm22"`).
    panel: String,
    /// Cell-type label within the panel.
    cell_type: String,
}

/// Resolved panel index. Built once at startup; consumed by [`aggregate`].
#[derive(Debug, Clone, Default)]
pub struct PanelIndex {
    /// Resolved gene_id → [(panel, cell_type)] memberships.
    by_gene_id: HashMap<String, Vec<(String, String)>>,
    /// Set of (panel, cell_type) pairs known to the index. Used to ensure
    /// every cell-type appears in the result even with zero reads.
    cells: HashSet<(String, String)>,
    /// One pre-resolved entry per loaded TSV row, including any rows that
    /// did not match a GTF gene. Surfaced for `panel_summary` so consumers
    /// can detect curation drift between the TSV and the input GTF.
    loaded_entries: u64,
    matched_entries: u64,
}

impl PanelIndex {
    /// Build the index from the user-selected bundled panels, optionally
    /// extended with one or more user-supplied TSV paths.
    ///
    /// `bundled_selection` is the list of panel names from `--panels`. If
    /// `None`, all bundled panels load. An empty `Some(vec![])` disables
    /// bundled panels entirely.
    pub fn build(
        genes: &IndexMap<String, Gene>,
        bundled_selection: Option<&[&str]>,
        user_tsv_paths: &[&Path],
    ) -> Result<Self> {
        let mut entries: Vec<PanelEntry> = Vec::new();

        let bundled_names: Vec<&str> = match bundled_selection {
            Some(s) => s.to_vec(),
            None => BUNDLED_PANEL_NAMES.to_vec(),
        };
        for name in &bundled_names {
            let text = match *name {
                "lm22" => LM22_TSV,
                "tabula_sapiens_cfrna" => TABULA_SAPIENS_CFRNA_TSV,
                "vorperian" => VORPERIAN_TSV,
                other => bail!(
                    "unknown bundled panel `{other}`; supported: {:?}",
                    BUNDLED_PANEL_NAMES
                ),
            };
            parse_tsv(text, &format!("bundled:{name}"), &mut entries)?;
        }
        for path in user_tsv_paths {
            let text = std::fs::read_to_string(path)
                .with_context(|| format!("failed to read panel TSV {}", path.display()))?;
            parse_tsv(&text, &format!("file:{}", path.display()), &mut entries)?;
        }

        let loaded = entries.len() as u64;

        // Build symbol → gene_id index from the GTF (for entries that supply
        // gene_symbol but not gene_id).
        let mut symbol_to_id: HashMap<&str, &str> = HashMap::new();
        for (gid, gene) in genes {
            if let Some(sym) = gene.attributes.get("gene_name") {
                symbol_to_id.insert(sym.as_str(), gid.as_str());
            }
        }

        let mut by_gene_id: HashMap<String, Vec<(String, String)>> = HashMap::new();
        let mut cells: HashSet<(String, String)> = HashSet::new();
        let mut matched: u64 = 0;
        for entry in entries {
            cells.insert((entry.panel.clone(), entry.cell_type.clone()));
            // Resolve gene_id either from the TSV directly or via the
            // gene_symbol → gene_id lookup against the GTF.
            let resolved_id: Option<String> = match entry.gene_id.as_deref() {
                Some(gid) if genes.contains_key(gid) => Some(gid.to_string()),
                _ => symbol_to_id
                    .get(entry.gene_symbol.as_str())
                    .map(|s| s.to_string()),
            };
            if let Some(gid) = resolved_id {
                by_gene_id
                    .entry(gid)
                    .or_default()
                    .push((entry.panel, entry.cell_type));
                matched += 1;
            }
        }

        Ok(Self {
            by_gene_id,
            cells,
            loaded_entries: loaded,
            matched_entries: matched,
        })
    }

    /// True if no bundled or user panels resolved to any GTF gene.
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Aggregate per-gene `featureCounts` reads into per-(panel, cell_type)
    /// totals + fractions of `denominator` (typically `fc_assigned`).
    ///
    /// Cell-types known to the index but with zero matched reads are still
    /// emitted with `reads = 0, fraction = 0.0` so cross-sample tables stay
    /// schema-stable.
    pub fn aggregate(
        &self,
        gene_counts: &IndexMap<String, GeneCounts>,
        denominator: u64,
    ) -> PanelsResult {
        let mut sums: HashMap<(String, String), u64> = HashMap::new();
        // Initialize known cells to zero so missing rows surface explicitly.
        for cell in &self.cells {
            sums.insert(cell.clone(), 0);
        }

        for (gid, members) in &self.by_gene_id {
            let reads = match gene_counts.get(gid) {
                Some(c) => c.fc_reads,
                None => continue,
            };
            for cell in members {
                if let Some(entry) = sums.get_mut(cell) {
                    *entry = entry.saturating_add(reads);
                }
            }
        }

        // Group results by panel for stable, schema-friendly output.
        let mut panels: BTreeMap<String, BTreeMap<String, PanelCellType>> = BTreeMap::new();
        for ((panel, cell_type), reads) in sums.into_iter() {
            let fraction = if denominator == 0 {
                0.0
            } else {
                reads as f64 / denominator as f64
            };
            panels
                .entry(panel)
                .or_default()
                .insert(cell_type, PanelCellType { reads, fraction });
        }

        PanelsResult {
            denominator,
            loaded_entries: self.loaded_entries,
            matched_entries: self.matched_entries,
            panels,
        }
    }
}

fn parse_tsv(text: &str, source_label: &str, into: &mut Vec<PanelEntry>) -> Result<()> {
    let mut header_seen = false;
    for (lineno_minus_one, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 5 {
            bail!(
                "{source_label}: line {} has {} tab-separated fields, expected 5 \
                 (gene_id<TAB>gene_symbol<TAB>panel<TAB>cell_type<TAB>weight)",
                lineno_minus_one + 1,
                fields.len()
            );
        }
        // First non-comment row is the header; skip it.
        if !header_seen {
            header_seen = true;
            if fields[0].eq_ignore_ascii_case("gene_id")
                && fields[1].eq_ignore_ascii_case("gene_symbol")
            {
                continue;
            }
            // No header row was supplied; fall through and treat this as a
            // data row.
        }
        let gene_id = if fields[0].is_empty() {
            None
        } else {
            Some(fields[0].to_string())
        };
        if fields[1].is_empty() {
            bail!(
                "{source_label}: line {} has empty gene_symbol",
                lineno_minus_one + 1
            );
        }
        let panel = fields[2].to_string();
        let cell_type = fields[3].to_string();
        if panel.is_empty() || cell_type.is_empty() {
            bail!(
                "{source_label}: line {} has empty panel or cell_type",
                lineno_minus_one + 1
            );
        }
        into.push(PanelEntry {
            gene_id,
            gene_symbol: fields[1].to_string(),
            panel,
            cell_type,
        });
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct PanelCellType {
    pub reads: u64,
    pub fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct PanelsResult {
    /// Total reads used as the denominator for `panels.<name>.<cell>.fraction`.
    pub denominator: u64,
    /// Total panel-TSV rows loaded across all selected sources.
    pub loaded_entries: u64,
    /// Subset of `loaded_entries` that resolved to a GTF gene.
    pub matched_entries: u64,
    /// Per-panel, per-cell-type tallies.
    pub panels: BTreeMap<String, BTreeMap<String, PanelCellType>>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn gene(symbol: &str) -> Gene {
        let mut attrs = HashMap::new();
        attrs.insert("gene_name".to_string(), symbol.to_string());
        Gene {
            gene_id: format!("ensg_{symbol}"),
            chrom: "chr1".to_string(),
            start: 0,
            end: 100,
            strand: '+',
            exons: vec![],
            effective_length: 0,
            attributes: attrs,
            transcripts: vec![],
        }
    }

    #[test]
    fn bundled_panels_parse_cleanly() {
        let mut entries = Vec::new();
        parse_tsv(LM22_TSV, "bundled:lm22", &mut entries).unwrap();
        assert!(!entries.is_empty());
        assert!(entries.iter().any(|e| e.panel == "lm22"));
        let mut entries2 = Vec::new();
        parse_tsv(
            TABULA_SAPIENS_CFRNA_TSV,
            "bundled:tabula_sapiens_cfrna",
            &mut entries2,
        )
        .unwrap();
        assert!(!entries2.is_empty());
        let mut entries3 = Vec::new();
        parse_tsv(VORPERIAN_TSV, "bundled:vorperian", &mut entries3).unwrap();
        assert!(!entries3.is_empty());
    }

    #[test]
    fn build_resolves_bundled_symbols_against_gtf() {
        // Populate a GTF that carries CD3D + ALB (one symbol from lm22, one
        // from tabula_sapiens_cfrna).
        let mut genes = IndexMap::new();
        genes.insert("ensg_CD3D".to_string(), gene("CD3D"));
        genes.insert("ensg_ALB".to_string(), gene("ALB"));
        let idx = PanelIndex::build(&genes, None, &[]).unwrap();
        assert!(idx.matched_entries > 0);
        // CD3D should map to lm22:T_cells in this curated subset.
        let mut found_lm22_t = false;
        for members in idx.by_gene_id.values() {
            for (panel, cell) in members {
                if panel == "lm22" && cell == "T_cells" {
                    found_lm22_t = true;
                }
            }
        }
        assert!(found_lm22_t);
    }

    #[test]
    fn aggregate_sums_reads_and_emits_zeroes_for_unmatched_cells() {
        let mut genes = IndexMap::new();
        genes.insert("ensg_CD3D".to_string(), gene("CD3D"));
        let idx = PanelIndex::build(&genes, Some(&["lm22"]), &[]).unwrap();

        let mut counts: IndexMap<String, GeneCounts> = IndexMap::new();
        counts.insert(
            "ensg_CD3D".to_string(),
            GeneCounts {
                fc_reads: 42,
                ..Default::default()
            },
        );

        let r = idx.aggregate(&counts, 1000);
        let lm22 = r.panels.get("lm22").unwrap();
        let t_cells = lm22.get("T_cells").unwrap();
        assert_eq!(t_cells.reads, 42);
        // T_cells_CD8 is in the bundled lm22 cell-type set even though the GTF
        // didn't carry CD8A/CD8B. Should still be present with reads=0.
        assert!(lm22.contains_key("T_cells_CD8"));
        assert_eq!(lm22.get("T_cells_CD8").unwrap().reads, 0);
    }

    #[test]
    fn build_rejects_unknown_bundled_panel() {
        let genes: IndexMap<String, Gene> = IndexMap::new();
        let r = PanelIndex::build(&genes, Some(&["does_not_exist"]), &[]);
        assert!(r.is_err());
    }

    #[test]
    fn malformed_tsv_is_rejected() {
        let mut entries = Vec::new();
        let r = parse_tsv("two\tcols", "test:bad", &mut entries);
        assert!(r.is_err());
    }

    #[test]
    fn empty_panels_selection_disables_bundled() {
        let genes: IndexMap<String, Gene> = IndexMap::new();
        let idx = PanelIndex::build(&genes, Some(&[]), &[]).unwrap();
        assert!(idx.is_empty());
        assert_eq!(idx.matched_entries, 0);
        assert_eq!(idx.loaded_entries, 0);
    }
}
