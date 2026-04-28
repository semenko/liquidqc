//! Hemoglobin / ribosomal-protein / apolipoprotein read fractions.
//!
//! Three bundled HGNC symbol lists (`data/gene_classes/`) are resolved against
//! `gene_name` once at GTF load time; finalize sums `fc_reads` per class.
//! Genes without a `gene_name` attribute never match (no fuzzy ENSG mapping).

use indexmap::IndexMap;
use serde::Serialize;
use std::collections::HashSet;

use crate::gtf::Gene;
use crate::rna::dupradar::counting::GeneCounts;

/// Bundled symbol list for hemoglobin subunit genes (HBA / HBB / HBG / HBM /
/// HBQ / HBZ). One symbol per line, `#` lines are comments.
const HEMOGLOBIN_LIST: &str = include_str!("../../data/gene_classes/hemoglobin.txt");

/// Bundled symbol list for ribosomal protein genes (HGNC RPS / RPL / MRPS /
/// MRPL families).
const RIBOSOMAL_PROTEIN_LIST: &str = include_str!("../../data/gene_classes/ribosomal_protein.txt");

/// Bundled symbol list for apolipoprotein genes (HGNC APO* family).
const APOLIPOPROTEIN_LIST: &str = include_str!("../../data/gene_classes/apolipoprotein.txt");

/// Pre-resolved per-class membership, keyed by gene_id. Built once at GTF
/// load time; consumed at envelope-build time when summing per-class reads.
#[derive(Debug, Clone, Default)]
pub struct GeneClassIndex {
    pub hemoglobin: HashSet<String>,
    pub ribosomal_protein: HashSet<String>,
    pub apolipoprotein: HashSet<String>,
    /// Number of bundled symbols for each class — surfaced in the result so
    /// consumers can sanity-check that the GTF actually carried the expected
    /// gene_name annotations (e.g., zero matches against a 10-gene list is a
    /// strong sign the GTF is missing the `gene_name` attribute).
    pub bundled_counts: GeneClassBundledCounts,
}

#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct GeneClassBundledCounts {
    pub hemoglobin: u64,
    pub ribosomal_protein: u64,
    pub apolipoprotein: u64,
}

impl GeneClassIndex {
    /// Resolve bundled symbols against the GTF gene set. `genes` is the
    /// `IndexMap<gene_id, Gene>` returned by [`crate::gtf::parse_gtf`]; the
    /// `gene_name` attribute is read from `gene.attributes["gene_name"]`.
    pub fn build(genes: &IndexMap<String, Gene>) -> Self {
        let hb_symbols = parse_symbol_list(HEMOGLOBIN_LIST);
        let rp_symbols = parse_symbol_list(RIBOSOMAL_PROTEIN_LIST);
        let apo_symbols = parse_symbol_list(APOLIPOPROTEIN_LIST);

        let bundled_counts = GeneClassBundledCounts {
            hemoglobin: hb_symbols.len() as u64,
            ribosomal_protein: rp_symbols.len() as u64,
            apolipoprotein: apo_symbols.len() as u64,
        };

        let mut hemoglobin = HashSet::new();
        let mut ribosomal_protein = HashSet::new();
        let mut apolipoprotein = HashSet::new();

        for (gene_id, gene) in genes {
            let Some(symbol) = gene.attributes.get("gene_name") else {
                continue;
            };
            if hb_symbols.contains(symbol.as_str()) {
                hemoglobin.insert(gene_id.clone());
            }
            if rp_symbols.contains(symbol.as_str()) {
                ribosomal_protein.insert(gene_id.clone());
            }
            if apo_symbols.contains(symbol.as_str()) {
                apolipoprotein.insert(gene_id.clone());
            }
        }
        Self {
            hemoglobin,
            ribosomal_protein,
            apolipoprotein,
            bundled_counts,
        }
    }

    /// Compute per-class read totals + fractions of the supplied denominator.
    /// `gene_counts` is the same `IndexMap<gene_id, GeneCounts>` carried by
    /// [`crate::rna::dupradar::counting::CountResult`]. Use
    /// `count_result.fc_assigned` as the denominator (matches the
    /// `featureCounts -g gene_id` assignment population).
    pub fn finalize(
        &self,
        gene_counts: &IndexMap<String, GeneCounts>,
        denominator: u64,
    ) -> GeneClassFractionsResult {
        let sum = |members: &HashSet<String>| -> u64 {
            let mut total: u64 = 0;
            for gid in members {
                if let Some(c) = gene_counts.get(gid) {
                    total = total.saturating_add(c.fc_reads);
                }
            }
            total
        };
        let hb_reads = sum(&self.hemoglobin);
        let rp_reads = sum(&self.ribosomal_protein);
        let apo_reads = sum(&self.apolipoprotein);
        use crate::rna::safe_fraction;
        GeneClassFractionsResult {
            hemoglobin: GeneClassEntry {
                bundled_symbols: self.bundled_counts.hemoglobin,
                matched_genes: self.hemoglobin.len() as u64,
                reads: hb_reads,
                fraction: safe_fraction(hb_reads, denominator),
            },
            ribosomal_protein: GeneClassEntry {
                bundled_symbols: self.bundled_counts.ribosomal_protein,
                matched_genes: self.ribosomal_protein.len() as u64,
                reads: rp_reads,
                fraction: safe_fraction(rp_reads, denominator),
            },
            apolipoprotein: GeneClassEntry {
                bundled_symbols: self.bundled_counts.apolipoprotein,
                matched_genes: self.apolipoprotein.len() as u64,
                reads: apo_reads,
                fraction: safe_fraction(apo_reads, denominator),
            },
            denominator,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct GeneClassEntry {
    /// Number of bundled symbols for this class (e.g. 10 for hemoglobin).
    pub bundled_symbols: u64,
    /// Number of GTF genes that resolved to a bundled symbol.
    pub matched_genes: u64,
    /// Sum of `fc_reads` across the matched genes.
    pub reads: u64,
    /// `reads / denominator`. 0.0 when the denominator is 0.
    pub fraction: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct GeneClassFractionsResult {
    pub hemoglobin: GeneClassEntry,
    pub ribosomal_protein: GeneClassEntry,
    pub apolipoprotein: GeneClassEntry,
    /// Denominator used for the `fraction` fields (typically
    /// `count_result.fc_assigned`).
    pub denominator: u64,
}

fn parse_symbol_list(text: &str) -> HashSet<&str> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn gene(name: &str) -> Gene {
        let mut attrs = HashMap::new();
        attrs.insert("gene_name".to_string(), name.to_string());
        Gene {
            gene_id: format!("ensg_{name}"),
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
    fn hemoglobin_list_is_non_empty_after_parse() {
        let symbols = parse_symbol_list(HEMOGLOBIN_LIST);
        assert!(symbols.contains("HBA1"));
        assert!(symbols.contains("HBB"));
        assert_eq!(symbols.len(), 10);
    }

    #[test]
    fn ribosomal_protein_list_includes_rps_rpl_mrps_mrpl() {
        let symbols = parse_symbol_list(RIBOSOMAL_PROTEIN_LIST);
        assert!(symbols.contains("RPS4Y1")); // Y-chromosome RPS — used by sex inference too
        assert!(symbols.contains("RPL10"));
        assert!(symbols.contains("MRPS2"));
        assert!(symbols.contains("MRPL1"));
    }

    #[test]
    fn apolipoprotein_list_excludes_apoptosis_genes() {
        let symbols = parse_symbol_list(APOLIPOPROTEIN_LIST);
        assert!(symbols.contains("APOE"));
        assert!(symbols.contains("APOB"));
        assert!(!symbols.contains("APOPT1"));
    }

    #[test]
    fn build_resolves_symbols_to_gene_ids() {
        let mut genes = IndexMap::new();
        genes.insert("g1".to_string(), gene("HBA1"));
        genes.insert("g2".to_string(), gene("RPL10"));
        genes.insert("g3".to_string(), gene("APOE"));
        genes.insert("g4".to_string(), gene("FOOBAR")); // unrelated
        let idx = GeneClassIndex::build(&genes);
        assert!(idx.hemoglobin.contains("g1"));
        assert!(idx.ribosomal_protein.contains("g2"));
        assert!(idx.apolipoprotein.contains("g3"));
        assert!(!idx.hemoglobin.contains("g4"));
    }

    #[test]
    fn build_skips_genes_without_gene_name_attribute() {
        let mut genes = IndexMap::new();
        let mut nameless = gene("HBA1");
        nameless.attributes.clear(); // strip gene_name
        genes.insert("g1".to_string(), nameless);
        let idx = GeneClassIndex::build(&genes);
        assert!(idx.hemoglobin.is_empty());
    }

    #[test]
    fn finalize_sums_fc_reads_and_divides_by_denominator() {
        let mut genes = IndexMap::new();
        genes.insert("hb1".to_string(), gene("HBB"));
        genes.insert("rp1".to_string(), gene("RPL10"));
        let idx = GeneClassIndex::build(&genes);

        let mut counts: IndexMap<String, GeneCounts> = IndexMap::new();
        counts.insert(
            "hb1".to_string(),
            GeneCounts {
                fc_reads: 100,
                ..Default::default()
            },
        );
        counts.insert(
            "rp1".to_string(),
            GeneCounts {
                fc_reads: 50,
                ..Default::default()
            },
        );

        let r = idx.finalize(&counts, 1000);
        assert_eq!(r.hemoglobin.reads, 100);
        assert_eq!(r.hemoglobin.matched_genes, 1);
        assert_eq!(r.hemoglobin.bundled_symbols, 10);
        assert!((r.hemoglobin.fraction - 0.1).abs() < 1e-12);
        assert_eq!(r.ribosomal_protein.reads, 50);
        assert!((r.ribosomal_protein.fraction - 0.05).abs() < 1e-12);
        assert_eq!(r.apolipoprotein.reads, 0);
        assert_eq!(r.apolipoprotein.fraction, 0.0);
    }

    #[test]
    fn finalize_zero_denominator_yields_zero_fractions() {
        let genes: IndexMap<String, Gene> = IndexMap::new();
        let idx = GeneClassIndex::build(&genes);
        let counts: IndexMap<String, GeneCounts> = IndexMap::new();
        let r = idx.finalize(&counts, 0);
        assert_eq!(r.hemoglobin.fraction, 0.0);
        assert_eq!(r.ribosomal_protein.fraction, 0.0);
        assert_eq!(r.apolipoprotein.fraction, 0.0);
    }
}
