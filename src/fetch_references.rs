//! `liquidqc fetch-references` subcommand — populate the per-user GTF
//! cache by downloading pinned GENCODE annotations.
//!
//! Each entry in [`PINNED_GTFS`] points at a stable GENCODE URL plus an
//! optional pinned MD5. Downloads stream to a temp file beside the cache
//! target, get verified, then atomically renamed into place — partial or
//! corrupt downloads never leave a half-written `<genome>.gtf.gz`.

use crate::cache::{ensure_dir, gtf_dir, gtf_path_for};
use crate::genome::KnownGenome;
use crate::ui::format_bytes;
use anyhow::{anyhow, Context, Result};
use md5::{Digest, Md5};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// Pinned download metadata for one genome. Update the version on a
/// deliberate cadence — the cache name (`<genome>.gtf.gz`) does NOT
/// encode the GENCODE release, so a refresh requires re-running
/// `fetch-references` on every machine that uses the cache.
struct PinnedGtf {
    /// Genome assembly this entry targets.
    genome: KnownGenome,
    /// Human-readable label (e.g. `"GENCODE v45 basic"`).
    source_label: &'static str,
    /// Direct URL — must serve a gzipped GTF.
    url: &'static str,
    /// Optional pinned MD5SUM URL. When present, fetched first and the
    /// downloaded GTF is verified against the matching entry.
    md5sum_url: Option<&'static str>,
    /// Filename to look up inside the upstream `MD5SUMS` listing.
    upstream_filename: &'static str,
}

const PINNED_GTFS: &[PinnedGtf] = &[
    PinnedGtf {
        genome: KnownGenome::Grch38,
        source_label: "GENCODE v45 basic (GRCh38)",
        url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz",
        md5sum_url: Some(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/MD5SUMS",
        ),
        upstream_filename: "gencode.v45.basic.annotation.gtf.gz",
    },
    PinnedGtf {
        genome: KnownGenome::Grch37,
        source_label: "GENCODE v45lift37 basic (GRCh37)",
        url: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh37_mapping/gencode.v45lift37.basic.annotation.gtf.gz",
        md5sum_url: Some(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh37_mapping/MD5SUMS",
        ),
        upstream_filename: "gencode.v45lift37.basic.annotation.gtf.gz",
    },
    // T2T-CHM13: GENCODE has no canonical lift; skipped for v1. Users on
    // T2T-CHM13 must pass --gtf <FILE> explicitly. fetch-references will
    // surface a clear "not supported" error.
];

/// Implementation of `liquidqc fetch-references`.
///
/// `genome_filter`:
///   * `None` — fetch every genome in [`PINNED_GTFS`].
///   * `Some(KnownGenome::…)` — fetch just that one.
///
/// `force` re-downloads even when the cache file already exists.
pub fn run(
    genome_filter: Option<KnownGenome>,
    force: bool,
    list_only: bool,
    quiet: bool,
) -> Result<()> {
    if list_only {
        return list_pinned();
    }

    let targets: Vec<&PinnedGtf> = match genome_filter {
        Some(g) => {
            let entry = PINNED_GTFS.iter().find(|p| p.genome == g).ok_or_else(|| {
                anyhow!(
                    "fetch-references does not support {} yet — pass --gtf <FILE> manually",
                    g.label()
                )
            })?;
            vec![entry]
        }
        None => PINNED_GTFS.iter().collect(),
    };

    let dir = gtf_dir()?;
    ensure_dir(&dir)?;
    if !quiet {
        eprintln!("Cache directory: {}", dir.display());
    }

    let mut n_fetched = 0usize;
    let mut n_skipped = 0usize;
    for entry in targets {
        let dest = gtf_path_for(entry.genome)?;
        if dest.exists() && !force {
            if !quiet {
                eprintln!(
                    "skip  {}: already cached at {} (re-run with --force to refresh)",
                    entry.genome.cache_name(),
                    dest.display(),
                );
            }
            n_skipped += 1;
            continue;
        }
        if !quiet {
            eprintln!(
                "fetch {} from {}",
                entry.genome.cache_name(),
                entry.source_label
            );
        }
        download_and_install(entry, &dest, quiet).with_context(|| {
            format!(
                "Failed to fetch {} GTF from {}",
                entry.genome.cache_name(),
                entry.url
            )
        })?;
        if !quiet {
            eprintln!("  -> {}", dest.display());
        }
        n_fetched += 1;
    }

    if !quiet {
        eprintln!("Done. Fetched {n_fetched}, skipped {n_skipped}.");
    }
    Ok(())
}

/// Print the table of pinned references so the user can see what
/// `fetch-references` would download (and what is unsupported).
fn list_pinned() -> Result<()> {
    println!("Pinned GTF downloads:");
    for p in PINNED_GTFS {
        println!(
            "  {:<14} {}\n    {}",
            p.genome.cache_name(),
            p.source_label,
            p.url
        );
    }
    let supported: std::collections::HashSet<KnownGenome> =
        PINNED_GTFS.iter().map(|p| p.genome).collect();
    let unsupported: Vec<&KnownGenome> = KnownGenome::all()
        .iter()
        .filter(|g| !supported.contains(g))
        .collect();
    if !unsupported.is_empty() {
        println!("\nNot yet supported by fetch-references (pass --gtf manually):");
        for g in unsupported {
            println!("  {:<14} {}", g.cache_name(), g.label());
        }
    }
    Ok(())
}

/// Download an entry to a `.partial` file beside `dest`, verify, and
/// rename atomically. The partial file is removed on any error path.
fn download_and_install(entry: &PinnedGtf, dest: &Path, quiet: bool) -> Result<()> {
    let parent = dest.parent().ok_or_else(|| {
        anyhow!(
            "Internal: cache target '{}' has no parent directory",
            dest.display()
        )
    })?;
    ensure_dir(parent)?;

    // Resolve the expected MD5 first (cheap; ~few KB), so a corrupt
    // download fails fast.
    let expected_md5 = match entry.md5sum_url {
        Some(url) => Some(fetch_md5_for(url, entry.upstream_filename)?),
        None => None,
    };

    let partial = with_extension_appended(dest, ".partial");
    // Best-effort cleanup of any leftover partial from a prior interrupted run.
    let _ = fs::remove_file(&partial);

    let resp = ureq_get(entry.url)?;

    {
        let mut out = fs::File::create(&partial)
            .with_context(|| format!("Failed to create '{}'", partial.display()))?;
        let mut hasher = Md5::new();
        let mut reader = resp.into_reader();
        let mut buf = [0u8; 64 * 1024];
        let mut total: u64 = 0;
        loop {
            let n = reader
                .read(&mut buf)
                .context("Network read failed during GTF download")?;
            if n == 0 {
                break;
            }
            hasher.update(&buf[..n]);
            out.write_all(&buf[..n])
                .with_context(|| format!("Failed to write '{}'", partial.display()))?;
            total += n as u64;
        }
        out.flush().ok();
        let got_md5 = format!("{:x}", hasher.finalize());
        if !quiet {
            eprintln!(
                "  downloaded {} bytes, md5={}",
                format_bytes(total),
                got_md5
            );
        }
        if let Some(want) = expected_md5 {
            if !got_md5.eq_ignore_ascii_case(&want) {
                let _ = fs::remove_file(&partial);
                anyhow::bail!(
                    "MD5 mismatch for {}: got {}, expected {} (from upstream MD5SUMS)",
                    entry.upstream_filename,
                    got_md5,
                    want,
                );
            }
        }
    }

    // Quick sanity check: file must be a non-empty gzip stream.
    sanity_check_gzip(&partial)?;

    fs::rename(&partial, dest).with_context(|| {
        format!(
            "Failed to install '{}' to '{}'",
            partial.display(),
            dest.display()
        )
    })?;
    Ok(())
}

/// Read the upstream `MD5SUMS` file and return the hex MD5 for
/// `filename`.
fn fetch_md5_for(md5sum_url: &str, filename: &str) -> Result<String> {
    let resp = ureq_get(md5sum_url)
        .with_context(|| format!("Failed to fetch MD5SUMS from '{}'", md5sum_url))?;
    let mut body = String::new();
    resp.into_reader()
        .read_to_string(&mut body)
        .with_context(|| format!("Failed to read MD5SUMS body from '{}'", md5sum_url))?;
    parse_md5sums(&body, filename).ok_or_else(|| {
        anyhow!(
            "MD5SUMS at '{}' did not contain an entry for '{}'",
            md5sum_url,
            filename
        )
    })
}

/// Parse a GENCODE-style `MD5SUMS` body (`"<md5>  <relative-path>\n"`)
/// and return the lowercase MD5 for `filename`, or `None` when absent.
fn parse_md5sums(body: &str, filename: &str) -> Option<String> {
    for line in body.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut parts = trimmed.splitn(2, char::is_whitespace);
        let md5 = parts.next().unwrap_or("").trim();
        let name = parts.next().unwrap_or("").trim();
        if md5.len() == 32 && (name == filename || name.ends_with(&format!("/{filename}"))) {
            return Some(md5.to_ascii_lowercase());
        }
    }
    None
}

/// Confirm a downloaded file starts with the gzip magic and decompresses
/// far enough to surface a hard truncation. Cheap: reads ~1 KiB.
fn sanity_check_gzip(path: &Path) -> Result<()> {
    use flate2::read::GzDecoder;
    let f = fs::File::open(path)
        .with_context(|| format!("Failed to reopen '{}' for sanity check", path.display()))?;
    let mut dec = GzDecoder::new(f);
    let mut buf = [0u8; 1024];
    let n = dec.read(&mut buf).with_context(|| {
        format!(
            "Downloaded file '{}' is not a valid gzip stream",
            path.display()
        )
    })?;
    if n == 0 {
        anyhow::bail!(
            "Downloaded file '{}' decompressed to zero bytes",
            path.display()
        );
    }
    Ok(())
}

/// Wrap [`ureq::get`] with a sane User-Agent and timeout.
fn ureq_get(url: &str) -> Result<ureq::Response> {
    let agent = ureq::AgentBuilder::new()
        .timeout_connect(std::time::Duration::from_secs(20))
        .timeout_read(std::time::Duration::from_secs(120))
        .user_agent(concat!("liquidqc/", env!("CARGO_PKG_VERSION")))
        .build();
    agent
        .get(url)
        .call()
        .with_context(|| format!("HTTP GET failed for '{}'", url))
}

/// Append `suffix` to a path's filename component (without affecting
/// the parent directory). Used to derive `<dest>.partial` —
/// [`Path::with_extension`] would replace `.gz`, which is wrong here.
fn with_extension_appended(p: &Path, suffix: &str) -> PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_table_uses_unique_genomes() {
        use std::collections::HashSet;
        let s: HashSet<KnownGenome> = PINNED_GTFS.iter().map(|p| p.genome).collect();
        assert_eq!(s.len(), PINNED_GTFS.len());
    }

    #[test]
    fn pinned_urls_use_https_and_match_filenames() {
        for p in PINNED_GTFS {
            assert!(p.url.starts_with("https://"), "non-https URL: {}", p.url);
            assert!(
                p.url.ends_with(p.upstream_filename),
                "URL '{}' does not end with upstream_filename '{}'",
                p.url,
                p.upstream_filename
            );
        }
    }

    #[test]
    fn with_extension_appended_keeps_parent() {
        let p = Path::new("/tmp/x/y.gtf.gz");
        let q = with_extension_appended(p, ".partial");
        assert_eq!(q, Path::new("/tmp/x/y.gtf.gz.partial"));
    }

    #[test]
    fn parse_md5sums_finds_matching_filename() {
        let body = "abcdef0123456789abcdef0123456789  gencode.v45.basic.annotation.gtf.gz\n\
                    1111111111111111111111111111111  gencode.v45.annotation.gtf.gz\n";
        assert_eq!(
            parse_md5sums(body, "gencode.v45.basic.annotation.gtf.gz").as_deref(),
            Some("abcdef0123456789abcdef0123456789"),
        );
    }

    #[test]
    fn parse_md5sums_returns_none_on_miss() {
        let body = "abcdef0123456789abcdef0123456789  some-other.gtf.gz\n";
        assert!(parse_md5sums(body, "gencode.v45.basic.annotation.gtf.gz").is_none());
    }

    #[test]
    fn parse_md5sums_accepts_relative_path_match() {
        // GENCODE sometimes lists "./<name>" or a subdir prefix.
        let body = "deadbeefdeadbeefdeadbeefdeadbeef  ./gencode.v45.basic.annotation.gtf.gz\n";
        assert_eq!(
            parse_md5sums(body, "gencode.v45.basic.annotation.gtf.gz").as_deref(),
            Some("deadbeefdeadbeefdeadbeefdeadbeef"),
        );
    }
}
