//! Per-user cache for reference assets (currently: GTF gene annotations).
//!
//! Resolution order for the cache root:
//!   1. `LIQUIDQC_CACHE_DIR` — explicit override (per-run or test).
//!   2. `XDG_CACHE_HOME/liquidqc` if `XDG_CACHE_HOME` is set and absolute.
//!   3. `dirs::cache_dir()/liquidqc` — platform default
//!      (`~/.cache/liquidqc` on Linux, `~/Library/Caches/liquidqc` on macOS).
//!
//! Population is handled by the `fetch-references` subcommand; this module
//! only provides lookup helpers and path resolution. Lookup is read-only —
//! `liquidqc rna` never writes to the cache.

use crate::genome::KnownGenome;
use crate::io::GZIP_MAGIC;
use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::{Path, PathBuf};

/// Subdirectory under the cache root holding GTF files.
const GTF_SUBDIR: &str = "gtf";

/// Resolve the cache root. Returns an error only if no candidate is set
/// (e.g. on a stripped-down container with no `HOME`).
pub fn cache_root() -> Result<PathBuf> {
    if let Ok(v) = std::env::var("LIQUIDQC_CACHE_DIR") {
        if !v.is_empty() {
            return Ok(PathBuf::from(v));
        }
    }
    if let Ok(v) = std::env::var("XDG_CACHE_HOME") {
        let p = PathBuf::from(&v);
        if p.is_absolute() {
            return Ok(p.join("liquidqc"));
        }
    }
    dirs::cache_dir()
        .map(|p| p.join("liquidqc"))
        .ok_or_else(|| {
            anyhow!(
                "Could not determine a cache directory. Set LIQUIDQC_CACHE_DIR or HOME, \
                 or pass --gtf <FILE> explicitly."
            )
        })
}

/// Directory holding cached GTF files.
pub fn gtf_dir() -> Result<PathBuf> {
    Ok(cache_root()?.join(GTF_SUBDIR))
}

/// Path at which a given genome's GTF lives once cached.
///
/// Files are named `<genome>.gtf.gz` (always gzip — the rest of the
/// pipeline auto-detects compression from magic bytes via [`crate::io`]).
pub fn gtf_path_for(genome: KnownGenome) -> Result<PathBuf> {
    Ok(gtf_dir()?.join(format!("{}.gtf.gz", genome.cache_name())))
}

/// Status of a cache lookup. `Stale` means a file exists at the cache path
/// but failed validation (e.g. truncated download from a pre-`sanity_check_gzip`
/// liquidqc release) and the caller should warn + treat the cache as missing.
#[derive(Debug, PartialEq, Eq)]
pub enum CacheLookup {
    Hit(PathBuf),
    Missing,
    Stale(PathBuf),
}

/// Look up a cached GTF for the given genome.
///
/// Returns `Hit` when the file exists and looks like a valid gzip stream
/// (first two bytes are `1f 8b`). A file that exists but fails the magic-byte
/// check is reported as `Stale` so the caller can prompt the user to re-fetch.
/// IO errors short-circuit to `Err`.
pub fn lookup_gtf(genome: KnownGenome) -> Result<CacheLookup> {
    let p = gtf_path_for(genome)?;
    match has_gzip_magic(&p) {
        Ok(true) => Ok(CacheLookup::Hit(p)),
        Ok(false) => Ok(CacheLookup::Stale(p)),
        Err(e) if e.kind() == ErrorKind::NotFound => Ok(CacheLookup::Missing),
        Err(e) => {
            Err(anyhow::Error::new(e).context(format!("opening cache file: {}", p.display())))
        }
    }
}

/// Read the first two bytes of `path` and return whether they match the
/// gzip magic number. A file shorter than 2 bytes is reported as not-gzip.
/// Surfaces the raw [`std::io::Error`] so callers can distinguish missing
/// files from other I/O failures without a separate `is_file` check.
fn has_gzip_magic(path: &Path) -> std::io::Result<bool> {
    let mut f = File::open(path)?;
    let mut buf = [0u8; 2];
    match f.read_exact(&mut buf) {
        Ok(()) => Ok(buf == GZIP_MAGIC),
        Err(e) if e.kind() == ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    }
}

/// Build a clear error message instructing the user to run
/// `liquidqc fetch-references`. Used when a BAM fingerprints to a known
/// genome but the cache is empty.
pub fn missing_gtf_error(genome: KnownGenome) -> String {
    let cache_path = gtf_path_for(genome)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| format!("<cache>/{}/{}.gtf.gz", GTF_SUBDIR, genome.cache_name()));
    format!(
        "BAM fingerprints to {} but no cached GTF was found at:\n  {}\n\n\
         Run `liquidqc fetch-references --genome {}` to populate the cache, \
         or pass `--gtf <FILE>` explicitly.",
        genome.label(),
        cache_path,
        genome.cache_name(),
    )
}

/// Ensure that `dir` exists, creating it (and parents) if not.
pub fn ensure_dir(dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir).map_err(|e| {
        anyhow!(
            "Failed to create cache directory '{}': {}",
            dir.display(),
            e
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::sync::{Mutex, MutexGuard, OnceLock};

    /// Cache-tests share a process-wide env var, so they must run serially.
    /// Each test acquires this lock and holds it through its env mutations.
    fn env_lock() -> MutexGuard<'static, ()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        // Poisoned lock from a panicked test still serialises subsequent
        // ones — recover the inner guard rather than re-panicking.
        LOCK.get_or_init(|| Mutex::new(()))
            .lock()
            .unwrap_or_else(|e| e.into_inner())
    }

    /// Set `LIQUIDQC_CACHE_DIR` for the duration of a test, restoring the
    /// previous value when dropped. Holds [`env_lock`] for its lifetime
    /// so concurrent tests don't see partial state.
    struct EnvGuard {
        key: &'static str,
        prev: Option<String>,
        _g: MutexGuard<'static, ()>,
    }
    impl EnvGuard {
        fn set(key: &'static str, val: &str) -> Self {
            let g = env_lock();
            let prev = std::env::var(key).ok();
            // SAFETY: env_lock serialises all tests in this module that
            // mutate the process environment.
            unsafe { std::env::set_var(key, val) };
            Self { key, prev, _g: g }
        }
    }
    impl Drop for EnvGuard {
        fn drop(&mut self) {
            unsafe {
                match &self.prev {
                    Some(v) => std::env::set_var(self.key, v),
                    None => std::env::remove_var(self.key),
                }
            }
        }
    }

    #[test]
    fn override_env_takes_precedence() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let root = cache_root().unwrap();
        assert_eq!(root, tmp.path());
        let gtf = gtf_path_for(KnownGenome::Grch38).unwrap();
        assert_eq!(gtf, tmp.path().join("gtf").join("GRCh38.gtf.gz"));
    }

    #[test]
    fn lookup_returns_missing_when_absent() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        assert_eq!(
            lookup_gtf(KnownGenome::Grch38).unwrap(),
            CacheLookup::Missing
        );
    }

    #[test]
    fn lookup_finds_present_file() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let p = gtf_path_for(KnownGenome::Grch38).unwrap();
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        // Write a 2-byte gzip header so the magic-byte check passes.
        std::fs::write(&p, GZIP_MAGIC).unwrap();
        let got = lookup_gtf(KnownGenome::Grch38).unwrap();
        assert_eq!(got, CacheLookup::Hit(p));
    }

    #[test]
    fn lookup_reports_stale_for_non_gzip() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let p = gtf_path_for(KnownGenome::Grch38).unwrap();
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        // 11 bytes of plain text — the exact failure mode observed in practice.
        std::fs::write(&p, b"not-a-gzip!").unwrap();
        let got = lookup_gtf(KnownGenome::Grch38).unwrap();
        assert_eq!(got, CacheLookup::Stale(p));
    }

    #[test]
    fn lookup_reports_stale_for_truncated_below_magic() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let p = gtf_path_for(KnownGenome::Grch38).unwrap();
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        // 1 byte: shorter than the 2-byte magic — must report Stale, not panic.
        std::fs::write(&p, [0x1fu8]).unwrap();
        let got = lookup_gtf(KnownGenome::Grch38).unwrap();
        assert_eq!(got, CacheLookup::Stale(p));
    }

    #[test]
    fn missing_gtf_error_mentions_subcommand_and_path() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let msg = missing_gtf_error(KnownGenome::Grch38);
        assert!(msg.contains("fetch-references"));
        assert!(msg.contains("GRCh38"));
        assert!(msg.contains("--genome"));
    }
}
