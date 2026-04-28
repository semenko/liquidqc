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
use anyhow::{anyhow, Result};
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

/// Look up a cached GTF for the given genome. Returns the path if it
/// exists, or `None` if the file is missing. An IO error short-circuits
/// to `Err`.
pub fn lookup_gtf(genome: KnownGenome) -> Result<Option<PathBuf>> {
    let p = gtf_path_for(genome)?;
    if p.is_file() {
        Ok(Some(p))
    } else {
        Ok(None)
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
    fn lookup_returns_none_when_missing() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        assert!(lookup_gtf(KnownGenome::Grch38).unwrap().is_none());
    }

    #[test]
    fn lookup_finds_present_file() {
        let tmp = tempfile::tempdir().unwrap();
        let _g = EnvGuard::set("LIQUIDQC_CACHE_DIR", tmp.path().to_str().unwrap());
        let p = gtf_path_for(KnownGenome::Grch38).unwrap();
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(&p, b"placeholder").unwrap();
        let got = lookup_gtf(KnownGenome::Grch38).unwrap();
        assert_eq!(got, Some(p));
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
