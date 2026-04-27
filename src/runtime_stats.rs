//! Process-level runtime statistics for the v1 envelope.
//!
//! Currently exposes peak resident-set-size via `getrusage(RUSAGE_SELF)`.
//! Linux reports `ru_maxrss` in KiB; macOS reports it in bytes — this module
//! normalizes both to MiB. On non-unix targets the helper returns 0.0.

/// Peak resident set size of the current process, in MiB.
///
/// `getrusage` reports the process-level high-water mark, not per-thread or
/// per-BAM. In multi-BAM runs every per-sample envelope therefore records
/// the same final value; this is honest but worth keeping in mind.
#[cfg(unix)]
pub fn peak_rss_mb() -> f64 {
    use std::mem::MaybeUninit;
    // SAFETY: `getrusage` is a thread-safe POSIX syscall; we pass a valid
    // pointer to a stack-allocated `rusage` struct it is allowed to fill.
    unsafe {
        let mut ru = MaybeUninit::<libc::rusage>::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, ru.as_mut_ptr()) != 0 {
            return 0.0;
        }
        let ru = ru.assume_init();
        let bytes = if cfg!(target_os = "macos") {
            // macOS reports bytes (XNU diverges from POSIX here).
            ru.ru_maxrss as f64
        } else {
            // Linux reports kibibytes.
            ru.ru_maxrss as f64 * 1024.0
        };
        bytes / (1024.0 * 1024.0)
    }
}

#[cfg(not(unix))]
pub fn peak_rss_mb() -> f64 {
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_rss_mb_is_non_negative() {
        let v = peak_rss_mb();
        assert!(v >= 0.0, "peak_rss_mb returned negative: {v}");
        // Sanity bound: no test process should be using terabytes.
        assert!(v < 1_000_000.0, "peak_rss_mb implausibly large: {v}");
    }
}
