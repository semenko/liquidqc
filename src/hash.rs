//! Streaming MD5 helpers for input-file fingerprinting.
//!
//! The schema v1 envelope requires `bam_md5`, `gtf_md5`, and (optionally)
//! `reference_fasta_md5`. These helpers compute MD5 in 1 MiB chunks so we
//! can hash multi-GB BAMs without buffering the whole file.
//!
//! Per the v1 contract:
//! - BAM hashes the file bytes as on disk (BGZF). [`md5_and_size`]
//! - GTF hashes the **uncompressed** bytes (gzip-transparent via [`crate::io::open_reader`]). [`md5_uncompressed`]
//! - FASTA hashes the file bytes as supplied. [`md5`]

use anyhow::{Context, Result};
use md5::{Digest, Md5};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Streaming chunk size for hashing large files.
const HASH_CHUNK_BYTES: usize = 1 << 20; // 1 MiB

/// Compute MD5 (lowercase hex) and byte size of a file in a single streaming pass.
pub fn md5_and_size<P: AsRef<Path>>(path: P) -> Result<(String, u64)> {
    let path = path.as_ref();
    let f =
        File::open(path).with_context(|| format!("opening file for md5: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(HASH_CHUNK_BYTES, f);
    let mut hasher = Md5::new();
    let mut buf = vec![0u8; HASH_CHUNK_BYTES];
    let mut total: u64 = 0;
    loop {
        let n = reader
            .read(&mut buf)
            .with_context(|| format!("reading file for md5: {}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
        total += n as u64;
    }
    Ok((format!("{:x}", hasher.finalize()), total))
}

/// Compute MD5 (lowercase hex) of a file's bytes as on disk.
pub fn md5<P: AsRef<Path>>(path: P) -> Result<String> {
    Ok(md5_and_size(path)?.0)
}

/// Compute MD5 (lowercase hex) of the **decompressed** contents of a file.
///
/// Uses [`crate::io::open_reader`] so the path is gzip-transparent: the hash
/// is over the uncompressed bytes regardless of whether the file on disk is
/// `.gtf` or `.gtf.gz`.
pub fn md5_uncompressed<P: AsRef<Path>>(path: P) -> Result<String> {
    let path = path.as_ref();
    let mut reader = crate::io::open_reader(path)
        .with_context(|| format!("opening for uncompressed md5: {}", path.display()))?;
    let mut hasher = Md5::new();
    let mut buf = vec![0u8; HASH_CHUNK_BYTES];
    loop {
        let n = reader
            .read(&mut buf)
            .with_context(|| format!("reading for uncompressed md5: {}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn write_temp(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    fn write_temp_gz(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        let f = File::create(&path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap();
        path
    }

    /// MD5("hello world") = 5eb63bbbe01eeed093cb22bb8f5acdc3 (well-known reference).
    const HELLO_WORLD_MD5: &str = "5eb63bbbe01eeed093cb22bb8f5acdc3";

    #[test]
    fn md5_and_size_known_vector() {
        let p = write_temp("liquidqc_hash_known.bin", b"hello world");
        let (h, n) = md5_and_size(&p).unwrap();
        assert_eq!(h, HELLO_WORLD_MD5);
        assert_eq!(n, 11);
    }

    #[test]
    fn md5_returns_lowercase_32_hex() {
        let p = write_temp("liquidqc_hash_lower.bin", b"hello world");
        let h = md5(&p).unwrap();
        assert_eq!(h.len(), 32);
        assert!(h
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()));
    }

    #[test]
    fn md5_uncompressed_matches_plain() {
        let plain = write_temp("liquidqc_hash_plain.txt", b"hello world");
        let gz = write_temp_gz("liquidqc_hash_plain.txt.gz", b"hello world");
        let h_plain = md5_uncompressed(&plain).unwrap();
        let h_gz = md5_uncompressed(&gz).unwrap();
        assert_eq!(h_plain, HELLO_WORLD_MD5);
        assert_eq!(h_gz, HELLO_WORLD_MD5);
    }

    #[test]
    fn md5_chunk_boundary() {
        // File larger than one HASH_CHUNK_BYTES window exercises the streaming loop.
        let bytes = vec![0xABu8; HASH_CHUNK_BYTES + 17];
        let p = write_temp("liquidqc_hash_big.bin", &bytes);
        let (h, n) = md5_and_size(&p).unwrap();
        assert_eq!(n, bytes.len() as u64);
        let mut hasher = Md5::new();
        hasher.update(&bytes);
        assert_eq!(h, format!("{:x}", hasher.finalize()));
    }
}
