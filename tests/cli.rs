//! CLI integration tests.
//!
//! These tests require the `cli` feature — run with:
//!   cargo test --features cli --test cli

use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn poa_bin() -> Command {
    let mut bin = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    bin.push("target/debug/poa-consensus");
    Command::new(bin)
}

fn build_cli() {
    let status = Command::new("cargo")
        .args(["build", "--features", "cli"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .expect("cargo build failed");
    assert!(status.success(), "CLI build failed");
}

// ── simple FASTA — single consensus ──────────────────────────────────────────

#[test]
fn cli_simple_fasta_consensus() {
    build_cli();
    let out = poa_bin()
        .arg(fixture("simple.fa"))
        .output()
        .expect("poa-consensus failed to run");

    assert!(out.status.success(), "exit code: {}", out.status);
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Header line
    assert!(
        lines[0].starts_with(">consensus"),
        "unexpected header: {}",
        lines[0]
    );
    assert!(
        lines[0].contains("reads=10"),
        "expected reads=10 in header: {}",
        lines[0]
    );

    // Sequence — majority wins; should be the ground truth
    let seq = lines[1];
    assert_eq!(
        seq, "ACGTACGTACGTACGTACGTACGTACGT",
        "consensus mismatch: {seq}"
    );
}

#[test]
fn cli_simple_fasta_correct_length() {
    build_cli();
    let out = poa_bin()
        .arg(fixture("simple.fa"))
        .output()
        .expect("poa-consensus failed to run");
    let stdout = String::from_utf8_lossy(&out.stdout);
    let seq = stdout.lines().nth(1).unwrap_or("");
    assert_eq!(seq.len(), 28, "wrong consensus length: {}", seq.len());
}

// ── noisy FASTA — SNPs don't corrupt consensus ────────────────────────────────

#[test]
fn cli_noisy_fasta_recovers_majority() {
    build_cli();
    let out = poa_bin()
        .arg(fixture("noisy.fa"))
        .output()
        .expect("poa-consensus failed to run");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let seq = stdout.lines().nth(1).unwrap_or("");
    assert_eq!(
        seq, "GCTAGCTAGCTAGCTA",
        "expected majority consensus: {seq}"
    );
}

// ── partial reads — semi-global mode ─────────────────────────────────────────

#[test]
fn cli_partial_reads_semi_global() {
    build_cli();
    let out = poa_bin()
        .args(["--semi-global", fixture("partial.fa").to_str().unwrap()])
        .output()
        .expect("poa-consensus failed to run");
    assert!(out.status.success(), "exit: {}", out.status);
    let stdout = String::from_utf8_lossy(&out.stdout);
    let seq = stdout.lines().nth(1).unwrap_or("");
    // Consensus must be the full-read ground truth; partial reads should not truncate it.
    assert_eq!(seq, "ATCGATCGATCGATCG", "semi-global consensus: {seq}");
}

// ── two-allele — multi mode ───────────────────────────────────────────────────

#[test]
fn cli_two_allele_multi_mode() {
    build_cli();
    let out = poa_bin()
        .args(["--multi", fixture("two_allele.fa").to_str().unwrap()])
        .output()
        .expect("poa-consensus failed to run");
    assert!(out.status.success(), "exit: {}", out.status);
    let stdout = String::from_utf8_lossy(&out.stdout);

    // Collect all FASTA records
    let mut seqs: Vec<&str> = Vec::new();
    for line in stdout.lines() {
        if !line.starts_with('>') {
            seqs.push(line);
        }
    }

    // Must detect exactly two alleles
    assert_eq!(seqs.len(), 2, "expected 2 alleles, got: {seqs:?}");

    // Both alleles are 22 bp; they differ at position 11 (A vs G).
    assert!(
        seqs.iter().all(|s| s.len() == 22),
        "allele lengths: {:?}",
        seqs
    );

    // One allele must have A and one must have G at position 11.
    let chars: Vec<char> = seqs.iter().map(|s| s.chars().nth(11).unwrap()).collect();
    assert!(
        chars.contains(&'A') && chars.contains(&'G'),
        "expected A+G split at pos 11, got: {chars:?}"
    );
}

// ── stdin input ───────────────────────────────────────────────────────────────

#[test]
fn cli_stdin_input() {
    build_cli();
    let fasta = std::fs::read(fixture("simple.fa")).expect("read fixture");
    let mut child = poa_bin()
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("spawn failed");

    use std::io::Write as _;
    child.stdin.take().unwrap().write_all(&fasta).unwrap();
    let out = child.wait_with_output().expect("wait failed");

    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let seq = stdout.lines().nth(1).unwrap_or("");
    assert_eq!(seq, "ACGTACGTACGTACGTACGTACGTACGT");
}

// ── single-record passthrough ─────────────────────────────────────────────────

#[test]
fn cli_single_record_passthrough() {
    build_cli();
    let fasta = b">solo\nACGTACGT\n";
    let mut child = poa_bin()
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("spawn failed");

    use std::io::Write as _;
    child.stdin.take().unwrap().write_all(fasta).unwrap();
    let out = child.wait_with_output().expect("wait failed");

    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let seq = stdout.lines().nth(1).unwrap_or("");
    assert_eq!(seq, "ACGTACGT");
}

// ── band flag propagated to header ───────────────────────────────────────────

#[test]
fn cli_band_header_fixed() {
    build_cli();
    let out = poa_bin()
        .args(["-b", "50", fixture("simple.fa").to_str().unwrap()])
        .output()
        .expect("poa-consensus failed to run");
    assert!(out.status.success());
    let header = String::from_utf8_lossy(&out.stdout)
        .lines()
        .next()
        .unwrap_or("")
        .to_string();
    assert!(header.contains("band=50"), "header: {header}");
}

#[test]
fn cli_band_header_adaptive() {
    build_cli();
    let out = poa_bin()
        .args(["--adaptive-band", fixture("simple.fa").to_str().unwrap()])
        .output()
        .expect("poa-consensus failed to run");
    assert!(out.status.success());
    let header = String::from_utf8_lossy(&out.stdout)
        .lines()
        .next()
        .unwrap_or("")
        .to_string();
    assert!(header.contains("band=adaptive"), "header: {header}");
}

// ── error handling ────────────────────────────────────────────────────────────

#[test]
fn cli_empty_input_exits_nonzero() {
    build_cli();
    let mut child = poa_bin()
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn");
    drop(child.stdin.take());
    let out = child.wait_with_output().expect("wait");
    assert!(
        !out.status.success(),
        "expected non-zero exit for empty input"
    );
}

#[test]
fn cli_bad_format_exits_nonzero() {
    build_cli();
    let bad = b"not fasta or fastq\n";
    let mut child = poa_bin()
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn");
    use std::io::Write as _;
    child.stdin.take().unwrap().write_all(bad).unwrap();
    let out = child.wait_with_output().expect("wait");
    assert!(!out.status.success());
}
