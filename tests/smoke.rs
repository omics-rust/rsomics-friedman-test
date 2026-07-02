//! End-to-end smoke tests against the built binary.

use std::io::Write;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rsomics-friedman-test")
}

fn write_temp(contents: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let dir = std::env::temp_dir();
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let path = dir.join(format!("friedman-smoke-{}-{id}.tsv", std::process::id()));
    let mut f = std::fs::File::create(&path).unwrap();
    f.write_all(contents.as_bytes()).unwrap();
    path
}

const MATRIX: &str = "1\t2\t3\n2\t1\t3\n3\t1\t2\n1\t3\t2\n";

#[test]
fn text_output_three_fields() {
    let path = write_temp(MATRIX);
    let out = Command::new(bin()).arg(&path).output().unwrap();
    std::fs::remove_file(&path).ok();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let parts: Vec<&str> = s.trim().split('\t').collect();
    assert_eq!(parts.len(), 3, "got {s:?}");
    assert_eq!(parts[1], "2");
}

#[test]
fn json_envelope() {
    let path = write_temp(MATRIX);
    let out = Command::new(bin())
        .arg(&path)
        .arg("--json")
        .output()
        .unwrap();
    std::fs::remove_file(&path).ok();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("\"status\":\"ok\""), "got {s:?}");
    assert!(s.contains("\"result\""), "got {s:?}");
    assert!(s.contains("\"Q\""), "got {s:?}");
    assert!(s.contains("\"df\""), "got {s:?}");
    assert!(s.contains("\"p\""), "got {s:?}");
    assert_eq!(
        s.trim().lines().count(),
        1,
        "json must be single line: {s:?}"
    );
}

#[test]
fn stdin_input() {
    let mut child = Command::new(bin())
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(MATRIX.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert_eq!(s.trim().split('\t').count(), 3);
}

#[test]
fn help_exits_zero() {
    let out = Command::new(bin()).arg("--help").output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("Friedman"), "got {s:?}");
}

#[test]
fn rejects_two_treatments() {
    let path = write_temp("1\t2\n2\t1\n");
    let out = Command::new(bin()).arg(&path).output().unwrap();
    std::fs::remove_file(&path).ok();
    assert!(!out.status.success());
}

/// Zero within-block variance used to hang forever in the survival-function
/// continued fraction. The binary must terminate and report NaN like SciPy.
#[test]
fn all_tied_terminates_with_nan() {
    let path = write_temp("5\t5\t5\n2\t2\t2\n9\t9\t9\n");
    let out = Command::new(bin()).arg(&path).output().unwrap();
    std::fs::remove_file(&path).ok();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let parts: Vec<&str> = s.trim().split('\t').collect();
    assert_eq!(parts.len(), 3, "got {s:?}");
    assert_eq!(parts[0], "NaN", "Q field should be NaN, got {s:?}");
    assert_eq!(parts[2], "NaN", "p field should be NaN, got {s:?}");
}
