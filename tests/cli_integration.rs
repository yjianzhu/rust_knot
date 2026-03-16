use std::fs;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::Command;

use rust_knot::io::read_data_xyz_frames;
use tempfile::tempdir;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

fn run_cli(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rust_knot"))
        .args(args)
        .output()
        .expect("failed to run rust_knot binary")
}

#[test]
fn help_includes_threads_and_no_command_mode() {
    let output = run_cli(&["--help"]);
    assert!(output.status.success());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--threads <n>"),
        "help missing --threads option:\n{stderr}"
    );
    assert!(
        !stderr.contains("<command>"),
        "help should not advertise command mode:\n{stderr}"
    );
}

#[test]
fn threads_zero_is_rejected() {
    let input = fixture_path("L400_knot4_1_open.xyz");
    let input_s = input.to_str().expect("non-UTF8 fixture path");

    let output = run_cli(&[input_s, "--threads", "0"]);
    assert!(
        !output.status.success(),
        "command should fail for --threads 0"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--threads must be at least 1"),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn unknown_flag_is_rejected() {
    let input = fixture_path("L400_knot4_1_open.xyz");
    let input_s = input.to_str().expect("non-UTF8 fixture path");

    let output = run_cli(&[input_s, "--threds", "2"]);
    assert!(
        !output.status.success(),
        "command should fail for unknown flag"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unknown flag '--threds'"),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn multi_frame_cli_respects_batch_and_threads() {
    let frame_path = fixture_path("L400_knot4_1_open.xyz");
    let frame = fs::read_to_string(&frame_path).expect("failed to read fixture");
    let expected_frames_per_copy = {
        let file = fs::File::open(&frame_path).expect("failed to reopen fixture");
        let mut reader = BufReader::new(file);
        read_data_xyz_frames(&mut reader)
            .expect("failed to parse fixture")
            .len()
    };
    let expected_total_frames = expected_frames_per_copy * 2;

    let temp = tempdir().expect("failed to create temp dir");
    let multi_frame_path = temp.path().join("multi_frame.xyz");
    let output_path = temp.path().join("knot_index.txt");
    fs::write(&multi_frame_path, format!("{frame}{frame}"))
        .expect("failed to write multi-frame xyz");

    let input_s = multi_frame_path.to_str().expect("non-UTF8 temp input path");
    let output_s = output_path.to_str().expect("non-UTF8 temp output path");
    let output = run_cli(&[
        input_s,
        "--batch",
        "1",
        "--threads",
        "2",
        "--output",
        output_s,
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let report = fs::read_to_string(&output_path).expect("failed to read output report");
    let lines: Vec<&str> = report.lines().collect();
    assert_eq!(
        lines.first().copied(),
        Some("# frame\tknottype\tknot_start\tknot_end\tknot_size"),
        "unexpected report header:\n{report}"
    );
    assert_eq!(
        lines.len(),
        expected_total_frames + 1,
        "expected header + {expected_total_frames} frame lines:\n{report}"
    );
    assert!(
        lines[1].starts_with("0\t4_1\t"),
        "unexpected first frame line:\n{report}"
    );
    assert!(
        lines[expected_total_frames].starts_with(&format!("{}\t4_1\t", expected_total_frames - 1)),
        "unexpected last frame line:\n{report}"
    );
}

#[test]
fn unknown_knot_writes_polynomial_to_knot_index() {
    let input = fixture_path("knot_10_100.xyz");
    let temp = tempdir().expect("failed to create temp dir");
    let output_path = temp.path().join("unknown_poly.txt");
    let input_s = input.to_str().expect("non-UTF8 fixture path");
    let output_s = output_path.to_str().expect("non-UTF8 output path");

    let output = run_cli(&[input_s, "--output", output_s]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let expected_poly = "1-4*t+9*t^2-12*t^3+13*t^4-12*t^5+9*t^6-4*t^7+t^8";
    let report = fs::read_to_string(&output_path).expect("failed to read output report");
    let lines: Vec<&str> = report.lines().collect();
    assert_eq!(
        lines.len(),
        2,
        "expected header + one frame line:\n{report}"
    );
    assert_eq!(
        lines[0],
        "# frame\tknottype\tknot_start\tknot_end\tknot_size"
    );
    assert_eq!(lines[1], format!("0\t{expected_poly}\t-1\t-1\t0"));

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("knot polynomial not found in table"),
        "expected notfound warning in stderr:\n{stderr}"
    );
}
