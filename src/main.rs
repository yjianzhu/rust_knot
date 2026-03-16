use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::time::Instant;

use rayon::ThreadPoolBuilder;
use rust_knot::alexander_table::AlexanderTable;
use rust_knot::batch::process_frames_streaming;
use rust_knot::config::KnotConfig;

fn print_usage(prog: &str) {
    eprintln!(
        "Usage: {prog} <xyz_file> [--table <path>] [target_type] [--ring] [--no-fast] [--debug] [--output <path>] [--batch <size>] [--threads <n>]"
    );
    eprintln!();
    eprintln!("  xyz_file:        path to XYZ coordinate file (single or multi-frame)");
    eprintln!("  --table <path>:  Alexander polynomial table (default: built-in ≤9 crossings)");
    eprintln!("  target_type:     (optional) knot type to search for core, e.g. '3_1'");
    eprintln!("  --ring:          treat chain as a closed ring");
    eprintln!("  --no-fast:       disable KMT simplification (enabled by default)");
    eprintln!("  --debug:         enable debug output");
    eprintln!("  --output <path>: write knot_index log (default: knot_index.txt)");
    eprintln!("  --batch <size>:  frames per batch (default: 64)");
    eprintln!("  --threads <n>:   rayon worker thread count (default: auto)");
    eprintln!("  -h, --help:      show this message");
}

fn require_arg(args: &[String], i: usize, flag: &str) -> String {
    if i >= args.len() {
        eprintln!("error: {flag} requires a value");
        std::process::exit(1);
    }
    args[i].clone()
}

fn run_cli(args: &[String], prog: &str) {
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage(prog);
        std::process::exit(if args.is_empty() { 1 } else { 0 });
    }

    let xyz_path = &args[0];

    let mut config = KnotConfig {
        faster: true,
        ..KnotConfig::default()
    };
    let mut target_type: Option<String> = None;
    let mut output_path = String::from("knot_index.txt");
    let mut batch_size: Option<usize> = None;
    let mut num_threads: Option<usize> = None;
    let mut table_path: Option<String> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--ring" => config.is_ring = true,
            "--no-fast" => config.faster = false,
            "--debug" => config.debug = true,
            "--table" => {
                i += 1;
                table_path = Some(require_arg(args, i, "--table"));
            }
            "--output" => {
                i += 1;
                output_path = require_arg(args, i, "--output");
            }
            "--batch" => {
                i += 1;
                let val = require_arg(args, i, "--batch");
                batch_size = Some(val.parse().unwrap_or_else(|_| {
                    eprintln!("error: --batch value '{val}' is not a valid integer");
                    std::process::exit(1);
                }));
            }
            "--threads" => {
                i += 1;
                let val = require_arg(args, i, "--threads");
                num_threads = Some(val.parse().unwrap_or_else(|_| {
                    eprintln!("error: --threads value '{val}' is not a valid integer");
                    std::process::exit(1);
                }));
            }
            s if !s.starts_with('-') => target_type = Some(s.to_string()),
            other => {
                eprintln!("error: unknown flag '{other}'");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if batch_size == Some(0) {
        eprintln!("error: --batch must be at least 1");
        std::process::exit(1);
    }

    if let Some(n) = num_threads {
        if n == 0 {
            eprintln!("error: --threads must be at least 1");
            std::process::exit(1);
        }
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap_or_else(|e| {
                eprintln!("error: failed to configure rayon thread pool: {e}");
                std::process::exit(1);
            });
    }

    let t0 = Instant::now();
    let table = match table_path {
        Some(ref path) => {
            eprintln!("Loading table from {path} (merged with built-in)...");
            AlexanderTable::builtin_with_file(path).expect("failed to load Alexander table")
        }
        None => {
            eprintln!("Using built-in table (≤9 crossings)");
            AlexanderTable::builtin()
        }
    };
    eprintln!(
        "Table ready: {} polynomial entries in {:?}",
        table.len(),
        t0.elapsed()
    );

    let file = File::open(xyz_path).expect("failed to open XYZ file");
    let reader = BufReader::new(file);

    eprintln!(
        "Config: is_ring={}, faster={}, extend_factor={}, num_rotations={}, batch_size={}, threads={}",
        config.is_ring,
        config.faster,
        config.extend_factor,
        config.num_rotations,
        batch_size.unwrap_or(64),
        num_threads
            .map(|n| n.to_string())
            .unwrap_or_else(|| "auto".to_string())
    );

    let out_file = File::create(&output_path).expect("failed to create output file");
    let mut writer = BufWriter::new(out_file);
    writeln!(writer, "# frame\tknottype\tknot_start\tknot_end\tknot_size")
        .expect("write header failed");

    let t0 = Instant::now();
    let mut n_knotted = 0usize;
    let mut n_errors = 0usize;
    let mut first_frame = None;

    let total_frames = process_frames_streaming(
        reader,
        &table,
        &config,
        target_type.as_deref(),
        batch_size,
        |batch_results| {
            for r in batch_results {
                if let Some(ref err) = r.error {
                    eprintln!("frame {}: error: {err}", r.frame);
                    n_errors += 1;
                } else if r.knot_type != "1" {
                    n_knotted += 1;
                }

                for w in &r.warnings {
                    eprintln!("frame {}: warning: {w}", r.frame);
                }

                let ktype = if r.error.is_some() {
                    "ERROR"
                } else {
                    &r.knot_type
                };
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    r.frame, ktype, r.knot_start, r.knot_end, r.knot_size
                )
                .expect("write failed");
            }

            if first_frame.is_none() && !batch_results.is_empty() && batch_results[0].frame == 0 {
                first_frame = Some(batch_results[0].clone());
            }
        },
    )
    .expect("processing failed");

    // Single-frame: print human-readable summary to stdout
    if total_frames == 1 {
        if let Some(ref r) = first_frame {
            if let Some(ref err) = r.error {
                eprintln!("Error: {err}");
            } else {
                println!("Knot type: {}", r.knot_type);
                if r.knot_start >= 0 {
                    println!(
                        "Knot core: [{}, {}], size = {}",
                        r.knot_start, r.knot_end, r.knot_size
                    );
                }
            }
        }
    }

    writer.flush().expect("flush failed");
    let compute_time = t0.elapsed();

    let knotted = n_knotted;
    let errors = n_errors;
    eprintln!(
        "Done: {} frames, {} knotted, {} unknot, {} errors — {:?}, written to {}",
        total_frames,
        knotted,
        total_frames - knotted - errors,
        errors,
        compute_time,
        output_path
    );
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    run_cli(&args[1..], &args[0]);
}
