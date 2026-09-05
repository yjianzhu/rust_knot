# rust_knot

[中文](README.md) | [English](README_EN.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22316028.svg)](https://doi.org/10.5281/zenodo.22316028)

A knot topology analysis library written in Rust. It identifies knot types using Alexander polynomials and locates the smallest knotted core.

This project is a complete rewrite of the C++ implementation (`knottype.cpp`, `knotsize.cpp`, `hull.cpp`, and `my_function.cpp`). It removes the GiNaC symbolic mathematics dependency, fixes known bugs, and provides a thread-safe, purely functional design.

## Algorithm

```
Input point sequence (3D coordinates)
    │
    ▼
KMT simplification (remove redundant points without changing topology)
    │
    ▼
Convex-hull endpoint extension (open chains: push both ends away from the knot)
    │
    ▼
O(n²) crossing detection (XY projection + Z-order comparison)
    │
    ▼
Alexander matrix construction (coefficients in Z[t])
    │
    ▼
Fraction-free Bareiss determinant (replaces GiNaC)
    │
    ▼
Polynomial table lookup
    │
    ▼
Knot type ("3_1", "4_1", ...)
    │
    ▼
Binary search for the smallest knotted core (knotsize)
    │
    ▼
KnotCoreResult { left, right, size }
```

## Project Structure

```
rust_knot/
├── Cargo.toml
├── src/
│   ├── lib.rs               # Crate entry point and public API re-exports
│   ├── main.rs              # CLI entry point (direct argument mode)
│   ├── config.rs            # KnotConfig — unified configuration
│   ├── batch.rs             # Batch processing: process_frame / process_frames_parallel / process_frames_streaming
│   ├── point.rs             # Point3 = [f64; 3], EPSILON = 1e-7
│   ├── error.rs             # KnotError variants (Io, PolynomialParse, DataParse, HullFailed, NotFound, EmptyChain)
│   ├── polynomial.rs        # Polynomial<i64> arithmetic + Bareiss determinant (~600 lines, replaces GiNaC)
│   ├── geometry.rs          # Crossing detection, triangle intersection, normal calculation, axis redirection
│   ├── alexander_table.rs   # Alexander polynomial lookup table (embedded ≤9 crossings + external files)
│   ├── kmt.rs               # KMT simplification for open and closed chains
│   ├── hull.rs              # Convex-hull endpoint extension based on the chull crate
│   ├── knottype.rs           # Knot identification: crossings → Alexander matrix → polynomial → lookup
│   ├── knotsize.rs           # Binary search for the knotted core
│   └── io.rs                 # XYZ/LAMMPS I/O + lazy XyzFrameIter frame iterator
├── tests/
│   └── integration.rs        # Six end-to-end integration tests
└── .github/
    └── workflows/
        ├── ci.yml            # GitHub Actions checks for pushes and pull requests
        └── release.yml       # Tag-triggered builds for three platforms
```

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| `chull` | 0.2 | 3D convex hull (QuickHull algorithm) |
| `thiserror` | 2 | Error type derivation |
| `rayon` | 1.10 | Parallel batch processing across frames |
| `approx` | 0.5 | Floating-point comparisons in tests (dev-dependency) |

All dependencies are pure Rust with no C/FFI bindings, enabling zero-configuration cross-compilation.

**Key design decision:** no symbolic mathematics library is used. Alexander matrix elements are limited to `{0, 1, -1, t, 1-t}`, so the project implements a lightweight `Polynomial { coeffs: Vec<i64> }` type (about 600 lines) together with the fraction-free Bareiss determinant algorithm. All computation stays within the integer polynomial ring `Z[t]`.

## Configuration (`KnotConfig`)

All hyperparameters and mode flags are managed by the `KnotConfig` struct:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `is_ring` | `bool` | `false` | Open chain (`false`) or closed ring (`true`) |
| `faster` | `bool` | `false` | Enable KMT simplification for faster execution without changing results |
| `debug` | `bool` | `false` | Write debugging information to stderr |
| `hull_plane_epsilon` | `f64` | `5e-3` | Convex-hull face tolerance; larger values are more permissive |
| `extend_factor` | `f64` | `100.0` | Endpoint extension scale; larger values push endpoints farther away |
| `num_rotations` | `u32` | `4` | Number of rotation searches in ring mode; more rotations improve precision but cost more time |

`EPSILON = 1e-7` in `point.rs` is the global geometric zero tolerance used for crossing detection and normal calculation. It is a low-level numerical constant and is therefore not part of `KnotConfig`.

## Build and Test

```bash
cd rust_knot

# Build
cargo build --release

# Run all 51 tests (45 unit + 6 integration)
cargo test

# Build the CLI tool
cargo build --release
```

## Usage

### CLI

```bash
# Minimal usage: embedded table (≤9 crossings), no external file required
cargo run -- input.xyz
# Or use the release binary
./target/release/rust_knot input.xyz

# Closed-ring mode
./target/release/rust_knot input.xyz --ring

# Add an external table (>9 crossings)
./target/release/rust_knot input.xyz --table extended_table.txt

# Multi-frame file with a custom batch size and output path
./target/release/rust_knot trajectory.xyz --ring --batch 128 --threads 8 --output result.txt

# All arguments
./target/release/rust_knot <xyz_file> [target_type] [--table <path>] [--ring] [--no-fast] [--debug] [--output <path>] [--batch <size>] [--threads <n>]
```

Output file `knot_index.txt`:

```
# frame	knottype	knot_start	knot_end	knot_size
0	3_1	114	145	32
1	3_1	107	158	52
2	3_1	108	171	64
...
```

### As a Library Dependency

Add the dependency to your project's `Cargo.toml`:

```toml
[dependencies]
# Option 1: local path (for development on the same machine)
rust_knot = { path = "../rust_knot" }

# Option 2: Git repository
rust_knot = { git = "https://github.com/yjianzhu/rust_knot.git" }

# Option 3: specific version or branch
rust_knot = { git = "https://github.com/yjianzhu/rust_knot.git", tag = "v0.2.4" }
```

#### Public API

| Type/Function | Description |
|---------------|-------------|
| `Point3` | Point type alias for `[f64; 3]` |
| `AlexanderTable` | Alexander polynomial lookup table |
| `KnotConfig` | Unified configuration struct |
| `get_knottype(&[Point3], &AlexanderTable, &KnotConfig) -> Result<String>` | Identify the knot type |
| `find_knot_core(&[Point3], &str, &AlexanderTable, &KnotConfig) -> Result<KnotCoreResult>` | Locate the smallest knotted core |
| `KnotCoreResult` | Core result containing `left`, `right`, `size`, `matched`, and `found_type` |
| `process_frame` / `process_frames_parallel` / `process_frames_streaming` | Batch-processing interfaces |
| `FrameResult` | Result for one frame in a batch |

#### Basic Example

```rust
use rust_knot::{AlexanderTable, KnotConfig, get_knottype, find_knot_core};

// Option 1: embedded table (≤9 crossings), no external file required
let table = AlexanderTable::builtin();

// Option 2: embedded table merged with an external file (covering ≥10 crossings)
let table = AlexanderTable::builtin_with_file("extended_table.txt")?;

// Option 3: external file only
let table = AlexanderTable::from_file("full_table.txt")?;

// Configuration
let config = KnotConfig {
    faster: true,
    is_ring: false,
    ..KnotConfig::default()
};

// Identify the knot type
let knot_type = get_knottype(&points, &table, &config)?;  // e.g. "3_1"

// Locate the smallest knotted core
let core = find_knot_core(&points, &knot_type, &table, &config)?;
println!("Core interval: [{}, {}], size: {}", core.left, core.right, core.size);
```

#### Multi-Frame Batch Processing

```rust
use rust_knot::batch::{process_frames_parallel, process_frames_streaming};
use std::io::BufReader;
use std::fs::File;

// Option 1: load everything and process in parallel (small files)
let results = process_frames_parallel(&frames, &table, &config, None);

// Option 2: streaming batches with bounded memory usage (large files)
let reader = BufReader::new(File::open("huge_trajectory.xyz")?);
let total = process_frames_streaming(
    reader, &table, &config, None,
    Some(128),           // batch_size
    |batch_results| {    // callback for each batch
        for r in batch_results {
            println!("frame {}: {}", r.frame, r.knot_type);
        }
    },
)?;
```

## Alexander Polynomial Table

### Embedded Table

The program embeds all 86 Alexander polynomials with crossing number ≤9, from the unknot and 3₁ through 9₄₉. This covers most common knot types. Use `AlexanderTable::builtin()` without any external file. Both sign variants are stored in the hash table so either form can be identified.

### External Extension Table

To identify knots with ≥10 crossings, load an external table using `--table` or `builtin_with_file()`. The format is:

```
knot_name	polynomial
10_1	4-9*t+4*t^2
10_2	-2+7*t-2*t^2
...
```

Embedded and external tables are merged and deduplicated automatically. When one polynomial maps to several knots, candidates are sorted by crossing number and the simplest one takes precedence.

## Cross-Platform Releases

GitHub Actions automatically builds releases for three platforms:

| Platform | Artifact |
|----------|----------|
| Linux x86_64 | `rust_knot-linux-x86_64` |
| Windows x86_64 | `rust_knot-windows-x86_64.exe` |
| macOS ARM | `rust_knot-macos-arm64` |

Publish a new version:

```bash
git tag v0.2.4
git push origin v0.2.4
# CI builds the binaries and publishes the GitHub release automatically
```

## Agent Skills

The repository includes two skills ready for use with Codex or Claude Code:

| Skill | Purpose |
|-------|---------|
| `rust-knot-cli` | Analyze XYZ trajectories with the global `rust_knot` CLI, construct arguments, and troubleshoot output |
| `export-knot-xyz` | Export initial XYZ conformations for a specified knot type and chain length from a local SQLite database |

### Install the Skills

```bash
git clone https://github.com/yjianzhu/rust_knot.git
cd rust_knot

# Codex
mkdir -p ~/.codex/skills
ln -sfn "$PWD/skills/rust-knot-cli" ~/.codex/skills/rust-knot-cli
ln -sfn "$PWD/skills/export-knot-xyz" ~/.codex/skills/export-knot-xyz

# Claude Code
mkdir -p ~/.claude/skills
ln -sfn "$PWD/skills/rust-knot-cli" ~/.claude/skills/rust-knot-cli
ln -sfn "$PWD/skills/export-knot-xyz" ~/.claude/skills/export-knot-xyz
```

### Install CLI Dependencies

`rust-knot-cli` uses these paths by default:

```text
~/.local/bin/rust_knot
~/.local/share/rust_knot/table_knot_Alexander_polynomial.txt
```

`export-knot-xyz` uses these paths by default:

```text
~/.local/bin/export_xyz
~/.local/share/rust_knot/knots_data.db
```

Download the `rust_knot` binary for your platform from GitHub Releases and place it in `~/.local/bin/`. For example, on Linux:

```bash
mkdir -p ~/.local/bin ~/.local/share/rust_knot
curl -L -o ~/.local/bin/rust_knot \
  https://github.com/yjianzhu/rust_knot/releases/latest/download/rust_knot-linux-x86_64
chmod +x ~/.local/bin/rust_knot
```

`export_xyz` and `knots_data.db` are packaged in the `export-knot-xyz-linux-x86_64.tar.gz` release asset:

```bash
curl -L -o /tmp/export-knot-xyz-linux-x86_64.tar.gz \
  https://github.com/yjianzhu/rust_knot/releases/latest/download/export-knot-xyz-linux-x86_64.tar.gz
tar -xzf /tmp/export-knot-xyz-linux-x86_64.tar.gz -C /tmp

cp /tmp/export-knot-xyz/bin/export_xyz ~/.local/bin/
cp /tmp/export-knot-xyz/share/rust_knot/knots_data.db ~/.local/share/rust_knot/
chmod +x ~/.local/bin/export_xyz
```

## Improvements over the C++ Version

### Bug Fixes

| Bug | C++ Location | Rust Fix |
|-----|--------------|----------|
| Character comparison used for knot complexity | `knotsize.cpp:259,279` — `temp[0] > target[0]` | Numeric comparison with `parse_knot_name("10_1") → (10,1)` |
| Unknown polynomial returned an empty string | `knottype.cpp:281` | Return `Err(KnotError::NotFound)` |
| NaN input caused a panic | `partial_cmp().unwrap()` | NaN-safe sorting with `total_cmp()` |
| Global GiNaC mutex blocked parallel execution | `knottype.cpp:193` | No shared mutable state; frames run in parallel with Rayon |
| Convex-hull failures were silently ignored | `hull.cpp catch(...)` | Return `Option` and emit a warning in debug mode |
| Malformed lookup-table entries were silently skipped | No validation | Strict mode reports errors with line numbers by default; lenient mode provides backward compatibility |

### Architectural Improvements

- **No GiNaC dependency:** integer polynomials and the Bareiss determinant fully replace symbolic computation.
- **Embedded polynomial table:** knots with ≤9 crossings are built in, so no external file is required at startup.
- **Streaming batch processing:** `XyzFrameIter` reads lazily and processes parallel batches with constant memory usage.
- **Thread safety:** all core functions are pure and use no global state.
- **Cross-platform:** pure-Rust dependencies and automated GitHub Actions releases for three platforms.
- **Strongly typed error handling:** the `KnotError` enum covers every error path.

## Test Coverage

51 tests (45 unit + 6 integration):

| Module | Tests | Coverage |
|--------|-------|----------|
| `polynomial` | 16 | Arithmetic, determinants, parsing, normalization, exact division, equality, and hashing |
| `alexander_table` | 7 | Table parsing, sign lookup, ambiguity handling, strict/lenient modes, embedded table, and merging |
| `geometry` | 4 | Crossing detection, normal calculation, and cross products |
| `hull` | 3 | Endpoint extension, degenerate geometry, and short chains |
| `kmt` | 2 | Straight-line simplification and endpoint preservation |
| `knotsize` | 3 | Knot-name parsing, complexity comparison, and point-sequence rotation |
| `io` | 4 | XYZ I/O, multi-frame reading, and round-trip consistency |
| `batch` | 2 | Single-frame processing and ordered parallel results |
| **Integration tests** | **6** | Trefoil ring, open figure-eight, unknot, core localization, and ambiguous table lookup |

## License

Same as the parent project.
