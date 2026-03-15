use std::io::BufReader;

use rust_knot::alexander_table::AlexanderTable;
use rust_knot::batch::process_frame;
use rust_knot::config::KnotConfig;
use rust_knot::io::read_data_xyz;
use rust_knot::knotsize::find_knot_core;
use rust_knot::knottype::get_knottype;

fn load_table() -> AlexanderTable {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/table_knot_Alexander_polynomial.txt"
    );
    AlexanderTable::from_file(path).expect("failed to load table")
}

fn load_xyz(name: &str) -> Vec<rust_knot::Point3> {
    let full = format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), name);
    let file = std::fs::File::open(&full).unwrap_or_else(|e| panic!("failed to open {full}: {e}"));
    let mut reader = BufReader::new(file);
    read_data_xyz(&mut reader).expect("failed to parse XYZ")
}

// ── Trefoil ring (3_1) ──

#[test]
fn test_trefoil_ring_identification() {
    let table = load_table();
    let points = load_xyz("L300_knot3_1_ring.xyz");

    let config = KnotConfig {
        is_ring: true,
        faster: true,
        ..KnotConfig::default()
    };

    let knot_type = get_knottype(&points, &table, &config).expect("identification failed");
    assert_eq!(knot_type, "3_1", "expected trefoil knot");
}

#[test]
fn test_trefoil_ring_core_found() {
    let table = load_table();
    let points = load_xyz("L300_knot3_1_ring.xyz");

    let config = KnotConfig {
        is_ring: true,
        faster: true,
        ..KnotConfig::default()
    };

    let core = find_knot_core(&points, "3_1", &table, &config).expect("core search failed");
    assert!(core.matched, "expected to find knot core");
    assert!(core.size > 0, "core size should be positive");
    assert!(
        core.size < points.len() as i32,
        "core should be smaller than full chain"
    );
}

// ── Figure-eight open (4_1) ──

#[test]
fn test_figure_eight_open_identification() {
    let table = load_table();
    let points = load_xyz("L400_knot4_1_open.xyz");

    let config = KnotConfig {
        is_ring: false,
        faster: true,
        ..KnotConfig::default()
    };

    let knot_type = get_knottype(&points, &table, &config).expect("identification failed");
    assert_eq!(knot_type, "4_1", "expected figure-eight knot");
}

#[test]
fn test_figure_eight_open_core_found() {
    let table = load_table();
    let points = load_xyz("L400_knot4_1_open.xyz");

    let config = KnotConfig {
        is_ring: false,
        faster: true,
        ..KnotConfig::default()
    };

    let core = find_knot_core(&points, "4_1", &table, &config).expect("core search failed");
    assert!(core.matched, "expected to find knot core");
    assert!(core.size > 0);
    assert!(core.size < points.len() as i32);
}

// ── Unknot detection ──

#[test]
fn test_straight_line_is_unknot() {
    let table = load_table();

    // Generate a simple straight line — should be unknot
    let points: Vec<rust_knot::Point3> = (0..50).map(|i| [i as f64, 0.0, 0.0]).collect();

    let config = KnotConfig {
        is_ring: false,
        faster: true,
        ..KnotConfig::default()
    };

    let knot_type = get_knottype(&points, &table, &config).expect("identification failed");
    assert_eq!(knot_type, "1", "straight line should be unknot");
}

// ── Ambiguous polynomial lookup ──

#[test]
fn test_ambiguous_lookup_returns_simplest() {
    let table = load_table();

    // 8_3 and 10_1 share the same polynomial; lookup should return 8_3 (lower crossing)
    let poly = rust_knot::polynomial::parse_polynomial("4-9*t+4*t^2").unwrap();
    assert_eq!(table.lookup(&poly), Some("8_3"));

    let all = table.lookup_all(&poly);
    assert!(all.len() >= 2, "expected at least 2 candidates");
    assert_eq!(all[0], "8_3");
}

#[test]
fn test_notfound_uses_polynomial_as_knot_type() {
    let table_data = "0_1\t1\n";
    let table = AlexanderTable::from_reader(std::io::Cursor::new(table_data)).unwrap();
    let points = load_xyz("L300_knot3_1_ring.xyz");

    let config = KnotConfig {
        is_ring: true,
        faster: true,
        ..KnotConfig::default()
    };

    let result = process_frame(0, &points, &table, &config, None);
    assert!(
        result.error.is_none(),
        "notfound should not be treated as fatal error"
    );
    assert_ne!(result.knot_type, "1");
    assert!(
        result.knot_type.contains('t'),
        "expected polynomial string in knot_type, got '{}'",
        result.knot_type
    );
    assert_eq!(result.knot_start, -1);
    assert_eq!(result.knot_end, -1);
    assert_eq!(result.knot_size, 0);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("polynomial not found in table")),
        "expected warning about missing polynomial, got {:?}",
        result.warnings
    );
}
