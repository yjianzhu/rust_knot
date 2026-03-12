use std::env;
use std::io::BufReader;
use std::time::Instant;

use rust_knot::alexander_table::AlexanderTable;
use rust_knot::io::read_data_xyz;
use rust_knot::knotsize::{find_knot_core, KnotSizeOptions};
use rust_knot::knottype::{get_knottype, KnotOptions};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <table_file> <xyz_file> [target_type]", args[0]);
        eprintln!("  table_file:  path to table_knot_Alexander_polynomial.txt");
        eprintln!("  xyz_file:    path to XYZ coordinate file");
        eprintln!("  target_type: (optional) knot type to search for core, e.g. '3_1'");
        std::process::exit(1);
    }

    let table_path = &args[1];
    let xyz_path = &args[2];
    let target_type = args.get(3).map(|s| s.as_str());

    // Load Alexander polynomial table
    let t0 = Instant::now();
    let table = AlexanderTable::from_file(table_path).expect("failed to load Alexander table");
    let table_time = t0.elapsed();
    println!(
        "Loaded Alexander table ({} entries) in {:?}",
        table.len(),
        table_time
    );

    // Read XYZ file
    let t0 = Instant::now();
    let file = std::fs::File::open(xyz_path).expect("failed to open XYZ file");
    let mut reader = BufReader::new(file);
    let points = read_data_xyz(&mut reader).expect("failed to parse XYZ file");
    let read_time = t0.elapsed();
    println!("Read {} points in {:?}", points.len(), read_time);

    // Identify knot type
    let opts = KnotOptions {
        faster: true,
        debug: false,
    };

    let t0 = Instant::now();
    match get_knottype(&points, &table, &opts) {
        Ok(knot_type) => {
            let type_time = t0.elapsed();
            println!("Knot type: {} (computed in {:?})", knot_type, type_time);

            // If target type specified or use detected type, find knot core
            let search_target = target_type.unwrap_or(&knot_type);
            if search_target != "1" {
                let size_opts = KnotSizeOptions {
                    is_ring: false,
                    knottype_options: opts.clone(),
                };

                let t0 = Instant::now();
                match find_knot_core(&points, search_target, &table, &size_opts) {
                    Ok(core) => {
                        let core_time = t0.elapsed();
                        if core.matched {
                            println!(
                                "Knot core: [{}, {}], size = {} (found in {:?})",
                                core.left, core.right, core.size, core_time
                            );
                        } else {
                            println!(
                                "No knot core found for type '{}' ({:?})",
                                search_target, core_time
                            );
                        }
                    }
                    Err(e) => {
                        eprintln!("Error finding knot core: {e}");
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("Error identifying knot type: {e}");
        }
    }
}
