use std::collections::BTreeMap;

use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::{KnotError, Result};
use crate::geometry::{cal_intersection, find_max_span, get_segment, Intersection};
use crate::hull::hull_ends;
use crate::kmt::{kmt_open_chain, kmt_ring};
use crate::point::Point3;
use crate::polynomial::Polynomial;

/// Identify the knot type of a chain (open or ring, based on `config.is_ring`).
pub fn get_knottype(
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String> {
    if config.is_ring {
        get_knottype_ring(points, table, config)
    } else {
        get_knottype_open(points, table, config)
    }
}

/// Identify the knot type of an open chain.
///
/// Pipeline: KMT simplification -> hull endpoint extension -> Alexander polynomial -> table lookup.
fn get_knottype_open(
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String> {
    let mut working = points.to_vec();

    if config.faster {
        kmt_open_chain(&mut working);
    }

    if config.debug && config.faster {
        for p in &working {
            eprintln!("1\t{:10.5}\t{:10.5}\t{:10.5}", p[0], p[1], p[2]);
        }
    }

    if working.len() > 4 {
        match hull_ends(&working, config.hull_plane_epsilon, config.extend_factor) {
            Some((start_ext, end_ext)) => {
                working.insert(0, start_ext);
                working.push(end_ext);
            }
            None => {
                if config.debug {
                    eprintln!(
                        "warning: hull_ends failed, using infinity-closure fallback (no head-tail segment)"
                    );
                }
            }
        }
    }

    let result = compute_knottype_from_chain(&working, table, config.debug)?;

    if config.debug {
        eprintln!("knot type: {result}");
    }

    Ok(result)
}

/// Identify the knot type of a ring (closed) chain.
fn get_knottype_ring(
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String> {
    let mut working = points.to_vec();

    if config.faster {
        kmt_ring(&mut working);
    }

    if working.is_empty() {
        return Ok("1".to_string());
    }

    working.push(working[0]);

    let result = compute_knottype_from_chain(&working, table, config.debug)?;
    if config.debug && result != "1" {
        eprintln!("ring knot type: {result}");
    }
    Ok(result)
}

/// Core algorithm: compute Alexander polynomial from a chain of points
/// and look up the knot type in the table.
///
/// Ports: knottype.cpp get_knottype_by_matrix_open_internal (lines 136-285)
#[allow(clippy::needless_range_loop)]
fn compute_knottype_from_chain(
    points: &[Point3],
    table: &AlexanderTable,
    debug: bool,
) -> Result<String> {
    let n = points.len();
    if n < 3 {
        return Ok("1".to_string());
    }

    let mut pts = points.to_vec();
    find_max_span(&mut pts);

    let num_segs = n - 1;
    let mut intersections = vec![vec![Intersection::invalid(); num_segs]; num_segs];
    for i in 0..num_segs {
        for j in 0..num_segs {
            if i == j {
                continue;
            }
            intersections[i][j] = cal_intersection(&pts[i], &pts[i + 1], &pts[j], &pts[j + 1]);
        }
    }

    let mut crossing: BTreeMap<(usize, usize), (usize, i32)> = BTreeMap::new();
    let mut segment: Vec<usize> = Vec::new();
    let mut count_crossing = 0usize;

    for i in 0..num_segs {
        let mut one_line: Vec<(usize, f64)> = Vec::new();
        for j in 0..num_segs {
            if intersections[i][j].valid {
                one_line.push((j, intersections[i][j].param.abs()));
            }
        }
        one_line.sort_by(|a, b| a.1.total_cmp(&b.1));

        for &(j, _) in &one_line {
            let key = if i < j { (i, j) } else { (j, i) };

            if let std::collections::btree_map::Entry::Vacant(e) = crossing.entry(key) {
                let is_over = if intersections[i][j].param > 0.0 {
                    1
                } else {
                    -1
                };
                e.insert((count_crossing, is_over));
                count_crossing += 1;
            }

            if !intersections[i][j].seg1_above {
                segment.push(i);
            }
        }
    }

    if count_crossing == 0 {
        return Ok("1".to_string());
    }

    let t_poly = Polynomial::t();
    let one_minus_t = Polynomial::one_minus_t();
    let neg_one = Polynomial::from_constant(-1);

    let mut matrix = vec![vec![Polynomial::zero(); count_crossing]; count_crossing];

    for i in 0..num_segs {
        let mut one_line: Vec<(usize, f64)> = Vec::new();
        for j in 0..num_segs {
            if intersections[i][j].valid {
                one_line.push((j, intersections[i][j].param.abs()));
            }
        }
        one_line.sort_by(|a, b| a.1.total_cmp(&b.1));

        let mut offset = 0usize;
        for &(j, _) in &one_line {
            let key = if i < j { (i, j) } else { (j, i) };
            let (crossing_idx, sign) = crossing[&key];

            if sign == 1 {
                if !intersections[i][j].seg1_above {
                    let i_seg = get_segment(&segment, offset, i);
                    let col1 = i_seg % count_crossing;
                    let col2 = (i_seg + 1) % count_crossing;
                    matrix[crossing_idx][col1] = &matrix[crossing_idx][col1] + &neg_one;
                    matrix[crossing_idx][col2] = &matrix[crossing_idx][col2] + &t_poly;
                    offset += 1;
                } else {
                    let i_seg = get_segment(&segment, offset, i);
                    let col = i_seg % count_crossing;
                    matrix[crossing_idx][col] = &matrix[crossing_idx][col] + &one_minus_t;
                }
            } else {
                if !intersections[i][j].seg1_above {
                    let i_seg = get_segment(&segment, offset, i);
                    let col1 = i_seg % count_crossing;
                    let col2 = (i_seg + 1) % count_crossing;
                    matrix[crossing_idx][col1] = &matrix[crossing_idx][col1] + &t_poly;
                    matrix[crossing_idx][col2] = &matrix[crossing_idx][col2] + &neg_one;
                    offset += 1;
                } else {
                    let i_seg = get_segment(&segment, offset, i);
                    let col = i_seg % count_crossing;
                    matrix[crossing_idx][col] = &matrix[crossing_idx][col] + &one_minus_t;
                }
            }
        }
    }

    let sub_size = count_crossing - 1;
    let mut sub_matrix = vec![vec![Polynomial::zero(); sub_size]; sub_size];
    for i in 0..sub_size {
        for j in 0..sub_size {
            sub_matrix[i][j] = matrix[i][j].clone();
        }
    }

    let poly = Polynomial::determinant(&sub_matrix);
    let normalized = poly.normalize();

    match table.lookup(&normalized) {
        Some(name) => Ok(name.to_string()),
        None => {
            if debug {
                eprintln!("knot type not found, with polynomial as: {normalized}");
            }
            Err(KnotError::NotFound)
        }
    }
}
