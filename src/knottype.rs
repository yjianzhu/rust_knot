use std::collections::BTreeMap;

use crate::alexander_table::AlexanderTable;
use crate::error::{KnotError, Result};
use crate::geometry::{cal_intersection, find_max_span, get_segment, Intersection};
use crate::hull::hull_ends;
use crate::kmt::{kmt_open_chain, kmt_ring};
use crate::point::Point3;
use crate::polynomial::Polynomial;

#[derive(Clone, Debug)]
pub struct KnotOptions {
    /// Apply KMT simplification before computing polynomial (faster but same result).
    pub faster: bool,
    /// Print debug information.
    pub debug: bool,
}

impl Default for KnotOptions {
    fn default() -> Self {
        KnotOptions {
            faster: false,
            debug: false,
        }
    }
}

/// Identify the knot type of an open chain.
///
/// Pipeline: KMT simplification → hull endpoint extension → Alexander polynomial → table lookup.
pub fn get_knottype(
    points: &[Point3],
    table: &AlexanderTable,
    options: &KnotOptions,
) -> Result<String> {
    let mut working = points.to_vec();

    if options.faster {
        kmt_open_chain(&mut working);
    }

    if options.debug && options.faster {
        for p in &working {
            eprintln!("1\t{:10.5}\t{:10.5}\t{:10.5}", p[0], p[1], p[2]);
        }
    }

    if working.len() > 4 {
        match hull_ends(&working) {
            Some((start_ext, end_ext)) => {
                working.insert(0, start_ext);
                working.push(end_ext);
            }
            None => {
                // Fallback: close the chain by connecting back to start
                if !working.is_empty() {
                    let first = working[0];
                    working.push(first);
                }
            }
        }
    }

    let result = compute_knottype_from_chain(&working, table, options.debug)?;

    if options.debug {
        eprintln!("knot type: {result}");
    }

    Ok(result)
}

/// Identify the knot type of a ring (closed) chain.
pub fn get_knottype_ring(
    points: &[Point3],
    table: &AlexanderTable,
    options: &KnotOptions,
) -> Result<String> {
    let mut working = points.to_vec();

    if options.faster {
        kmt_ring(&mut working);
    }

    if working.is_empty() {
        return Ok("1".to_string());
    }

    // Close the ring by appending the first point
    let first = working[0];
    working.push(first);

    let result = compute_knottype_from_chain(&working, table, options.debug)?;

    if options.debug && result != "1" {
        eprintln!("ring knot type: {result}");
    }

    Ok(result)
}

/// Core algorithm: compute Alexander polynomial from a chain of points
/// and look up the knot type in the table.
///
/// Ports: knottype.cpp get_knottype_by_matrix_open_internal (lines 136-285)
fn compute_knottype_from_chain(
    points: &[Point3],
    table: &AlexanderTable,
    debug: bool,
) -> Result<String> {
    let n = points.len();
    if n < 3 {
        return Ok("1".to_string());
    }

    // Reorient for best 2D projection
    let mut pts = points.to_vec();
    find_max_span(&mut pts);

    // Compute all pairwise intersections (O(n^2))
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

    // Build crossing map and segment list
    // crossing: (i, j) -> (crossing_index, sign)
    // where i < j always
    let mut crossing: BTreeMap<(usize, usize), (usize, i32)> = BTreeMap::new();
    let mut segment: Vec<usize> = Vec::new();
    let mut count_crossing = 0usize;

    for i in 0..num_segs {
        // Collect all valid intersections for segment i, sorted by |param|
        let mut one_line: Vec<(usize, f64)> = Vec::new();
        for j in 0..num_segs {
            if intersections[i][j].valid {
                one_line.push((j, intersections[i][j].param.abs()));
            }
        }
        one_line.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for &(j, _) in &one_line {
            let key = if i < j { (i, j) } else { (j, i) };

            if !crossing.contains_key(&key) {
                let is_over = if intersections[i][j].param > 0.0 {
                    1
                } else {
                    -1
                };
                crossing.insert(key, (count_crossing, is_over));
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

    // Build Alexander matrix
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
        one_line.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let mut offset = 0usize;
        for &(j, _) in &one_line {
            let key = if i < j { (i, j) } else { (j, i) };
            let (crossing_idx, sign) = crossing[&key];

            if sign == 1 {
                if !intersections[i][j].seg1_above {
                    // Under-crossing arc: -1 at current segment, +t at next
                    let i_seg = get_segment(&segment, offset, i);
                    let col1 = i_seg % count_crossing;
                    let col2 = (i_seg + 1) % count_crossing;
                    matrix[crossing_idx][col1] = &matrix[crossing_idx][col1] + &neg_one;
                    matrix[crossing_idx][col2] = &matrix[crossing_idx][col2] + &t_poly;
                    offset += 1;
                } else {
                    // Over-crossing arc: +(1-t)
                    let i_seg = get_segment(&segment, offset, i);
                    let col = i_seg % count_crossing;
                    matrix[crossing_idx][col] = &matrix[crossing_idx][col] + &one_minus_t;
                }
            } else {
                if !intersections[i][j].seg1_above {
                    // Under-crossing arc: +t at current, -1 at next
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

    // Take (n-1)x(n-1) submatrix and compute determinant
    let sub_size = count_crossing - 1;
    let mut sub_matrix = vec![vec![Polynomial::zero(); sub_size]; sub_size];
    for i in 0..sub_size {
        for j in 0..sub_size {
            sub_matrix[i][j] = matrix[i][j].clone();
        }
    }

    let poly = Polynomial::determinant(&sub_matrix);

    // Normalize: divide by t^ldegree, ensure positive leading coefficient
    let normalized = poly.normalize();

    // Look up in table
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
