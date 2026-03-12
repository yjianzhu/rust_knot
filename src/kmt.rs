use crate::geometry::{cal_normals, judge_triangle};
use crate::point::Point3;

/// KMT simplification for open chains.
/// Iteratively removes interior points whose removal doesn't change the knot type
/// (i.e., the triangle formed by the point and its neighbors doesn't intersect
/// any other segment of the chain).
///
/// Ports: knottype.cpp KMT_open_chain (lines 53-92)
pub fn kmt_open_chain(points: &mut Vec<Point3>) {
    loop {
        let initial_count = points.len();

        let mut i = 1;
        while i < points.len().saturating_sub(1) {
            let plane = cal_normals(&points[i - 1], &points[i], &points[i + 1]);

            match plane {
                None => {
                    // Collinear points — safe to remove
                    points.remove(i);
                    continue; // don't increment i
                }
                Some(plane) => {
                    // Check if removing this point would change topology
                    let mut intersects = false;
                    for j in 0..points.len().saturating_sub(1) {
                        if j == i - 1 || j == i {
                            continue;
                        }
                        if judge_triangle(
                            &points[i - 1],
                            &points[i],
                            &points[i + 1],
                            &plane,
                            &points[j],
                            &points[j + 1],
                        ) {
                            intersects = true;
                            break;
                        }
                    }
                    if !intersects {
                        points.remove(i);
                        continue;
                    }
                }
            }
            i += 1;
        }

        if points.len() == initial_count {
            break;
        }
    }
}

/// KMT simplification for ring (closed) chains.
/// Same logic as open chain but with modular indexing.
///
/// Ports: knottype.cpp KMT_ring (lines 95-133)
pub fn kmt_ring(points: &mut Vec<Point3>) {
    loop {
        let initial_count = points.len();
        if initial_count < 3 {
            break;
        }

        let mut i: usize = 0;
        let mut checked = 0;
        while checked < points.len() && points.len() >= 3 {
            let n = points.len();
            let prev = if i == 0 { n - 1 } else { i - 1 };
            let curr = i % n;
            let next = (i + 1) % n;

            let plane = cal_normals(&points[prev], &points[curr], &points[next]);

            match plane {
                None => {
                    points.remove(curr);
                    if curr < points.len() {
                        // Stay at current position since elements shifted
                    } else {
                        i = 0;
                    }
                    checked = 0;
                    continue;
                }
                Some(plane) => {
                    let mut intersects = false;
                    let n = points.len();
                    // Check all segments except those adjacent to the triangle
                    for offset in 2..n - 1 {
                        let j = (i + offset) % n;
                        let j_next = (j + 1) % n;
                        if judge_triangle(
                            &points[prev],
                            &points[curr],
                            &points[next],
                            &plane,
                            &points[j],
                            &points[j_next],
                        ) {
                            intersects = true;
                            break;
                        }
                    }
                    if !intersects {
                        points.remove(curr);
                        if points.is_empty() {
                            break;
                        }
                        if curr >= points.len() {
                            i = 0;
                        }
                        checked = 0;
                        continue;
                    }
                }
            }
            i = (i + 1) % points.len().max(1);
            checked += 1;
        }

        if points.len() == initial_count {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmt_open_straight_line() {
        // A straight line should be simplified to just endpoints
        let mut points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        kmt_open_chain(&mut points);
        assert_eq!(points.len(), 2);
    }

    #[test]
    fn test_kmt_open_preserves_endpoints() {
        let mut points = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let first = points.first().unwrap().clone();
        let last = points.last().unwrap().clone();
        kmt_open_chain(&mut points);
        assert_eq!(*points.first().unwrap(), first);
        assert_eq!(*points.last().unwrap(), last);
    }
}
