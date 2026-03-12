use crate::point::Point3;
use chull::ConvexHullWrapper;

const HULL_PLANE_EPSILON: f64 = 5e-3;
const EXTEND_FACTOR: f64 = 100.0;

fn vector_length(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Compute extended endpoints for an open chain by pushing them outward
/// through the convex hull. This ensures the closure doesn't interfere
/// with the knot structure.
///
/// Returns `Some((start_ext, end_ext))` — two new points to prepend/append.
/// Returns `None` if the hull computation fails.
///
/// Ports: hull.cpp hull_ends_internal
pub fn hull_ends(points: &[Point3]) -> Option<(Point3, Point3)> {
    if points.len() < 4 {
        return None;
    }

    let pts: Vec<Vec<f64>> = points.iter().map(|p| vec![p[0], p[1], p[2]]).collect();
    let hull = ConvexHullWrapper::try_new(&pts, None).ok()?;
    let (vertices, indices) = hull.vertices_indices();
    if indices.is_empty() {
        return None;
    }

    // Compute center of all points
    let n = points.len() as f64;
    let mut center = [0.0; 3];
    for p in points {
        for d in 0..3 {
            center[d] += p[d];
        }
    }
    for d in 0..3 {
        center[d] /= n;
    }

    // Compute max radius from center
    let max_radius = points
        .iter()
        .map(|p| {
            let delta = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            vector_length(&delta)
        })
        .fold(0.0f64, f64::max);

    if max_radius <= f64::EPSILON {
        return None;
    }

    // Build list of triangular faces with normals and offsets
    struct Face {
        normal: [f64; 3],
        offset: f64,
    }

    let mut faces = Vec::new();
    // chull returns flat list of vertex indices, every 3 consecutive indices = 1 face
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let a = [
            vertices[tri[0]][0],
            vertices[tri[0]][1],
            vertices[tri[0]][2],
        ];
        let b = [
            vertices[tri[1]][0],
            vertices[tri[1]][1],
            vertices[tri[1]][2],
        ];
        let c = [
            vertices[tri[2]][0],
            vertices[tri[2]][1],
            vertices[tri[2]][2],
        ];

        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal = crate::geometry::cross_product(&e1, &e2);
        let len = vector_length(&normal);
        if len < f64::EPSILON {
            continue;
        }
        let normal = [normal[0] / len, normal[1] / len, normal[2] / len];
        let offset = normal[0] * a[0] + normal[1] * a[1] + normal[2] * a[2];
        faces.push(Face { normal, offset });
    }

    if faces.is_empty() {
        return None;
    }

    let first_pt = points.first()?;
    let last_pt = points.last()?;

    // Find closest face to first point and last point
    let dist_to_face = |pt: &Point3, face: &Face| -> f64 {
        let d =
            pt[0] * face.normal[0] + pt[1] * face.normal[1] + pt[2] * face.normal[2] - face.offset;
        d.abs()
    };

    let first_face = faces
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| dist_to_face(first_pt, a).total_cmp(&dist_to_face(first_pt, b)))
        .map(|(i, _)| i)?;

    let last_face = faces
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| dist_to_face(last_pt, a).total_cmp(&dist_to_face(last_pt, b)))
        .map(|(i, _)| i)?;

    // Compute extension vectors
    let compute_extension = |pt: &Point3, face_idx: usize| -> [f64; 3] {
        let face = &faces[face_idx];
        let dist = dist_to_face(pt, face);

        if dist < HULL_PLANE_EPSILON {
            // Point is on the hull face — extend away from center
            [pt[0] - center[0], pt[1] - center[1], pt[2] - center[2]]
        } else {
            // Project point onto face plane and use normal direction
            let dot = pt[0] * face.normal[0] + pt[1] * face.normal[1] + pt[2] * face.normal[2]
                - face.offset;
            [
                -dot * face.normal[0],
                -dot * face.normal[1],
                -dot * face.normal[2],
            ]
        }
    };

    let ext_first = compute_extension(first_pt, first_face);
    let ext_last = compute_extension(last_pt, last_face);

    let len_first = vector_length(&ext_first);
    let len_last = vector_length(&ext_last);

    if len_first <= f64::EPSILON || len_last <= f64::EPSILON {
        return None;
    }

    let scale_first = max_radius / len_first;
    let scale_last = max_radius / len_last;

    let point1 = [
        first_pt[0] + EXTEND_FACTOR * scale_first * ext_first[0],
        first_pt[1] + EXTEND_FACTOR * scale_first * ext_first[1],
        first_pt[2] + EXTEND_FACTOR * scale_first * ext_first[2],
    ];

    let point2 = [
        last_pt[0] + EXTEND_FACTOR * scale_last * ext_last[0],
        last_pt[1] + EXTEND_FACTOR * scale_last * ext_last[1],
        last_pt[2] + EXTEND_FACTOR * scale_last * ext_last[2],
    ];

    Some((point1, point2))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_finite_point(p: &Point3) -> bool {
        p[0].is_finite() && p[1].is_finite() && p[2].is_finite()
    }

    #[test]
    fn hull_ends_returns_none_for_short_chain() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        assert!(hull_ends(&points).is_none());
    }

    #[test]
    fn hull_ends_returns_none_for_degenerate_geometry() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        assert!(hull_ends(&points).is_none());
    }

    #[test]
    fn hull_ends_returns_extended_endpoints_for_valid_chain() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.25, 0.5, 1.2],
            [0.0, 0.7, 0.6],
        ];

        let (start_ext, end_ext) = hull_ends(&points).expect("expected non-degenerate hull");
        assert!(is_finite_point(&start_ext));
        assert!(is_finite_point(&end_ext));
        assert_ne!(start_ext, points[0]);
        assert_ne!(end_ext, points[points.len() - 1]);
    }
}
