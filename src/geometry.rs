use crate::point::{Point3, EPSILON};

/// Crossing information from 2D projected segment intersection.
#[derive(Clone, Debug)]
pub struct Intersection {
    /// Signed parameter on segment 1. Sign encodes crossing orientation.
    pub param: f64,
    /// Parameter on segment 2 [0, 1].
    pub param2: f64,
    /// True if segment 1 is above segment 2 at the crossing point.
    pub seg1_above: bool,
    /// True if the intersection is valid (segments cross in their interiors).
    pub valid: bool,
}

impl Intersection {
    pub fn invalid() -> Self {
        Intersection {
            param: 0.0,
            param2: 0.0,
            seg1_above: false,
            valid: false,
        }
    }
}

/// Calculate the 2D-projected intersection of two 3D line segments.
/// Uses Cramer's rule on the XY projection, then compares Z values
/// to determine which segment is above.
///
/// Ports: my_function.cpp cal_interSection (lines 195-243)
pub fn cal_intersection(
    seg1_start: &Point3,
    seg1_end: &Point3,
    seg2_start: &Point3,
    seg2_end: &Point3,
) -> Intersection {
    let a1 = seg1_end[0] - seg1_start[0];
    let b1 = seg2_start[0] - seg2_end[0];
    let a2 = seg1_end[1] - seg1_start[1];
    let b2 = seg2_start[1] - seg2_end[1];

    let det = a1 * b2 - b1 * a2;
    if det.abs() < EPSILON {
        return Intersection::invalid();
    }

    let c1 = seg2_start[0] - seg1_start[0];
    let c2 = seg2_start[1] - seg1_start[1];

    let u = (c1 * b2 - b1 * c2) / det;
    let v = (a1 * c2 - c1 * a2) / det;

    if u <= EPSILON || u >= 1.0 - EPSILON || v <= EPSILON || v >= 1.0 - EPSILON {
        return Intersection::invalid();
    }

    let z1 = seg1_start[2] + u * (seg1_end[2] - seg1_start[2]);
    let z2 = seg2_start[2] + v * (seg2_end[2] - seg2_start[2]);

    let seg1_above = z1 > z2;
    let det_positive = det > 0.0;

    let param = if seg1_above {
        if det_positive {
            u
        } else {
            -u
        }
    } else if det_positive {
        -u
    } else {
        u
    };

    Intersection {
        param,
        param2: v,
        seg1_above,
        valid: true,
    }
}

/// Cross product of two 3-vectors.
pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the plane equation (normal + offset) for three points.
/// Returns `Some([nx, ny, nz, d])` where nx*x + ny*y + nz*z + d = 0,
/// or `None` if the points are collinear.
///
/// Ports: my_function.cpp cal_normals (lines 264-278)
pub fn cal_normals(i: &Point3, j: &Point3, k: &Point3) -> Option<[f64; 4]> {
    let vij = [j[0] - i[0], j[1] - i[1], j[2] - i[2]];
    let vik = [k[0] - i[0], k[1] - i[1], k[2] - i[2]];
    let n = cross_product(&vij, &vik);
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();

    if len.abs() < EPSILON {
        return None;
    }

    let nx = n[0] / len;
    let ny = n[1] / len;
    let nz = n[2] / len;
    let d = -(i[0] * nx + i[1] * ny + i[2] * nz);

    Some([nx, ny, nz, d])
}

/// Test if a line segment (line_1 to line_2) intersects triangle (a, b, c)
/// with known plane equation.
///
/// Uses Moller-Trumbore-like ray-triangle intersection.
/// Returns true if the segment passes through the triangle interior.
///
/// Ports: my_function.cpp judge_triangle (lines 161-193)
pub fn judge_triangle(
    a: &Point3,
    b: &Point3,
    c: &Point3,
    plane: &[f64; 4],
    line_1: &Point3,
    line_2: &Point3,
) -> bool {
    let dis1 = line_1[0] * plane[0] + line_1[1] * plane[1] + line_1[2] * plane[2] + plane[3];
    let dis2 = line_2[0] * plane[0] + line_2[1] * plane[1] + line_2[2] * plane[2] + plane[3];

    if dis1 * dis2 > EPSILON {
        return false;
    }
    if dis1.abs() < EPSILON && dis2.abs() < EPSILON {
        return true;
    }

    let big_t = [line_1[0] - a[0], line_1[1] - a[1], line_1[2] - a[2]];
    let d = [
        line_2[0] - line_1[0],
        line_2[1] - line_1[1],
        line_2[2] - line_1[2],
    ];
    let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];

    let m = cross_product(&d, &e2);
    let det = m[0] * e1[0] + m[1] * e1[1] + m[2] * e1[2];

    if det.abs() < EPSILON {
        return false;
    }

    let kk = cross_product(&big_t, &e1);
    let t_param = (kk[0] * e2[0] + kk[1] * e2[1] + kk[2] * e2[2]) / det;
    let u = (m[0] * big_t[0] + m[1] * big_t[1] + m[2] * big_t[2]) / det;
    let v = (kk[0] * d[0] + kk[1] * d[1] + kk[2] * d[2]) / det;

    !(u <= 0.0 || v <= 0.0 || u + v >= 1.0 || t_param >= 1.0 || t_param <= 0.0)
}

/// Reorient coordinates so the axis with the smallest span becomes Z.
/// This maximizes the 2D projection area for crossing detection.
///
/// Ports: my_function.cpp find_max_span (lines 8-42)
pub fn find_max_span(points: &mut [Point3]) {
    if points.is_empty() {
        return;
    }

    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];

    for p in points.iter() {
        for d in 0..3 {
            if p[d] < min[d] {
                min[d] = p[d];
            }
            if p[d] > max[d] {
                max[d] = p[d];
            }
        }
    }

    let span = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];

    if span[0] < span[2] && span[0] < span[1] {
        // Swap X and Z
        for p in points.iter_mut() {
            p.swap(0, 2);
        }
    } else if span[1] < span[0] && span[1] < span[2] {
        // Swap Y and Z
        for p in points.iter_mut() {
            p.swap(1, 2);
        }
    }
}

/// Map a line segment index to the arc-segment index.
///
/// Ports: my_function.cpp get_segment (lines 144-154)
pub fn get_segment(segment: &[usize], offset: usize, line_segment: usize) -> usize {
    for (i, &seg) in segment.iter().enumerate() {
        if line_segment <= seg {
            return offset + i;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross_product(&a, &b);
        assert!((c[0]).abs() < EPSILON);
        assert!((c[1]).abs() < EPSILON);
        assert!((c[2] - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_cal_normals_basic() {
        let i = [0.0, 0.0, 0.0];
        let j = [1.0, 0.0, 0.0];
        let k = [0.0, 1.0, 0.0];
        let plane = cal_normals(&i, &j, &k).unwrap();
        // Normal should be (0, 0, 1) or (0, 0, -1)
        assert!(plane[0].abs() < EPSILON);
        assert!(plane[1].abs() < EPSILON);
        assert!((plane[2].abs() - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_cal_normals_collinear() {
        let i = [0.0, 0.0, 0.0];
        let j = [1.0, 0.0, 0.0];
        let k = [2.0, 0.0, 0.0];
        assert!(cal_normals(&i, &j, &k).is_none());
    }

    #[test]
    fn test_intersection_crossing() {
        // Two segments that cross in XY plane
        let s1a = [0.0, 0.0, 0.0];
        let s1b = [1.0, 1.0, 1.0];
        let s2a = [0.0, 1.0, 0.0];
        let s2b = [1.0, 0.0, 0.0];
        let ix = cal_intersection(&s1a, &s1b, &s2a, &s2b);
        assert!(ix.valid);
    }

    #[test]
    fn test_intersection_parallel() {
        let s1a = [0.0, 0.0, 0.0];
        let s1b = [1.0, 0.0, 0.0];
        let s2a = [0.0, 1.0, 0.0];
        let s2b = [1.0, 1.0, 0.0];
        let ix = cal_intersection(&s1a, &s1b, &s2a, &s2b);
        assert!(!ix.valid);
    }
}
