use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::Result;
use crate::knottype::get_knottype;
use crate::point::Point3;

#[derive(Clone, Debug)]
pub struct KnotCoreResult {
    /// Inclusive left index (0-based).
    pub left: i32,
    /// Inclusive right index (0-based).
    pub right: i32,
    /// Number of points in the knot core.
    pub size: i32,
    /// Whether a matching knot core was found.
    pub matched: bool,
    /// The knot type of the found core.
    pub found_type: String,
}

impl Default for KnotCoreResult {
    fn default() -> Self {
        KnotCoreResult {
            left: -1,
            right: -1,
            size: 0,
            matched: false,
            found_type: String::new(),
        }
    }
}

/// Parse a knot name like "3_1" into (crossing_number, index).
fn parse_knot_name(name: &str) -> Option<(u32, u32)> {
    let parts: Vec<&str> = name.split('_').collect();
    if parts.len() != 2 {
        return None;
    }
    let crossing: u32 = parts[0].parse().ok()?;
    let index: u32 = parts[1].parse().ok()?;
    Some((crossing, index))
}

/// Returns true if knot `a` is at least as complex as knot `b`.
/// Bug fix: C++ used character comparison which failed for "10_1" vs "3_1".
fn knot_at_least_as_complex(a: &str, b: &str) -> bool {
    match (parse_knot_name(a), parse_knot_name(b)) {
        (Some((ac, ai)), Some((bc, bi))) => (ac, ai) >= (bc, bi),
        _ => false,
    }
}

/// Check if a sub-window [l, r] matches the target knot type.
fn window_is_type(
    points: &[Point3],
    l: usize,
    r: usize,
    target: &str,
    table: &AlexanderTable,
    config: &KnotConfig,
) -> bool {
    if r < l || r >= points.len() {
        return false;
    }
    let sub: Vec<Point3> = points[l..=r].to_vec();
    match get_knottype(&sub, table, config) {
        Ok(t) => t == target,
        Err(e) => {
            if config.debug {
                eprintln!("debug: window_is_type [{l},{r}] failed: {e}");
            }
            false
        }
    }
}

/// Get knot type for a sub-window as open chain (returns empty string on error).
fn window_knottype_open(
    points: &[Point3],
    l: usize,
    r: usize,
    table: &AlexanderTable,
    config: &KnotConfig,
) -> String {
    if r < l || r >= points.len() {
        return String::new();
    }
    let sub: Vec<Point3> = points[l..=r].to_vec();
    let open_cfg = KnotConfig {
        is_ring: false,
        ..*config
    };
    match get_knottype(&sub, table, &open_cfg) {
        Ok(t) => t,
        Err(e) => {
            if config.debug {
                eprintln!("debug: window_knottype_open [{l},{r}] failed: {e}");
            }
            String::new()
        }
    }
}

/// Find the minimal knot core in a chain of points.
///
/// Uses `config.is_ring` to select open-chain or ring-mode search.
/// Uses `config.num_rotations` for ring rotation count.
pub fn find_knot_core(
    points: &[Point3],
    target: &str,
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<KnotCoreResult> {
    if points.is_empty() {
        return Ok(KnotCoreResult::default());
    }

    if !config.is_ring {
        Ok(find_knot_core_open(points, target, table, config))
    } else {
        Ok(find_knot_core_ring(points, target, table, config))
    }
}

fn find_knot_core_open(
    points: &[Point3],
    target: &str,
    table: &AlexanderTable,
    config: &KnotConfig,
) -> KnotCoreResult {
    let mut res = KnotCoreResult::default();
    let n = points.len() as i32;

    // Use open-chain config for all sub-window checks
    let open_cfg = KnotConfig {
        is_ring: false,
        ..*config
    };

    // Verify full sequence matches target
    if !window_is_type(points, 0, (n - 1) as usize, target, table, &open_cfg) {
        return res;
    }

    // Binary search for minimum right boundary
    let mut left = 0i32;
    let mut right = n - 1;
    while left < right {
        let mid = (left + right) / 2;
        if window_is_type(points, 0, mid as usize, target, table, &open_cfg) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    let mut knot_right = left + 1;

    // Binary search for maximum left boundary
    left = 0;
    right = n - 1;
    while left < right - 1 {
        let mid = (left + right) / 2;
        if window_is_type(
            points,
            mid as usize,
            (n - 1) as usize,
            target,
            table,
            &open_cfg,
        ) {
            left = mid;
        } else {
            right = mid;
        }
    }
    let mut knot_left = left + 1;

    // Handle slipknot case
    if knot_left >= knot_right {
        left = knot_left - 1;
        right = n - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if window_is_type(
                points,
                (knot_left - 1) as usize,
                mid as usize,
                target,
                table,
                &open_cfg,
            ) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        let knot_right_temp = left + 1;

        left = 0;
        right = knot_right - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if window_is_type(
                points,
                mid as usize,
                (knot_right - 1) as usize,
                target,
                table,
                &open_cfg,
            ) {
                left = mid;
            } else {
                right = mid;
            }
        }
        let knot_left_temp = left + 1;

        if knot_right_temp - knot_left > knot_right - knot_left_temp {
            knot_left = knot_left_temp;
        } else {
            knot_right = knot_right_temp;
        }
    }

    // Alternating contraction
    let mut flag_left = 0;
    while knot_left > 1 && knot_right < n && knot_left < knot_right {
        if window_is_type(
            points,
            (knot_left - 1) as usize,
            (knot_right - 1) as usize,
            target,
            table,
            &open_cfg,
        ) {
            break;
        } else if flag_left == 1 {
            knot_right += 1;
            flag_left = 0;
        } else {
            knot_left -= 1;
            flag_left = 1;
        }
    }

    res.left = knot_left - 1;
    res.right = knot_right - 1;
    res.size = knot_right - knot_left + 1;
    res.matched = true;
    res.found_type = target.to_string();
    res
}

fn find_knot_core_ring(
    points: &[Point3],
    target: &str,
    table: &AlexanderTable,
    config: &KnotConfig,
) -> KnotCoreResult {
    let n = points.len() as i32;
    let num_shifts = (config.num_rotations as i32).max(1);
    let shift_size = n / num_shifts;

    let ring_cfg = KnotConfig {
        is_ring: true,
        ..*config
    };
    let open_cfg = KnotConfig {
        is_ring: false,
        ..*config
    };

    struct Candidate {
        left: i32,
        right: i32,
        size: i32,
    }
    let mut candidates: Vec<Candidate> = Vec::new();

    for rot in 0..num_shifts {
        let rotated = rotate_points(points, (rot * shift_size) as usize);

        // Verify with ring detection
        if !window_is_type(&rotated, 0, (n - 1) as usize, target, table, &ring_cfg) {
            continue;
        }

        // Binary search min right boundary (open chain)
        let mut left = 0i32;
        let mut right = n - 1;
        while left < right {
            let mid = (left + right) / 2;
            if window_is_type(&rotated, 0, mid as usize, target, table, &open_cfg) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        let mut knot_right = left + 1;

        // Binary search max left boundary (open chain)
        left = 0;
        right = n - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if window_is_type(
                &rotated,
                mid as usize,
                (n - 1) as usize,
                target,
                table,
                &open_cfg,
            ) {
                left = mid;
            } else {
                right = mid;
            }
        }
        let mut knot_left = left + 1;

        // Handle slipknot
        if knot_left >= knot_right {
            left = knot_left - 1;
            right = n - 1;
            while left < right - 1 {
                let mid = (left + right) / 2;
                if window_is_type(
                    &rotated,
                    (knot_left - 1) as usize,
                    mid as usize,
                    target,
                    table,
                    &open_cfg,
                ) {
                    right = mid;
                } else {
                    left = mid + 1;
                }
            }
            let knot_right_temp = left + 1;

            left = 0;
            right = knot_right - 1;
            while left < right - 1 {
                let mid = (left + right) / 2;
                if window_is_type(
                    &rotated,
                    mid as usize,
                    (knot_right - 1) as usize,
                    target,
                    table,
                    &open_cfg,
                ) {
                    left = mid;
                } else {
                    right = mid;
                }
            }
            let knot_left_temp = left + 1;

            if knot_right_temp - knot_left > knot_right - knot_left_temp {
                knot_left = knot_left_temp;
            } else {
                knot_right = knot_right_temp;
            }
        }

        // Alternating contraction
        let mut flag_left = 0;
        while knot_left > 1 && knot_right < n && knot_left <= knot_right {
            if window_is_type(
                &rotated,
                (knot_left - 1) as usize,
                (knot_right - 1) as usize,
                target,
                table,
                &open_cfg,
            ) {
                break;
            } else if flag_left == 1 {
                knot_right += 1;
                flag_left = 0;
            } else {
                knot_left -= 1;
                flag_left = 1;
            }
        }

        // Right tail contraction
        left = knot_left - 1;
        right = knot_right - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if window_is_type(
                &rotated,
                mid as usize,
                (knot_right - 1) as usize,
                target,
                table,
                &open_cfg,
            ) {
                left = mid;
            } else {
                let temp_type = window_knottype_open(
                    &rotated,
                    mid as usize,
                    (knot_right - 1) as usize,
                    table,
                    config,
                );
                if knot_at_least_as_complex(&temp_type, target) {
                    left = mid;
                } else {
                    right = mid;
                }
            }
        }
        knot_left = left + 1;

        // Left tail contraction
        left = knot_left - 1;
        right = knot_right - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if window_is_type(
                &rotated,
                (knot_left - 1) as usize,
                mid as usize,
                target,
                table,
                &open_cfg,
            ) {
                right = mid;
            } else {
                let temp_type = window_knottype_open(
                    &rotated,
                    (knot_left - 1) as usize,
                    mid as usize,
                    table,
                    config,
                );
                if knot_at_least_as_complex(&temp_type, target) {
                    right = mid;
                } else {
                    left = mid;
                }
            }
        }
        knot_right = right + 1;

        // Map back to original indices
        let offset = rot * shift_size;
        let orig_left = (knot_left - 1 + offset) % n;
        let orig_right = (knot_right - 1 + offset) % n;
        candidates.push(Candidate {
            left: orig_left,
            right: orig_right,
            size: knot_right - knot_left + 1,
        });
    }

    if candidates.is_empty() {
        return KnotCoreResult::default();
    }

    let best = candidates
        .iter()
        .filter(|c| c.size > 0)
        .min_by_key(|c| c.size)
        .unwrap_or(&candidates[0]);

    KnotCoreResult {
        left: best.left,
        right: best.right,
        size: best.size,
        matched: true,
        found_type: target.to_string(),
    }
}

fn rotate_points(points: &[Point3], k: usize) -> Vec<Point3> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let k = k % n;
    let mut rotated = Vec::with_capacity(n);
    for i in 0..n {
        rotated.push(points[(i + k) % n]);
    }
    rotated
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_knot_name() {
        assert_eq!(parse_knot_name("3_1"), Some((3, 1)));
        assert_eq!(parse_knot_name("10_1"), Some((10, 1)));
        assert_eq!(parse_knot_name("1"), None);
    }

    #[test]
    fn test_knot_complexity_comparison() {
        assert!(knot_at_least_as_complex("10_1", "3_1"));
        assert!(!knot_at_least_as_complex("3_1", "10_1"));
        assert!(knot_at_least_as_complex("3_1", "3_1"));
    }

    #[test]
    fn test_rotate_points() {
        let pts = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let rotated = rotate_points(&pts, 1);
        assert_eq!(rotated[0], [2.0, 0.0, 0.0]);
        assert_eq!(rotated[1], [3.0, 0.0, 0.0]);
        assert_eq!(rotated[2], [1.0, 0.0, 0.0]);
    }
}
