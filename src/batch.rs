use std::io::BufRead;

use rayon::prelude::*;

use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::KnotError;
use crate::io::XyzFrameIter;
use crate::knotsize::find_knot_core;
use crate::knottype::get_knottype;
use crate::point::Point3;

/// Result of processing a single frame.
#[derive(Clone, Debug)]
pub struct FrameResult {
    pub frame: usize,
    pub knot_type: String,
    pub knot_start: i32,
    pub knot_end: i32,
    pub knot_size: i32,
    pub error: Option<String>,
    /// Non-fatal warnings (e.g. knot core search failed but type was identified).
    pub warnings: Vec<String>,
}

/// Process a single frame: identify knot type, then find knot core.
///
/// If `target_type` is provided, the core search uses it instead of the
/// detected knot type.
pub fn process_frame(
    frame: usize,
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
) -> FrameResult {
    let mut warnings = Vec::new();
    let (knot_type, unknown_poly) = match get_knottype(points, table, config) {
        Ok(t) => (t, false),
        Err(KnotError::NotFound(poly)) => {
            warnings.push(format!("knot polynomial not found in table: {poly}"));
            (poly, true)
        }
        Err(e) => {
            return FrameResult {
                frame,
                knot_type: String::new(),
                knot_start: -1,
                knot_end: -1,
                knot_size: 0,
                error: Some(e.to_string()),
                warnings: Vec::new(),
            };
        }
    };

    let search_target = if knot_type == "1" {
        None
    } else if unknown_poly {
        target_type
    } else {
        Some(target_type.unwrap_or(&knot_type))
    };

    let (start, end, size) = if let Some(search) = search_target {
        match find_knot_core(points, search, table, config) {
            Ok(core) if core.matched => (core.left, core.right, core.size),
            Ok(_) => {
                warnings.push(format!("knot core not found for type '{search}'"));
                (-1, -1, 0)
            }
            Err(e) => {
                warnings.push(format!("knot core search failed: {e}"));
                (-1, -1, 0)
            }
        }
    } else {
        (-1, -1, 0)
    };

    FrameResult {
        frame,
        knot_type,
        knot_start: start,
        knot_end: end,
        knot_size: size,
        error: None,
        warnings,
    }
}

/// Process all frames in parallel using rayon, returning results in frame order.
pub fn process_frames_parallel(
    frames: &[Vec<Point3>],
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
) -> Vec<FrameResult> {
    frames
        .par_iter()
        .enumerate()
        .map(|(i, points)| process_frame(i, points, table, config, target_type))
        .collect()
}

/// Default batch size for streaming processing.
const DEFAULT_BATCH_SIZE: usize = 64;

/// Streaming batch processor: reads `batch_size` frames at a time from an XYZ
/// reader, processes each batch in parallel, and calls `on_batch` with results.
///
/// Memory usage is bounded to `batch_size * points_per_frame`.
/// Returns total number of frames processed.
pub fn process_frames_streaming<R, F>(
    reader: R,
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
    batch_size: Option<usize>,
    mut on_batch: F,
) -> crate::error::Result<usize>
where
    R: BufRead,
    F: FnMut(&[FrameResult]),
{
    let batch_size = batch_size.unwrap_or(DEFAULT_BATCH_SIZE);
    let mut frame_offset = 0usize;
    let mut batch = Vec::with_capacity(batch_size);
    let iter = XyzFrameIter::new(reader);

    for result in iter {
        let points = result?;
        batch.push(points);

        if batch.len() >= batch_size {
            let results = process_batch(&batch, frame_offset, table, config, target_type);
            on_batch(&results);
            frame_offset += batch.len();
            batch.clear();
        }
    }

    // Process remaining frames
    if !batch.is_empty() {
        let results = process_batch(&batch, frame_offset, table, config, target_type);
        on_batch(&results);
        frame_offset += batch.len();
    }

    Ok(frame_offset)
}

fn process_batch(
    batch: &[Vec<Point3>],
    offset: usize,
    table: &AlexanderTable,
    config: &KnotConfig,
    target_type: Option<&str>,
) -> Vec<FrameResult> {
    batch
        .par_iter()
        .enumerate()
        .map(|(i, points)| process_frame(offset + i, points, table, config, target_type))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point3;

    fn straight_line(n: usize) -> Vec<Point3> {
        (0..n).map(|i| [i as f64, 0.0, 0.0]).collect()
    }

    #[test]
    fn test_process_frame_unknot() {
        let table_data = "0_1\t1\n3_1\t-1+1*t\n";
        let table = AlexanderTable::from_reader(std::io::Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let pts = straight_line(50);
        let result = process_frame(0, &pts, &table, &config, None);
        assert_eq!(result.knot_type, "1");
        assert_eq!(result.knot_start, -1);
        assert!(result.error.is_none());
    }

    #[test]
    fn test_process_frames_parallel_ordering() {
        let table_data = "0_1\t1\n";
        let table = AlexanderTable::from_reader(std::io::Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let frames: Vec<Vec<Point3>> = (0..5).map(|_| straight_line(50)).collect();
        let results = process_frames_parallel(&frames, &table, &config, None);
        assert_eq!(results.len(), 5);
        for (i, r) in results.iter().enumerate() {
            assert_eq!(r.frame, i);
        }
    }
}
