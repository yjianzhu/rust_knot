use std::io::{BufRead, Write};

use crate::error::{KnotError, Result};
use crate::point::Point3;

/// Read points from XYZ format.
/// Format:
/// ```text
/// N
/// (comment line)
/// atom_type  x  y  z
/// ...
/// ```
pub fn read_data_xyz<R: BufRead>(reader: &mut R) -> Result<Vec<Point3>> {
    let mut line = String::new();

    // Read N
    reader.read_line(&mut line).map_err(KnotError::Io)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| KnotError::DataParse(format!("invalid point count: {e}")))?;

    // Skip comment line
    line.clear();
    reader.read_line(&mut line).map_err(KnotError::Io)?;

    let mut points = Vec::with_capacity(n);
    for i in 0..n {
        line.clear();
        reader.read_line(&mut line).map_err(KnotError::Io)?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(KnotError::DataParse(format!(
                "line {}: expected 'type x y z', got '{}'",
                i + 3,
                line.trim()
            )));
        }
        let x: f64 = parts[1]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad x at line {}: {e}", i + 3)))?;
        let y: f64 = parts[2]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad y at line {}: {e}", i + 3)))?;
        let z: f64 = parts[3]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad z at line {}: {e}", i + 3)))?;
        points.push([x, y, z]);
    }

    Ok(points)
}

/// Read all frames from an XYZ stream.
/// Each frame uses the standard XYZ layout:
/// ```text
/// N
/// (comment line)
/// atom_type  x  y  z
/// ...
/// ```
/// Returns an empty vector for empty input.
pub fn read_data_xyz_frames<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<Point3>>> {
    let mut frames = Vec::new();
    let mut line = String::new();

    loop {
        // Find the next non-empty frame header line (point count).
        let n = loop {
            line.clear();
            let bytes = reader.read_line(&mut line).map_err(KnotError::Io)?;
            if bytes == 0 {
                return Ok(frames);
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            break trimmed
                .parse::<usize>()
                .map_err(|e| KnotError::DataParse(format!("invalid point count: {e}")))?;
        };

        // Skip comment line
        line.clear();
        if reader.read_line(&mut line).map_err(KnotError::Io)? == 0 {
            return Err(KnotError::DataParse(
                "unexpected EOF after XYZ point count".into(),
            ));
        }

        let mut points = Vec::with_capacity(n);
        for i in 0..n {
            line.clear();
            if reader.read_line(&mut line).map_err(KnotError::Io)? == 0 {
                return Err(KnotError::DataParse(format!(
                    "unexpected EOF in XYZ data line {}",
                    i + 1
                )));
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(KnotError::DataParse(format!(
                    "line {}: expected 'type x y z', got '{}'",
                    i + 3,
                    line.trim()
                )));
            }
            let x: f64 = parts[1]
                .parse()
                .map_err(|e| KnotError::DataParse(format!("bad x at line {}: {e}", i + 3)))?;
            let y: f64 = parts[2]
                .parse()
                .map_err(|e| KnotError::DataParse(format!("bad y at line {}: {e}", i + 3)))?;
            let z: f64 = parts[3]
                .parse()
                .map_err(|e| KnotError::DataParse(format!("bad z at line {}: {e}", i + 3)))?;
            points.push([x, y, z]);
        }

        frames.push(points);
    }
}

/// Read points from LAMMPS dump format.
/// Expects: 3 header lines, then N, then 5 more header lines, then data.
pub fn read_data_lammps<R: BufRead>(reader: &mut R) -> Result<Vec<Point3>> {
    let mut line = String::new();

    // Skip first 3 lines
    for _ in 0..3 {
        line.clear();
        if reader.read_line(&mut line).map_err(KnotError::Io)? == 0 {
            return Err(KnotError::DataParse(
                "unexpected EOF in LAMMPS header".into(),
            ));
        }
    }

    // Read N
    line.clear();
    reader.read_line(&mut line).map_err(KnotError::Io)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| KnotError::DataParse(format!("invalid atom count: {e}")))?;

    // Skip 5 more header lines (with the remaining part of the line after N + 5 more)
    // Actually, the C++ reads one line that has N, then skips 6 lines total
    // Let's match C++ which does: read N, then getline 6 times
    for _ in 0..5 {
        line.clear();
        if reader.read_line(&mut line).map_err(KnotError::Io)? == 0 {
            return Err(KnotError::DataParse(
                "unexpected EOF in LAMMPS header".into(),
            ));
        }
    }

    let mut points = Vec::with_capacity(n);
    for _ in 0..n {
        line.clear();
        reader.read_line(&mut line).map_err(KnotError::Io)?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            return Err(KnotError::DataParse(format!(
                "LAMMPS data line too short: '{}'",
                line.trim()
            )));
        }
        // Format: id atom_type x y z
        let x: f64 = parts[2]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad x: {e}")))?;
        let y: f64 = parts[3]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad y: {e}")))?;
        let z: f64 = parts[4]
            .parse()
            .map_err(|e| KnotError::DataParse(format!("bad z: {e}")))?;
        points.push([x, y, z]);
    }

    Ok(points)
}

/// Write points in XYZ format.
pub fn write_data_xyz<W: Write>(points: &[Point3], writer: &mut W) -> Result<()> {
    writeln!(writer, "{}", points.len()).map_err(KnotError::Io)?;
    writeln!(writer).map_err(KnotError::Io)?;
    for p in points {
        writeln!(writer, "1\t{:10.5}\t{:10.5}\t{:10.5}", p[0], p[1], p[2])
            .map_err(KnotError::Io)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_xyz() {
        let data = "3\ncomment\n1 1.0 2.0 3.0\n1 4.0 5.0 6.0\n1 7.0 8.0 9.0\n";
        let mut reader = Cursor::new(data);
        let points = read_data_xyz(&mut reader).unwrap();
        assert_eq!(points.len(), 3);
        assert_eq!(points[0], [1.0, 2.0, 3.0]);
        assert_eq!(points[2], [7.0, 8.0, 9.0]);
    }

    #[test]
    fn test_write_read_roundtrip() {
        let points = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut buf = Vec::new();
        write_data_xyz(&points, &mut buf).unwrap();

        let mut reader = Cursor::new(buf);
        let read_back = read_data_xyz(&mut reader).unwrap();
        assert_eq!(read_back.len(), 2);
        // Check approximate equality due to formatting
        for (a, b) in points.iter().zip(read_back.iter()) {
            for d in 0..3 {
                assert!((a[d] - b[d]).abs() < 0.001);
            }
        }
    }

    #[test]
    fn test_read_xyz_frames() {
        let data = concat!(
            "2\n",
            "frame0\n",
            "1 1.0 2.0 3.0\n",
            "1 4.0 5.0 6.0\n",
            "3\n",
            "frame1\n",
            "1 7.0 8.0 9.0\n",
            "1 10.0 11.0 12.0\n",
            "1 13.0 14.0 15.0\n"
        );
        let mut reader = Cursor::new(data);
        let frames = read_data_xyz_frames(&mut reader).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].len(), 2);
        assert_eq!(frames[1].len(), 3);
        assert_eq!(frames[0][0], [1.0, 2.0, 3.0]);
        assert_eq!(frames[1][2], [13.0, 14.0, 15.0]);
    }

    #[test]
    fn test_read_xyz_frames_with_blank_lines_and_empty_input() {
        let data = "\n\n1\nframe0\n1 1.0 2.0 3.0\n\n";
        let mut reader = Cursor::new(data);
        let frames = read_data_xyz_frames(&mut reader).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0][0], [1.0, 2.0, 3.0]);

        let mut empty_reader = Cursor::new("");
        let empty_frames = read_data_xyz_frames(&mut empty_reader).unwrap();
        assert!(empty_frames.is_empty());
    }
}
