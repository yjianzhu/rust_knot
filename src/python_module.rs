use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};

use ndarray::{Array2, Array3};
use numpy::{IntoPyArray, PyArray2, PyArray3, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::exceptions::{PyFileNotFoundError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::alexander_table::AlexanderTable;
use crate::config::KnotConfig;
use crate::error::KnotError;
use crate::io::{read_data_xyz_frames, write_data_xyz};
use crate::kmt::{kmt_open_chain, kmt_ring};
use crate::knotsize::find_knot_core;
use crate::knottype::get_knottype;
use crate::point::Point3;

const DEFAULT_TABLE_FILENAME: &str = "table_knot_Alexander_polynomial.txt";

fn to_py_runtime_error(message: impl Into<String>) -> PyErr {
    PyRuntimeError::new_err(message.into())
}

fn resolve_table_path() -> Option<PathBuf> {
    if let Ok(path) = env::var("PYTHONKNOT_ALEXANDER_TABLE") {
        return Some(PathBuf::from(path));
    }

    if let Ok(home) = env::var("HOME") {
        let default_path = PathBuf::from(home)
            .join(".local")
            .join("share")
            .join("rust_knot")
            .join(DEFAULT_TABLE_FILENAME);
        if default_path.exists() {
            return Some(default_path);
        }
    }

    None
}

fn load_table() -> PyResult<AlexanderTable> {
    if let Some(path) = resolve_table_path() {
        if !path.exists() {
            return Err(PyFileNotFoundError::new_err(format!(
                "Alexander table file not found: {}",
                path.display()
            )));
        }
        AlexanderTable::builtin_with_file(path.to_string_lossy().as_ref())
            .map_err(|e| to_py_runtime_error(e.to_string()))
    } else {
        Ok(AlexanderTable::builtin())
    }
}

fn config_from_chain_type(chain_type: &str) -> PyResult<KnotConfig> {
    match chain_type {
        "ring" => Ok(KnotConfig {
            is_ring: true,
            faster: true,
            ..KnotConfig::default()
        }),
        "open" => Ok(KnotConfig {
            is_ring: false,
            faster: true,
            ..KnotConfig::default()
        }),
        _ => Err(PyRuntimeError::new_err(
            "Invalid chain_type. Expected 'ring' or 'open'.",
        )),
    }
}

fn points_from_array2(array: PyReadonlyArray2<f64>) -> PyResult<Vec<Point3>> {
    let view = array.as_array();
    let shape = view.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyRuntimeError::new_err(
            "Expected array with shape (N, 3) or (F, N, 3).",
        ));
    }

    let mut points = Vec::with_capacity(shape[0]);
    for row in view.outer_iter() {
        points.push([row[0], row[1], row[2]]);
    }
    Ok(points)
}

fn frames_from_array3(array: PyReadonlyArray3<f64>) -> PyResult<Vec<Vec<Point3>>> {
    let view = array.as_array();
    let shape = view.shape();
    if shape.len() != 3 || shape[2] != 3 {
        return Err(PyRuntimeError::new_err(
            "Expected array with shape (N, 3) or (F, N, 3).",
        ));
    }

    let mut frames = Vec::with_capacity(shape[0]);
    for frame_view in view.outer_iter() {
        let mut points = Vec::with_capacity(shape[1]);
        for row in frame_view.outer_iter() {
            points.push([row[0], row[1], row[2]]);
        }
        frames.push(points);
    }

    Ok(frames)
}

fn frames_from_path(path: &Path) -> PyResult<Vec<Vec<Point3>>> {
    if !path.exists() {
        return Err(PyFileNotFoundError::new_err(format!(
            "Cannot open file: {}",
            path.display()
        )));
    }
    let file = File::open(path).map_err(|e| to_py_runtime_error(e.to_string()))?;
    let mut reader = BufReader::new(file);
    read_data_xyz_frames(&mut reader).map_err(|e| to_py_runtime_error(e.to_string()))
}

fn parse_pdb_model_index(line: &str, fallback: usize) -> usize {
    line.get(5..)
        .unwrap_or("")
        .trim()
        .split_whitespace()
        .next()
        .and_then(|s| s.parse::<usize>().ok())
        .map(|n| n.saturating_sub(1))
        .unwrap_or(fallback)
}

fn parse_pdb_coord(
    line: &str,
    start: usize,
    end: usize,
    field: &str,
    line_no: usize,
) -> Result<f64, String> {
    line.get(start..end)
        .unwrap_or("")
        .trim()
        .parse::<f64>()
        .map_err(|e| format!("Invalid {field} at PDB line {line_no}: {e}"))
}

fn frames_from_pdb_path(path: &Path, atom_filter: &str) -> PyResult<Vec<Vec<Point3>>> {
    if !path.exists() {
        return Err(PyFileNotFoundError::new_err(format!(
            "Cannot open file: {}",
            path.display()
        )));
    }

    enum AtomFilter {
        All,
        Ca,
    }

    let atom_filter = match atom_filter.to_ascii_lowercase().as_str() {
        "all" => AtomFilter::All,
        "ca" => AtomFilter::Ca,
        _ => {
            return Err(PyValueError::new_err(
                "Invalid atom_filter. Expected 'all' or 'ca'.",
            ))
        }
    };

    let file = File::open(path).map_err(|e| to_py_runtime_error(e.to_string()))?;
    let mut reader = BufReader::new(file);
    let mut frames: Vec<Vec<Point3>> = Vec::new();
    let mut current_frame = 0usize;
    let mut line = String::new();
    let mut line_no = 0usize;
    let mut saw_atom = false;

    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| to_py_runtime_error(e.to_string()))?;
        if n == 0 {
            break;
        }
        line_no += 1;

        if line.starts_with("MODEL") {
            current_frame = parse_pdb_model_index(&line, frames.len());
            if current_frame >= frames.len() {
                frames.resize_with(current_frame + 1, Vec::new);
            }
            continue;
        }
        if line.starts_with("ENDMDL") {
            continue;
        }
        if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
            continue;
        }

        if frames.is_empty() {
            frames.push(Vec::new());
            current_frame = 0;
        }
        if current_frame >= frames.len() {
            frames.resize_with(current_frame + 1, Vec::new);
        }

        if matches!(atom_filter, AtomFilter::Ca) {
            let atom_name = line.get(12..16).unwrap_or("").trim();
            if atom_name != "CA" {
                continue;
            }
        }

        let x = parse_pdb_coord(&line, 30, 38, "x", line_no).map_err(to_py_runtime_error)?;
        let y = parse_pdb_coord(&line, 38, 46, "y", line_no).map_err(to_py_runtime_error)?;
        let z = parse_pdb_coord(&line, 46, 54, "z", line_no).map_err(to_py_runtime_error)?;
        frames[current_frame].push([x, y, z]);
        saw_atom = true;
    }

    if !saw_atom {
        return Err(PyRuntimeError::new_err(
            "No ATOM/HETATM records found in PDB file.",
        ));
    }
    Ok(frames)
}

fn frames_from_input(input_data: &PyAny) -> PyResult<Vec<Vec<Point3>>> {
    if let Ok(path) = input_data.extract::<String>() {
        return frames_from_path(Path::new(&path));
    }
    if let Ok(array3) = input_data.extract::<PyReadonlyArray3<f64>>() {
        return frames_from_array3(array3);
    }
    if let Ok(array2) = input_data.extract::<PyReadonlyArray2<f64>>() {
        return Ok(vec![points_from_array2(array2)?]);
    }
    Err(PyRuntimeError::new_err(
        "Expected input to be XYZ filename or NumPy array with shape (N, 3) / (F, N, 3).",
    ))
}

fn frames_to_pyarray<'py>(py: Python<'py>, frames: &[Vec<Point3>]) -> PyResult<&'py PyArray3<f64>> {
    let n_frames = frames.len();
    let max_atoms = frames.iter().map(|f| f.len()).max().unwrap_or(0);
    // Ragged frames are padded to a dense (F, N, 3) tensor with NaN sentinels.
    let mut data = vec![f64::NAN; n_frames * max_atoms * 3];

    for (frame_idx, frame) in frames.iter().enumerate() {
        for (atom_idx, point) in frame.iter().enumerate() {
            let base = (frame_idx * max_atoms + atom_idx) * 3;
            data[base] = point[0];
            data[base + 1] = point[1];
            data[base + 2] = point[2];
        }
    }

    let array = Array3::from_shape_vec((n_frames, max_atoms, 3), data)
        .map_err(|e| to_py_runtime_error(e.to_string()))?;
    Ok(array.into_pyarray(py))
}

fn points_to_pyarray<'py>(py: Python<'py>, points: &[Point3]) -> PyResult<&'py PyArray2<f64>> {
    let mut data = Vec::with_capacity(points.len() * 3);
    for point in points {
        data.extend_from_slice(point);
    }
    let array = Array2::from_shape_vec((points.len(), 3), data)
        .map_err(|e| to_py_runtime_error(e.to_string()))?;
    Ok(array.into_pyarray(py))
}

fn classify_frame(
    points: &[Point3],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> Result<String, String> {
    match get_knottype(points, table, config) {
        Ok(t) => Ok(t),
        Err(KnotError::NotFound(poly)) => Ok(poly),
        Err(e) => Err(e.to_string()),
    }
}

fn classify_frames(
    frames: &[Vec<Point3>],
    table: &AlexanderTable,
    config: &KnotConfig,
    num_threads: Option<usize>,
) -> Result<Vec<String>, String> {
    let results: Vec<Result<String, String>> = if let Some(threads) = num_threads {
        if threads == 0 {
            return Err("threads must be >= 1".to_string());
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?;
        pool.install(|| {
            frames
                .par_iter()
                .map(|points| classify_frame(points, table, config))
                .collect()
        })
    } else {
        frames
            .par_iter()
            .map(|points| classify_frame(points, table, config))
            .collect()
    };

    let mut out = Vec::with_capacity(results.len());
    for result in results {
        match result {
            Ok(value) => out.push(value),
            Err(message) => return Err(message),
        }
    }
    Ok(out)
}

fn compute_knot_size(
    frames: &[Vec<Point3>],
    table: &AlexanderTable,
    config: &KnotConfig,
    num_threads: Option<usize>,
) -> Result<(Vec<String>, Vec<Vec<i32>>), String> {
    let compute_one = |points: &[Point3]| -> Result<(String, Vec<i32>), String> {
        let knot_type = classify_frame(points, table, config)?;
        let size = if knot_type == "1" || !knot_type.contains('_') {
            vec![-1, -1, 0]
        } else {
            let core =
                find_knot_core(points, &knot_type, table, config).map_err(|e| e.to_string())?;
            if core.matched {
                vec![core.left, core.right, core.size]
            } else {
                vec![-1, -1, 0]
            }
        };
        Ok((knot_type, size))
    };

    let results: Vec<Result<(String, Vec<i32>), String>> = if let Some(threads) = num_threads {
        if threads == 0 {
            return Err("threads must be >= 1".to_string());
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?;
        pool.install(|| {
            frames
                .par_iter()
                .map(|points| compute_one(points))
                .collect()
        })
    } else {
        frames
            .par_iter()
            .map(|points| compute_one(points))
            .collect()
    };

    let mut knot_types = Vec::with_capacity(results.len());
    let mut knot_sizes = Vec::with_capacity(results.len());
    for result in results {
        match result {
            Ok((knot_type, knot_size)) => {
                knot_types.push(knot_type);
                knot_sizes.push(knot_size);
            }
            Err(message) => return Err(message),
        }
    }

    Ok((knot_types, knot_sizes))
}

#[pyfunction]
/// Read XYZ file into a NumPy array with shape `(frames, atoms, 3)`.
///
/// If frames have different atom counts, shorter frames are padded with `NaN`.
fn read_xyz(py: Python<'_>, filename: &str) -> PyResult<Py<PyArray3<f64>>> {
    let frames = frames_from_path(Path::new(filename))?;
    Ok(frames_to_pyarray(py, &frames)?.to_owned())
}

#[pyfunction(signature = (filename, atom_filter = "all"))]
/// Read PDB file into a NumPy array with shape `(frames, atoms, 3)`.
///
/// If frames have different atom counts, shorter frames are padded with `NaN`.
fn read_pdb(py: Python<'_>, filename: &str, atom_filter: &str) -> PyResult<Py<PyArray3<f64>>> {
    let frames = frames_from_pdb_path(Path::new(filename), atom_filter)?;
    Ok(frames_to_pyarray(py, &frames)?.to_owned())
}

#[pyfunction(signature = (filename, input_data, append = false))]
fn write_xyz(filename: &str, input_data: &PyAny, append: bool) -> PyResult<()> {
    let frames = frames_from_input(input_data)?;
    let mut options = OpenOptions::new();
    options.create(true).write(true);
    if append {
        options.append(true);
    } else {
        options.truncate(true);
    }
    let file = options
        .open(filename)
        .map_err(|e| to_py_runtime_error(e.to_string()))?;
    let mut writer = BufWriter::new(file);

    for frame in &frames {
        write_data_xyz(frame, &mut writer).map_err(|e| to_py_runtime_error(e.to_string()))?;
    }

    Ok(())
}

#[pyfunction(signature = (input_data, chain_type = "ring", threads = None))]
fn knot_type(
    py: Python<'_>,
    input_data: &PyAny,
    chain_type: &str,
    threads: Option<usize>,
) -> PyResult<Vec<String>> {
    if matches!(threads, Some(0)) {
        return Err(PyValueError::new_err("threads must be >= 1"));
    }
    let frames = frames_from_input(input_data)?;
    let table = load_table()?;
    let config = config_from_chain_type(chain_type)?;
    py.allow_threads(|| classify_frames(&frames, &table, &config, threads))
        .map_err(to_py_runtime_error)
}

#[pyfunction(signature = (input_data, chain_type = "ring", threads = None))]
fn knot_size(
    py: Python<'_>,
    input_data: &PyAny,
    chain_type: &str,
    threads: Option<usize>,
) -> PyResult<(Vec<String>, Vec<Vec<i32>>)> {
    if matches!(threads, Some(0)) {
        return Err(PyValueError::new_err("threads must be >= 1"));
    }
    let frames = frames_from_input(input_data)?;
    let table = load_table()?;
    let config = config_from_chain_type(chain_type)?;
    py.allow_threads(|| compute_knot_size(&frames, &table, &config, threads))
        .map_err(to_py_runtime_error)
}

#[pyfunction(signature = (input_data, chain_type = "ring"))]
fn kmt(py: Python<'_>, input_data: &PyAny, chain_type: &str) -> PyResult<Py<PyArray2<f64>>> {
    let mut points = if let Ok(array2) = input_data.extract::<PyReadonlyArray2<f64>>() {
        points_from_array2(array2)?
    } else {
        return Err(PyRuntimeError::new_err(
            "kmt expects a NumPy array with shape (N, 3).",
        ));
    };

    match chain_type {
        "ring" => kmt_ring(&mut points),
        "open" => kmt_open_chain(&mut points),
        _ => {
            return Err(PyRuntimeError::new_err(
                "Invalid chain_type. Expected 'ring' or 'open'.",
            ))
        }
    }

    Ok(points_to_pyarray(py, &points)?.to_owned())
}

#[pymodule]
fn alexander_poly(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(read_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(write_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(knot_type, m)?)?;
    m.add_function(wrap_pyfunction!(knot_size, m)?)?;
    m.add_function(wrap_pyfunction!(kmt, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::sync::Once;

    use super::*;

    fn init_python() {
        static INIT: Once = Once::new();
        INIT.call_once(pyo3::prepare_freethreaded_python);
    }

    fn straight_line(n: usize) -> Vec<Point3> {
        (0..n).map(|i| [i as f64, 0.0, 0.0]).collect()
    }

    #[test]
    fn frames_from_array3_parses_expected_points() {
        init_python();
        Python::with_gil(|py| {
            let data = vec![
                0.0, 1.0, 2.0, 3.0, 4.0, 5.0, // frame 0
                6.0, 7.0, 8.0, 9.0, 10.0, 11.0, // frame 1
            ];
            let array = Array3::from_shape_vec((2, 2, 3), data).unwrap();
            let py_array = array.into_pyarray(py);
            let frames = frames_from_array3(py_array.readonly()).unwrap();
            assert_eq!(frames.len(), 2);
            assert_eq!(frames[0], vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]);
            assert_eq!(frames[1], vec![[6.0, 7.0, 8.0], [9.0, 10.0, 11.0]]);
        });
    }

    #[test]
    fn frames_from_array3_rejects_invalid_shape() {
        init_python();
        Python::with_gil(|py| {
            let array = Array3::from_shape_vec((1, 2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();
            let py_array = array.into_pyarray(py);
            let err = frames_from_array3(py_array.readonly()).unwrap_err();
            assert!(err.to_string().contains("Expected array with shape"));
        });
    }

    #[test]
    fn classify_frames_defaults_to_parallel_when_threads_none() {
        let table_data = "0_1\t1\n";
        let table = AlexanderTable::from_reader(Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let frames: Vec<Vec<Point3>> = (0..4).map(|_| straight_line(50)).collect();

        let knot_types = classify_frames(&frames, &table, &config, None).unwrap();
        assert_eq!(knot_types, vec!["1".to_string(); 4]);
    }

    #[test]
    fn compute_knot_size_defaults_to_parallel_when_threads_none() {
        let table_data = "0_1\t1\n";
        let table = AlexanderTable::from_reader(Cursor::new(table_data)).unwrap();
        let config = KnotConfig::default();
        let frames: Vec<Vec<Point3>> = (0..3).map(|_| straight_line(50)).collect();

        let (knot_types, sizes) = compute_knot_size(&frames, &table, &config, None).unwrap();
        assert_eq!(knot_types, vec!["1".to_string(); 3]);
        assert_eq!(sizes, vec![vec![-1, -1, 0]; 3]);
    }
}
