use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter};
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

const DEFAULT_TABLE_PATH: &str =
    "/home/yongjian/.local/share/rust_knot/table_knot_Alexander_polynomial.txt";

fn to_py_runtime_error(message: impl Into<String>) -> PyErr {
    PyRuntimeError::new_err(message.into())
}

fn resolve_table_path() -> Option<PathBuf> {
    if let Ok(path) = env::var("PYTHONKNOT_ALEXANDER_TABLE") {
        return Some(PathBuf::from(path));
    }

    let default_path = PathBuf::from(DEFAULT_TABLE_PATH);
    if default_path.exists() {
        return Some(default_path);
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
) -> PyResult<Vec<String>> {
    let results: Vec<Result<String, String>> = if let Some(threads) = num_threads {
        if threads == 0 {
            return Err(PyValueError::new_err("num_threads must be >= 1"));
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| to_py_runtime_error(e.to_string()))?;
        pool.install(|| {
            frames
                .par_iter()
                .map(|points| classify_frame(points, table, config))
                .collect()
        })
    } else {
        frames
            .iter()
            .map(|points| classify_frame(points, table, config))
            .collect()
    };

    let mut out = Vec::with_capacity(results.len());
    for result in results {
        match result {
            Ok(value) => out.push(value),
            Err(message) => return Err(to_py_runtime_error(message)),
        }
    }
    Ok(out)
}

fn compute_knot_size(
    frames: &[Vec<Point3>],
    table: &AlexanderTable,
    config: &KnotConfig,
) -> PyResult<(Vec<String>, Vec<Vec<i32>>)> {
    let mut knot_types = Vec::with_capacity(frames.len());
    let mut knot_sizes = Vec::with_capacity(frames.len());

    for points in frames {
        let knot_type = classify_frame(points, table, config).map_err(to_py_runtime_error)?;
        let size = if knot_type == "1" || !knot_type.contains('_') {
            vec![-1, -1, 0]
        } else {
            let core = find_knot_core(points, &knot_type, table, config)
                .map_err(|e| to_py_runtime_error(e.to_string()))?;
            if core.matched {
                vec![core.left, core.right, core.size]
            } else {
                vec![-1, -1, 0]
            }
        };
        knot_types.push(knot_type);
        knot_sizes.push(size);
    }

    Ok((knot_types, knot_sizes))
}

#[pyfunction]
fn read_xyz(py: Python<'_>, filename: &str) -> PyResult<Py<PyArray3<f64>>> {
    let frames = frames_from_path(Path::new(filename))?;
    Ok(frames_to_pyarray(py, &frames)?.to_owned())
}

#[pyfunction]
fn write_xyz(filename: &str, input_data: &PyAny) -> PyResult<()> {
    let frames = frames_from_input(input_data)?;
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(filename)
        .map_err(|e| to_py_runtime_error(e.to_string()))?;
    let mut writer = BufWriter::new(file);

    for frame in &frames {
        write_data_xyz(frame, &mut writer).map_err(|e| to_py_runtime_error(e.to_string()))?;
    }

    Ok(())
}

#[pyfunction(signature = (input_data, chain_type = "ring", num_threads = None))]
fn knot_type(
    input_data: &PyAny,
    chain_type: &str,
    num_threads: Option<usize>,
) -> PyResult<Vec<String>> {
    let frames = frames_from_input(input_data)?;
    let table = load_table()?;
    let config = config_from_chain_type(chain_type)?;
    classify_frames(&frames, &table, &config, num_threads)
}

#[pyfunction(signature = (input_data, chain_type = "ring"))]
fn knot_size(input_data: &PyAny, chain_type: &str) -> PyResult<(Vec<String>, Vec<Vec<i32>>)> {
    let frames = frames_from_input(input_data)?;
    let table = load_table()?;
    let config = config_from_chain_type(chain_type)?;
    compute_knot_size(&frames, &table, &config)
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
fn alexander_poly(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(write_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(knot_type, m)?)?;
    m.add_function(wrap_pyfunction!(knot_size, m)?)?;
    m.add_function(wrap_pyfunction!(kmt, m)?)?;

    let _ = py;
    Ok(())
}
