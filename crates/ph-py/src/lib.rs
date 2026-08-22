//! Python bindings for the analysis layer.
//!
//! # What is and is not exposed
//!
//! **Exposed:** analytic constituent phases, elastic stress conversion, the
//! Schuster statistic, and the block-shift null. These are pure computation with
//! no SPICE state, so they bind cleanly and are what a Python forecasting or ML
//! stack actually calls in a loop.
//!
//! **Not exposed:** ephemeris and tidal-tensor computation. Those need a SPICE
//! session, which holds `Rc`-based state and is not `Send`. They are also batch
//! work rather than interactive — use the `ph-features` CLI, which emits CSV with
//! a provenance header (`docs/14` §5).
//!
//! The split is deliberate: heavy geometry through the CLI, analysis primitives
//! through here.
//!
//! # Build
//!
//! ```text
//! cd crates/ph-py && maturin develop --release
//! ```

use ph_core::{doodson, love, stats};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Names of the constituents this module knows.
#[pyfunction]
fn constituents() -> Vec<String> {
    doodson::CONSTITUENTS.iter().map(|c| c.name.to_string()).collect()
}

/// Period of a named constituent, in days.
#[pyfunction]
fn constituent_period(name: &str) -> PyResult<f64> {
    doodson::constituent(name)
        .map(|c| c.period_days())
        .ok_or_else(|| PyValueError::new_err(format!("unknown constituent: {name}")))
}

/// Analytic phase of a constituent at each time, radians in `[0, 2π)`.
///
/// Times are **days since 2000-01-01T00:00 UTC**. The phase is an integer
/// combination of the six Doodson fundamental arguments, so it is uniform over
/// long spans by construction — the property a Schuster null requires.
#[pyfunction]
fn constituent_phases(name: &str, days: Vec<f64>) -> PyResult<Vec<f64>> {
    let c = doodson::constituent(name)
        .ok_or_else(|| PyValueError::new_err(format!("unknown constituent: {name}")))?;
    Ok(c.phases(&days))
}

/// Convert a tidal tensor component (s⁻²) to stress (Pa) for Earth.
///
/// Degree-2 scalar calibration via Love numbers; good to roughly a factor of 2.
/// Not a per-component elastic solution.
#[pyfunction]
fn stress_from_tensor(component: f64) -> f64 {
    love::Elastic::EARTH.stress(component)
}

/// Strain corresponding to a tidal tensor component, dimensionless.
#[pyfunction]
fn strain_from_tensor(component: f64) -> f64 {
    love::Elastic::EARTH.strain(component)
}

/// Critical period `T_a = 2π Aσ₀ / τ̇`, in the units of `stressing_rate`.
#[pyfunction]
fn critical_period(a_sigma: f64, stressing_rate: f64) -> f64 {
    love::critical_period(a_sigma, stressing_rate)
}

/// Schuster statistics for orders `1..=max_order`.
///
/// Returns a list of `(order, d_squared, p_value, phase)`.
///
/// ⚠ The analytic `p_value` assumes independent, uniformly sampled phases.
/// Earthquake catalogues violate both. Use it to screen, and establish
/// significance with [`block_shift_null`].
#[pyfunction]
#[pyo3(signature = (phases, max_order = 1))]
fn schuster(phases: Vec<f64>, max_order: usize) -> PyResult<Vec<(usize, f64, f64, f64)>> {
    if max_order == 0 {
        return Err(PyValueError::new_err("max_order must be at least 1"));
    }
    if phases.is_empty() {
        return Err(PyValueError::new_err("phases must not be empty"));
    }
    Ok(stats::schuster(&phases, max_order)
        .into_iter()
        .map(|s| (s.order, s.d_squared, s.p_value, s.phase))
        .collect())
}

/// Normalised Schuster power `D²/N` for a constituent at the given event times.
#[pyfunction]
fn constituent_power(name: &str, days: Vec<f64>) -> PyResult<f64> {
    let c = doodson::constituent(name)
        .ok_or_else(|| PyValueError::new_err(format!("unknown constituent: {name}")))?;
    if days.is_empty() {
        return Err(PyValueError::new_err("days must not be empty"));
    }
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &d in &days {
        let (s, co) = c.phase_at(d).sin_cos();
        a += co;
        b += s;
    }
    Ok((a * a + b * b) / days.len() as f64)
}

/// Block-shift null for a constituent: observed power and the null distribution.
///
/// Returns `(observed, sorted_null_samples)`.
///
/// A **global** time shift cannot work for a single constituent — `D²` is
/// invariant under rotation and a global shift *is* a rotation. Each block is
/// therefore shifted independently, preserving within-block clustering while
/// randomising alignment between blocks. Block length defaults to
/// `max(4 × period, 30 d)`.
#[pyfunction]
#[pyo3(signature = (name, days, trials = 400, seed = 0x0D00D5, block_days = None))]
fn block_shift_null(
    name: &str,
    days: Vec<f64>,
    trials: usize,
    seed: u64,
    block_days: Option<f64>,
) -> PyResult<(f64, Vec<f64>)> {
    let c = doodson::constituent(name)
        .ok_or_else(|| PyValueError::new_err(format!("unknown constituent: {name}")))?;
    if days.len() < 2 {
        return Err(PyValueError::new_err("need at least two event times"));
    }
    let period = c.period_days();
    let block = block_days.unwrap_or_else(|| (4.0 * period).max(30.0));
    if block <= 0.0 {
        return Err(PyValueError::new_err("block_days must be positive"));
    }

    let t0 = days.iter().cloned().fold(f64::MAX, f64::min);
    let t1 = days.iter().cloned().fold(f64::MIN, f64::max);
    let n_blocks = (((t1 - t0) / block).floor() as usize) + 2;

    let power = |ph: &dyn Fn(f64) -> f64| -> f64 {
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for &d in &days {
            let (s, co) = ph(d).sin_cos();
            a += co;
            b += s;
        }
        (a * a + b * b) / days.len() as f64
    };

    let observed = power(&|d| c.phase_at(d));

    let mut state = seed | 1;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / ((1u64 << 53) as f64)
    };

    let mut null = Vec::with_capacity(trials);
    for _ in 0..trials {
        let offs: Vec<f64> = (0..n_blocks).map(|_| next() * period).collect();
        null.push(power(&|d| {
            let b = (((d - t0) / block).floor() as usize).min(offs.len() - 1);
            c.phase_at(d + offs[b])
        }));
    }
    null.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Ok((observed, null))
}

#[pymodule]
fn planetary_harmonics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", "Analysis primitives for celestial-terrestrial harmonic feature analysis.")?;
    m.add_function(wrap_pyfunction!(constituents, m)?)?;
    m.add_function(wrap_pyfunction!(constituent_period, m)?)?;
    m.add_function(wrap_pyfunction!(constituent_phases, m)?)?;
    m.add_function(wrap_pyfunction!(constituent_power, m)?)?;
    m.add_function(wrap_pyfunction!(stress_from_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(strain_from_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(critical_period, m)?)?;
    m.add_function(wrap_pyfunction!(schuster, m)?)?;
    m.add_function(wrap_pyfunction!(block_shift_null, m)?)?;
    Ok(())
}
