//! WebAssembly and TypeScript surface for PlanetaryHarmonics.
//!
//! This is the leaf of the composition chain. `ph-core` and `rustspice-core` are
//! consumed as **Rust libraries**, so there is exactly one WASM boundary — the one
//! JavaScript actually crosses.
//!
//! # Boundary discipline
//!
//! Boundary crossings dominate cost, so every API here is **columnar**: arrays in,
//! typed arrays out. There is no scalar per-call entry point, and none should be
//! added.
//!
//! That principle also decides what is exposed. [`ph_core::events`] takes a closure
//! over the angle, and calling a JavaScript closure per Newton iteration would put
//! the boundary inside the hot loop. So event finding is **not** exposed as a
//! generic root-finder: `aspectTimes` takes the bodies and the target and does the
//! whole search in Rust, returning only the answers.
//!
//! # Two families
//!
//! **Analytic** — constituent phases and commensurabilities need no kernels and
//! are pure arithmetic. Usable immediately.
//!
//! **Ephemeris-backed** — tidal tensors and aspect times need SPICE kernels, loaded
//! as bytes via [`Harmonics::load_kernel`].

use ph_core::{commensurability as com, doodson, events, fault, love, tidal::TidalTensor};
use rustspice_core::{Aberration, Et, KernelSet, Session};
use wasm_bindgen::prelude::*;

fn err(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// The library entry point.
#[wasm_bindgen]
pub struct Harmonics {
    kernels: KernelSet,
    session: Option<Session>,
}

impl Default for Harmonics {
    fn default() -> Self {
        Self::new()
    }
}

#[wasm_bindgen]
impl Harmonics {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Harmonics {
        Harmonics {
            kernels: KernelSet::new(),
            session: None,
        }
    }

    /// Add a SPICE kernel from bytes. Browsers have no filesystem, so this is the
    /// only way in. Invalidates any open session.
    #[wasm_bindgen(js_name = loadKernel)]
    pub fn load_kernel(&mut self, name: &str, data: &[u8]) {
        self.kernels.add(name, data.to_vec());
        self.session = None;
    }

    fn session(&mut self) -> Result<&mut Session, JsValue> {
        if self.session.is_none() {
            self.session = Some(self.kernels.open().map_err(err)?);
        }
        Ok(self.session.as_mut().unwrap())
    }

    // ---- analytic: no kernels required ----

    /// Names of the tidal constituents this build knows.
    #[wasm_bindgen(js_name = constituentNames)]
    pub fn constituent_names(&self) -> Vec<String> {
        doodson::CONSTITUENTS
            .iter()
            .map(|c| c.name.to_string())
            .collect()
    }

    /// Period of a named constituent, in days.
    #[wasm_bindgen(js_name = constituentPeriod)]
    pub fn constituent_period(&self, name: &str) -> Result<f64, JsValue> {
        doodson::constituent(name)
            .map(|c| c.period_days())
            .ok_or_else(|| err(format!("unknown constituent: {name}")))
    }

    /// Constituent phase at each time, radians in `[0, 2π)`.
    ///
    /// Times are days since 2000-01-01T00:00 UTC. `lonDeg` is east longitude —
    /// tidal phase is local, and for semidiurnal constituents a 180° error is a
    /// whole cycle.
    #[wasm_bindgen(js_name = constituentPhases)]
    pub fn constituent_phases(
        &self,
        name: &str,
        days: &[f64],
        lon_deg: f64,
    ) -> Result<Vec<f64>, JsValue> {
        let c = doodson::constituent(name)
            .ok_or_else(|| err(format!("unknown constituent: {name}")))?;
        Ok(days
            .iter()
            .map(|&d| c.phase_at_longitude(d, lon_deg))
            .collect())
    }

    /// Enumerate multi-body commensurabilities satisfying `Σ kᵢ = 0`.
    ///
    /// Returns the coefficients flattened row-major, `nBodies` per combination.
    #[wasm_bindgen(js_name = enumerateCommensurabilities)]
    pub fn enumerate_commensurabilities(&self, n_bodies: usize, max_coeff: i32) -> Vec<i32> {
        com::enumerate(n_bodies, max_coeff)
            .into_iter()
            .flat_map(|c| c.k)
            .collect()
    }

    /// Period of each commensurability, given mean motions in radians per day.
    ///
    /// `k` is flattened row-major with `rates.len()` entries per combination.
    /// Degenerate combinations, whose combined rate vanishes, yield `NaN`.
    #[wasm_bindgen(js_name = commensurabilityPeriods)]
    pub fn commensurability_periods(&self, k: &[i32], rates: &[f64]) -> Result<Vec<f64>, JsValue> {
        let n = rates.len();
        if n == 0 || k.len() % n != 0 {
            return Err(err("k length must be a multiple of rates length"));
        }
        Ok(k.chunks(n)
            .map(|row| {
                com::Commensurability::new(row.to_vec())
                    .and_then(|c| c.period(rates))
                    .unwrap_or(f64::NAN)
            })
            .collect())
    }

    /// Convert tidal tensor components (s⁻²) to stress (Pa) for Earth.
    ///
    /// Degree-2 scalar calibration via Love numbers; good to roughly a factor of 2.
    #[wasm_bindgen(js_name = stressFromTensor)]
    pub fn stress_from_tensor(&self, components: &[f64]) -> Vec<f64> {
        let e = love::Elastic::EARTH;
        components.iter().map(|&c| e.stress(c)).collect()
    }

    // ---- ephemeris-backed: kernels required ----

    /// Parse a time string to days since 2000-01-01T00:00 UTC.
    #[wasm_bindgen(js_name = parseTime)]
    pub fn parse_time(&mut self, s: &str) -> Result<f64, JsValue> {
        let sess = self.session()?;
        let epoch = sess.parse_time("2000-01-01T00:00:00").map_err(err)?;
        let t = sess.parse_time(s).map_err(err)?;
        Ok((t.0 - epoch.0) / 86400.0)
    }

    /// Times at which two bodies reach a given angular separation in ecliptic
    /// longitude — conjunctions, oppositions, and every aspect between.
    ///
    /// The entire search runs in Rust. Cost is proportional to the **number of
    /// events**, not to the precision requested, so millisecond resolution over
    /// decades is inexpensive.
    ///
    /// `meanPeriodDays` is the expected recurrence interval, seeding the linear
    /// predictor — the synodic period of the pair. Returns days since 2000-01-01.
    ///
    /// ⚠ Uses **apparent** positions (light-time and stellar aberration), the
    /// convention for astronomical event times. Tidal work wants geometric
    /// positions instead; the two differ by ~40 s for lunar phase.
    #[wasm_bindgen(js_name = aspectTimes)]
    #[allow(clippy::too_many_arguments)]
    pub fn aspect_times(
        &mut self,
        body_a: &str,
        body_b: &str,
        observer: &str,
        target_deg: f64,
        start_day: f64,
        end_day: f64,
        mean_period_days: f64,
        tol_seconds: f64,
    ) -> Result<Vec<f64>, JsValue> {
        if mean_period_days <= 0.0 || end_day <= start_day || tol_seconds <= 0.0 {
            return Err(err("invalid span, period, or tolerance"));
        }
        let epoch = {
            let s = self.session()?;
            s.parse_time("2000-01-01T00:00:00").map_err(err)?
        };
        let rate = std::f64::consts::TAU / mean_period_days;
        let target = target_deg.to_radians();

        // Interior mutability so the closure can borrow the session.
        let cell = std::cell::RefCell::new(self.session()?);
        let theta = |days: f64| -> f64 {
            let mut s = cell.borrow_mut();
            let et = Et(epoch.0 + days * 86400.0);
            let a = s
                .position(body_a, et, "ECLIPJ2000", observer, Aberration::LightTimeStellar)
                .map(|p| p.y.atan2(p.x))
                .unwrap_or(f64::NAN);
            let b = s
                .position(body_b, et, "ECLIPJ2000", observer, Aberration::LightTimeStellar)
                .map(|p| p.y.atan2(p.x))
                .unwrap_or(f64::NAN);
            let d = (a - b).rem_euclid(std::f64::consts::TAU);
            // Re-add the mean advance so the angle grows monotonically and the
            // linear predictor applies.
            d + std::f64::consts::TAU
                * (days / mean_period_days - d / std::f64::consts::TAU).round()
        };

        let tol = rate * (tol_seconds / 86400.0);
        Ok(events::crossings(theta, rate, target, start_day, end_day, tol)
            .into_iter()
            .map(|c| c.time)
            .collect())
    }

    /// Combined tidal tensor at each epoch, six components per epoch in
    /// `(xx, yy, zz, xy, xz, yz)` order, flattened.
    ///
    /// Geometric positions — tidal force acts on the instantaneous configuration.
    #[wasm_bindgen(js_name = tidalTensors)]
    pub fn tidal_tensors(
        &mut self,
        bodies: Vec<String>,
        days: &[f64],
        frame: &str,
        observer: &str,
    ) -> Result<Vec<f64>, JsValue> {
        let sess = self.session()?;
        let epoch = sess.parse_time("2000-01-01T00:00:00").map_err(err)?;
        let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch.0 + d * 86400.0)).collect();

        let mut acc = vec![TidalTensor::default(); epochs.len()];
        for name in &bodies {
            let gm = sess.constant(name, "GM").map_err(err)?[0];
            let ps = sess
                .positions(name, &epochs, frame, observer, Aberration::None)
                .map_err(err)?;
            for (a, p) in acc.iter_mut().zip(ps) {
                let t = TidalTensor::from_body(gm, [p.x, p.y, p.z]);
                for i in 0..3 {
                    for j in 0..3 {
                        a.m[i][j] += t.m[i][j];
                    }
                }
            }
        }
        Ok(acc.iter().flat_map(|t| fault::components(t)).collect())
    }
}
