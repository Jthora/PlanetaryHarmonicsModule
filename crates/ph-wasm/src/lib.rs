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

use ph_core::chart::{self, Chart, Frame};
use ph_core::{
    chart_cycles, chart_features, chart_local, commensurability as com, doodson, events, fault,
    love, tidal::TidalTensor,
};
use rustspice_core::{Aberration, Et, KernelSet, Session};
use wasm_bindgen::prelude::*;

fn err(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// Parse a frame name. Accepts the short forms the feature names use.
fn parse_frame(s: &str) -> Result<Frame, JsValue> {
    match s.to_ascii_lowercase().as_str() {
        "geo" | "geocentric" => Ok(Frame::Geocentric),
        "helio" | "heliocentric" => Ok(Frame::Heliocentric),
        "bary" | "barycentric" => Ok(Frame::Barycentric),
        other => Err(err(format!(
            "unknown frame '{other}' (expected geocentric, heliocentric or barycentric)"
        ))),
    }
}

/// How a feature vector is laid out. Shared by name and value production.
struct Spec {
    frames: Vec<Frame>,
    max_harmonic: usize,
    max_base: usize,
    cycle_harmonic: usize,
    site: Option<chart_local::Site>,
    site_harmonic: usize,
}

/// Build one epoch's feature vector from one chart per frame.
///
/// **This is the single source of column order.** `featureNames` runs it over
/// placeholder charts and `features` runs it over real ones, so the two cannot
/// disagree about what column 4,912 is — a guarantee of construction rather than
/// of discipline. Getting that wrong would mislabel every column silently, which
/// is the worst kind of bug this API could have.
fn assemble(charts: &[Chart], spec: &Spec) -> chart_features::FeatureSet {
    let mut out = chart_features::FeatureSet::default();
    for (frame, c) in spec.frames.iter().zip(charts) {
        out.extend(chart_features::all(c, spec.max_harmonic, spec.max_base));
        // Lunar and eclipse features are geocentric quantities by nature; in other
        // frames the Sun or Earth is the observer and those families come back
        // empty of their own accord, which keeps this branch-free.
        let cyc = chart_cycles::all(c, spec.cycle_harmonic);
        let prefix = match frame {
            Frame::Geocentric => "geo",
            Frame::Heliocentric => "helio",
            Frame::Barycentric => "bary",
        };
        for (n, v) in cyc.names.iter().zip(&cyc.values) {
            out.push(format!("{prefix}.{n}"), *v);
        }
    }
    if let Some(site) = spec.site {
        // Site angles need the apparent sky, so they come from the geocentric
        // chart. `features` rejects a spec that asks for them without one.
        if let Some(i) = spec.frames.iter().position(|f| *f == Frame::Geocentric) {
            out.extend(chart_local::all(&charts[i], site, spec.site_harmonic));
        }
    }
    out
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

    #[allow(clippy::too_many_arguments)]
    fn spec(
        frames: Vec<String>,
        max_harmonic: usize,
        max_base: usize,
        cycle_harmonic: usize,
        with_site: bool,
        site_lat_deg: f64,
        site_lon_deg: f64,
        site_harmonic: usize,
    ) -> Result<Spec, JsValue> {
        if frames.is_empty() {
            return Err(err("at least one frame is required"));
        }
        let fs: Result<Vec<Frame>, JsValue> = frames.iter().map(|f| parse_frame(f)).collect();
        let fs = fs?;
        if with_site && !fs.contains(&Frame::Geocentric) {
            return Err(err(
                "site-local angles need the geocentric frame; add 'geocentric' to frames",
            ));
        }
        if !(-90.0..=90.0).contains(&site_lat_deg) && with_site {
            return Err(err("site latitude must be between -90 and 90"));
        }
        Ok(Spec {
            frames: fs,
            max_harmonic,
            max_base,
            cycle_harmonic,
            site: if with_site {
                Some(chart_local::Site::from_degrees(site_lat_deg, site_lon_deg))
            } else {
                None
            },
            site_harmonic,
        })
    }

    fn session(&mut self) -> Result<&mut Session, JsValue> {
        if self.session.is_none() {
            self.session = Some(self.kernels.open().map_err(err)?);
        }
        Ok(self.session.as_mut().unwrap())
    }

    // ---- charts and derived features ----

    /// Bodies carried in every chart, in the order [`charts`] returns them.
    #[wasm_bindgen(js_name = bodyNames)]
    pub fn body_names(&self) -> Vec<String> {
        chart::BODIES.iter().map(|s| s.to_string()).collect()
    }

    /// Chart primitives: eight numbers per body per epoch, flattened.
    ///
    /// Per body, in order: ecliptic longitude and latitude (radians), distance
    /// (km), their three time derivatives (per day), then right ascension and
    /// declination (radians). Row-major by epoch, then by body.
    ///
    /// **Geometric positions, no aberration** — these are dynamical quantities.
    /// [`aspect_times`] uses apparent positions instead, which is the right
    /// convention for event timing and the wrong one for feature generation; the
    /// two differ by about 40 seconds of lunar phase.
    #[wasm_bindgen(js_name = charts)]
    pub fn charts(&mut self, days: &[f64], frame: &str) -> Result<Vec<f64>, JsValue> {
        let f = parse_frame(frame)?;
        let sess = self.session()?;
        let epoch = sess.parse_time("2000-01-01T00:00:00").map_err(err)?;
        let cs = chart::charts(sess, days, f, epoch).map_err(err)?;
        let mut out = Vec::with_capacity(days.len() * chart::BODIES.len() * 8);
        for c in &cs {
            for s in &c.states {
                out.extend_from_slice(&[
                    s.lon, s.lat, s.dist, s.lon_speed, s.lat_speed, s.dist_speed, s.ra, s.dec,
                ]);
            }
        }
        Ok(out)
    }

    // ---- analytic: no kernels required ----

    /// Column names for a feature spec, in the exact order [`features`] returns.
    ///
    /// Needs no kernels, so a caller can fetch and cache the names before any
    /// ephemeris is loaded. **The order is an API contract**: cache it once and
    /// index into the float32 matrix by position.
    #[wasm_bindgen(js_name = featureNames)]
    #[allow(clippy::too_many_arguments)]
    pub fn feature_names(
        &self,
        frames: Vec<String>,
        max_harmonic: usize,
        max_base: usize,
        cycle_harmonic: usize,
        with_site: bool,
        site_lat_deg: f64,
        site_lon_deg: f64,
        site_harmonic: usize,
    ) -> Result<Vec<String>, JsValue> {
        let spec = Self::spec(
            frames,
            max_harmonic,
            max_base,
            cycle_harmonic,
            with_site,
            site_lat_deg,
            site_lon_deg,
            site_harmonic,
        )?;
        let placeholders: Vec<Chart> = spec.frames.iter().map(|f| Chart::placeholder(*f)).collect();
        Ok(assemble(&placeholders, &spec).names)
    }

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

    /// Derived feature matrix: one row per epoch, `featureNames().length` columns,
    /// row-major, **float32**.
    ///
    /// f32 because every feature is a trigonometric quantity or a circular
    /// statistic whose seventh digit is far below anything a consumer can use, and
    /// halving the matrix decides whether a long span fits in memory.
    ///
    /// Pass `withSite = false` and any latitude to skip site-local angles. With it
    /// enabled, `frames` **must** include geocentric — site angles are properties
    /// of the apparent sky, and computing them from a heliocentric chart would
    /// return numbers that mean nothing rather than an error.
    #[wasm_bindgen(js_name = features)]
    #[allow(clippy::too_many_arguments)]
    pub fn features(
        &mut self,
        days: &[f64],
        frames: Vec<String>,
        max_harmonic: usize,
        max_base: usize,
        cycle_harmonic: usize,
        with_site: bool,
        site_lat_deg: f64,
        site_lon_deg: f64,
        site_harmonic: usize,
    ) -> Result<Vec<f32>, JsValue> {
        let spec = Self::spec(
            frames,
            max_harmonic,
            max_base,
            cycle_harmonic,
            with_site,
            site_lat_deg,
            site_lon_deg,
            site_harmonic,
        )?;
        let epoch = {
            let s = self.session()?;
            s.parse_time("2000-01-01T00:00:00").map_err(err)?
        };
        let mut per_frame = Vec::with_capacity(spec.frames.len());
        for f in &spec.frames {
            let sess = self.session()?;
            per_frame.push(chart::charts(sess, days, *f, epoch).map_err(err)?);
        }

        let mut out: Vec<f32> = Vec::new();
        let mut width = 0usize;
        for i in 0..days.len() {
            let row: Vec<Chart> = per_frame.iter().map(|cs| cs[i].clone()).collect();
            let fs = assemble(&row, &spec);
            if width == 0 {
                width = fs.len();
                out.reserve(width * days.len());
            }
            out.extend(fs.values.iter().map(|v| *v as f32));
        }
        Ok(out)
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
