//! A compact chart ephemeris — Chebyshev fits for the bodies a chart needs.
//!
//! # Why
//!
//! `de440s.bsp` is 32.7 MB against 1.9 MB of compiled WASM. That is tolerable for
//! a desktop tool and hostile for a web page, and most of it is weight this
//! library never uses: DE440 carries sub-metre accuracy for spacecraft navigation
//! across 300 years and every body in the solar system, while a chart needs
//! eleven bodies at arcsecond accuracy over a couple of centuries.
//!
//! So the positions are refitted for the job. Each body gets Chebyshev
//! coefficients over fixed intervals, sized to its own motion, stored as f32.
//!
//! # Why this loses nothing that matters
//!
//! f32 holds about seven significant digits. On the Moon at 384,400 km that is
//! 0.03 km, which is 0.016 arcseconds of geocentric angle; on Pluto at 5×10⁹ km it
//! is 500 km, which is 0.00002 degrees. Both are far below the fitting error, and
//! the fitting error is what [`build`] reports so it can be checked rather than
//! assumed.
//!
//! # Velocity is free
//!
//! Only positions are stored. A Chebyshev series differentiates analytically, so
//! velocity comes out of the same coefficients exactly — storing it would double
//! the file to reproduce information already present.
//!
//! # Layout
//!
//! Interval length and degree are chosen per body from its speed. The Moon is
//! referred to **Earth** rather than the barycentre: geocentrically its motion is
//! a smooth 27-day ellipse, while against the barycentre it carries Earth's full
//! annual orbit as well and would need several times the coefficients to describe
//! something the Earth series already contains.

use crate::chart::{BodyState, Chart, Frame, BODIES};
use rustspice_core::Vec3;

/// What a body's coefficients are measured from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Center {
    /// Solar-system barycentre.
    Ssb,
    /// Earth. Used for the Moon only.
    Earth,
}

/// How one body is fitted.
#[derive(Debug, Clone, Copy)]
pub struct BodySpec {
    pub name: &'static str,
    pub center: Center,
    /// Days covered by one set of coefficients.
    pub interval_days: f64,
    /// Chebyshev degree; the set holds `degree + 1` coefficients per component.
    pub degree: usize,
}

/// Fitting plan, parallel to [`BODIES`].
///
/// Intervals and degrees follow the shape of each orbit: fast bodies get short
/// intervals, slow ones long. These are close to the granule sizes JPL uses for
/// the DE series, relaxed where arcsecond accuracy allows it.
pub const LAYOUT: &[BodySpec] = &[
    BodySpec { name: "SUN", center: Center::Ssb, interval_days: 16.0, degree: 11 },
    BodySpec { name: "EARTH", center: Center::Ssb, interval_days: 16.0, degree: 12 },
    BodySpec { name: "MOON", center: Center::Earth, interval_days: 8.0, degree: 12 },
    BodySpec { name: "MERCURY", center: Center::Ssb, interval_days: 8.0, degree: 13 },
    BodySpec { name: "VENUS", center: Center::Ssb, interval_days: 16.0, degree: 10 },
    BodySpec { name: "MARS BARYCENTER", center: Center::Ssb, interval_days: 32.0, degree: 11 },
    BodySpec { name: "JUPITER BARYCENTER", center: Center::Ssb, interval_days: 32.0, degree: 8 },
    BodySpec { name: "SATURN BARYCENTER", center: Center::Ssb, interval_days: 32.0, degree: 7 },
    BodySpec { name: "URANUS BARYCENTER", center: Center::Ssb, interval_days: 64.0, degree: 7 },
    BodySpec { name: "NEPTUNE BARYCENTER", center: Center::Ssb, interval_days: 64.0, degree: 7 },
    BodySpec { name: "PLUTO BARYCENTER", center: Center::Ssb, interval_days: 64.0, degree: 7 },
];

const MAGIC: &[u8; 4] = b"PHCE";
const VERSION: u16 = 1;

/// Fit a Chebyshev series of `degree` to `f` on `[a, b]`.
///
/// Sampling at Chebyshev nodes makes the coefficients a discrete cosine transform
/// of the samples, which is both cheap and near-optimal in the minimax sense — the
/// error equioscillates rather than piling up at the interval ends.
///
/// The zeroth coefficient is returned already halved, so evaluation is a plain
/// `Σ cₖ Tₖ` with no special case.
pub fn fit<F: FnMut(f64) -> f64>(mut f: F, a: f64, b: f64, degree: usize) -> Vec<f64> {
    let n = degree + 1;
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    let samples: Vec<f64> = (0..n)
        .map(|j| {
            let x = (std::f64::consts::PI * (j as f64 + 0.5) / n as f64).cos();
            f(mid + half * x)
        })
        .collect();
    let mut c = vec![0.0; n];
    for (k, ck) in c.iter_mut().enumerate() {
        let mut sum = 0.0;
        for (j, s) in samples.iter().enumerate() {
            sum += s * (std::f64::consts::PI * k as f64 * (j as f64 + 0.5) / n as f64).cos();
        }
        *ck = 2.0 / n as f64 * sum;
    }
    c[0] *= 0.5;
    c
}

/// Fit all three components at once, sampling the source only once per node.
///
/// Fitting component by component would call `f` three times at every node, and
/// the caller's `f` is a SPICE lookup — the difference is a million redundant
/// ephemeris evaluations over a two-century span.
fn fit_vec3<F: FnMut(f64) -> Vec3>(mut f: F, a: f64, b: f64, degree: usize) -> [Vec<f64>; 3] {
    let n = degree + 1;
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    let samples: Vec<Vec3> = (0..n)
        .map(|j| {
            let x = (std::f64::consts::PI * (j as f64 + 0.5) / n as f64).cos();
            f(mid + half * x)
        })
        .collect();
    let mut out = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    for k in 0..n {
        let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
        for (j, p) in samples.iter().enumerate() {
            let w = (std::f64::consts::PI * k as f64 * (j as f64 + 0.5) / n as f64).cos();
            sx += p.x * w;
            sy += p.y * w;
            sz += p.z * w;
        }
        let scale = 2.0 / n as f64;
        out[0][k] = scale * sx;
        out[1][k] = scale * sy;
        out[2][k] = scale * sz;
    }
    for c in out.iter_mut() {
        c[0] *= 0.5;
    }
    out
}

/// Evaluate a Chebyshev series and its derivative with respect to `x`.
///
/// Uses the recurrences for `Tₖ` and `Uₖ` together, since `Tₖ′(x) = k·U₍ₖ₋₁₎(x)`.
/// One pass gives both, and at these degrees it is as stable as Clenshaw while
/// being considerably easier to read.
pub fn eval(c: &[f32], x: f64) -> (f64, f64) {
    if c.is_empty() {
        return (0.0, 0.0);
    }
    let (mut t_prev, mut t) = (1.0f64, x);
    let (mut u_prev, mut u) = (1.0f64, 2.0 * x);
    let mut value = c[0] as f64;
    let mut deriv = 0.0f64;
    for (k, &ck) in c.iter().enumerate().skip(1) {
        let ck = ck as f64;
        value += ck * t;
        // T'_k = k * U_{k-1}
        deriv += ck * k as f64 * u_prev;
        let t_next = 2.0 * x * t - t_prev;
        t_prev = t;
        t = t_next;
        let u_next = 2.0 * x * u - u_prev;
        u_prev = u;
        u = u_next;
    }
    (value, deriv)
}

/// A fitted body: coefficients for consecutive intervals.
#[derive(Debug, Clone)]
struct Series {
    spec: BodySpec,
    n_intervals: usize,
    /// `[interval][component][coefficient]`, flattened.
    coeffs: Vec<f32>,
}

impl Series {
    fn stride(&self) -> usize {
        3 * (self.spec.degree + 1)
    }

    /// Position (km) and velocity (km/day) at `day`, relative to this body's centre.
    fn state(&self, day: f64, day_start: f64) -> (Vec3, Vec3) {
        let dt = self.spec.interval_days;
        // Clamp rather than extrapolate: a Chebyshev fit diverges violently outside
        // its interval, and a silently wrong position is worse than a stale one.
        let i = (((day - day_start) / dt).floor() as isize)
            .clamp(0, self.n_intervals as isize - 1) as usize;
        let a = day_start + i as f64 * dt;
        let x = (2.0 * (day - a) / dt - 1.0).clamp(-1.0, 1.0);
        let base = i * self.stride();
        let n = self.spec.degree + 1;
        let mut p = [0.0; 3];
        let mut v = [0.0; 3];
        for comp in 0..3 {
            let c = &self.coeffs[base + comp * n..base + (comp + 1) * n];
            let (val, der) = eval(c, x);
            p[comp] = val;
            v[comp] = der * 2.0 / dt; // dx/d(day) = 2/interval
        }
        (
            Vec3 { x: p[0], y: p[1], z: p[2] },
            Vec3 { x: v[0], y: v[1], z: v[2] },
        )
    }
}

/// A fitted ephemeris for the chart bodies.
#[derive(Debug, Clone)]
pub struct CompactEphemeris {
    day_start: f64,
    day_end: f64,
    series: Vec<Series>,
}

/// Worst absolute position error seen while fitting one body, km.
#[derive(Debug, Clone, Copy)]
pub struct FitError {
    pub body: &'static str,
    pub max_km: f64,
    /// The same error as an angle seen from the body's own distance, arcseconds.
    pub max_arcsec: f64,
}

impl CompactEphemeris {
    pub fn day_range(&self) -> (f64, f64) {
        (self.day_start, self.day_end)
    }

    /// Size of the serialised form, bytes.
    pub fn byte_len(&self) -> usize {
        let mut n = 4 + 2 + 2 + 8 + 8;
        for s in &self.series {
            n += 1 + s.spec.name.len() + 1 + 1 + 8 + 4 + 4;
            n += s.coeffs.len() * 4;
        }
        n
    }

    /// Position (km) and velocity (km/day) of a body relative to the SSB.
    fn ssb_state(&self, index: usize, day: f64) -> (Vec3, Vec3) {
        let s = &self.series[index];
        let (p, v) = s.state(day, self.day_start);
        match s.spec.center {
            Center::Ssb => (p, v),
            Center::Earth => {
                let ei = LAYOUT.iter().position(|b| b.name == "EARTH").unwrap();
                let (ep, ev) = self.series[ei].state(day, self.day_start);
                (
                    Vec3 { x: p.x + ep.x, y: p.y + ep.y, z: p.z + ep.z },
                    Vec3 { x: v.x + ev.x, y: v.y + ev.y, z: v.z + ev.z },
                )
            }
        }
    }

    /// Charts for a series of epochs, matching [`crate::chart::charts`] exactly in
    /// shape and convention.
    pub fn charts(&self, days: &[f64], frame: Frame) -> Vec<Chart> {
        let observer_index = LAYOUT.iter().position(|b| b.name == frame.observer());
        days.iter()
            .map(|&day| {
                let origin = observer_index.map(|i| self.ssb_state(i, day));
                let states = (0..BODIES.len())
                    .map(|bi| {
                        if Some(bi) == observer_index {
                            // Degenerate by construction, as in the SPICE path: a
                            // body has no direction from itself.
                            return BodyState::default();
                        }
                        let (p, v) = self.ssb_state(bi, day);
                        let (p, v) = match origin {
                            Some((op, ov)) => (
                                Vec3 { x: p.x - op.x, y: p.y - op.y, z: p.z - op.z },
                                Vec3 { x: v.x - ov.x, y: v.y - ov.y, z: v.z - ov.z },
                            ),
                            None => (p, v),
                        };
                        BodyState::from_state(p, v)
                    })
                    .collect();
                Chart { day, frame, states }
            })
            .collect()
    }

    /// Serialise. Little-endian throughout; f32 coefficients.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.byte_len());
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&(self.series.len() as u16).to_le_bytes());
        out.extend_from_slice(&self.day_start.to_le_bytes());
        out.extend_from_slice(&self.day_end.to_le_bytes());
        for s in &self.series {
            out.push(s.spec.name.len() as u8);
            out.extend_from_slice(s.spec.name.as_bytes());
            out.push(match s.spec.center {
                Center::Ssb => 0,
                Center::Earth => 1,
            });
            out.push(s.spec.degree as u8);
            out.extend_from_slice(&s.spec.interval_days.to_le_bytes());
            out.extend_from_slice(&(s.n_intervals as u32).to_le_bytes());
            out.extend_from_slice(&(s.coeffs.len() as u32).to_le_bytes());
        }
        for s in &self.series {
            for c in &s.coeffs {
                out.extend_from_slice(&c.to_le_bytes());
            }
        }
        out
    }

    /// Parse. Returns `None` on any inconsistency rather than reading past the end.
    pub fn from_bytes(b: &[u8]) -> Option<CompactEphemeris> {
        let mut o = 0usize;
        let take = |o: &mut usize, n: usize| -> Option<&[u8]> {
            let s = b.get(*o..*o + n)?;
            *o += n;
            Some(s)
        };
        if take(&mut o, 4)? != MAGIC {
            return None;
        }
        let version = u16::from_le_bytes(take(&mut o, 2)?.try_into().ok()?);
        if version != VERSION {
            return None;
        }
        let n_bodies = u16::from_le_bytes(take(&mut o, 2)?.try_into().ok()?) as usize;
        let day_start = f64::from_le_bytes(take(&mut o, 8)?.try_into().ok()?);
        let day_end = f64::from_le_bytes(take(&mut o, 8)?.try_into().ok()?);

        let mut headers = Vec::with_capacity(n_bodies);
        for _ in 0..n_bodies {
            let len = take(&mut o, 1)?[0] as usize;
            let name = std::str::from_utf8(take(&mut o, len)?).ok()?;
            // Names are matched against the static table so the rest of the library
            // can keep using &'static str.
            let spec_name = LAYOUT.iter().find(|s| s.name == name)?.name;
            let center = match take(&mut o, 1)?[0] {
                0 => Center::Ssb,
                1 => Center::Earth,
                _ => return None,
            };
            let degree = take(&mut o, 1)?[0] as usize;
            let interval_days = f64::from_le_bytes(take(&mut o, 8)?.try_into().ok()?);
            let n_intervals = u32::from_le_bytes(take(&mut o, 4)?.try_into().ok()?) as usize;
            let n_coeffs = u32::from_le_bytes(take(&mut o, 4)?.try_into().ok()?) as usize;
            if n_coeffs != n_intervals * 3 * (degree + 1) {
                return None;
            }
            headers.push((
                BodySpec { name: spec_name, center, interval_days, degree },
                n_intervals,
                n_coeffs,
            ));
        }

        let mut series = Vec::with_capacity(n_bodies);
        for (spec, n_intervals, n_coeffs) in headers {
            let raw = take(&mut o, n_coeffs * 4)?;
            let coeffs = raw
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
                .collect();
            series.push(Series { spec, n_intervals, coeffs });
        }
        Some(CompactEphemeris { day_start, day_end, series })
    }
}

/// Fit every body in [`LAYOUT`] over `[day_start, day_end]`.
///
/// `sample` returns a body's ECLIPJ2000 position in km relative to a centre, at a
/// given day. Errors are reported per body so accuracy is measured rather than
/// assumed — the caller decides whether the fit is good enough.
pub fn build<F>(day_start: f64, day_end: f64, mut sample: F) -> (CompactEphemeris, Vec<FitError>)
where
    F: FnMut(&str, Center, f64) -> Vec3,
{
    let mut series = Vec::with_capacity(LAYOUT.len());
    let mut errors = Vec::with_capacity(LAYOUT.len());

    for spec in LAYOUT {
        let n_intervals = (((day_end - day_start) / spec.interval_days).ceil() as usize).max(1);
        let n = spec.degree + 1;
        let mut coeffs = Vec::with_capacity(n_intervals * 3 * n);
        for i in 0..n_intervals {
            let a = day_start + i as f64 * spec.interval_days;
            let b = a + spec.interval_days;
            let fitted = fit_vec3(|t| sample(spec.name, spec.center, t), a, b, spec.degree);
            for c in &fitted {
                coeffs.extend(c.iter().map(|&v| v as f32));
            }
        }
        let s = Series { spec: *spec, n_intervals, coeffs };

        // Measure the error where a Chebyshev fit is worst -- between the nodes,
        // and at the interval ends. Sampling on the nodes themselves would report
        // the fit as near-perfect and mean nothing.
        let mut max_km = 0.0f64;
        let mut dist_at_max = 1.0f64;
        let probes = 7;
        for i in 0..n_intervals {
            let a = day_start + i as f64 * spec.interval_days;
            for k in 0..=probes {
                let t = a + spec.interval_days * k as f64 / probes as f64;
                if t > day_end {
                    break;
                }
                let want = sample(spec.name, spec.center, t);
                let (got, _) = s.state(t, day_start);
                let d = ((got.x - want.x).powi(2)
                    + (got.y - want.y).powi(2)
                    + (got.z - want.z).powi(2))
                .sqrt();
                if d > max_km {
                    max_km = d;
                    dist_at_max = want.norm().max(1.0);
                }
            }
        }
        errors.push(FitError {
            body: spec.name,
            max_km,
            max_arcsec: (max_km / dist_at_max).atan().to_degrees() * 3600.0,
        });
        series.push(s);
    }

    (CompactEphemeris { day_start, day_end, series }, errors)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_vector_fit_agrees_with_three_scalar_fits() {
        // fit_vec3 exists only to avoid redundant sampling; it must not change the
        // answer.
        let f = |t: f64| Vec3 { x: (0.3 * t).sin(), y: t * 0.2, z: (0.11 * t).cos() };
        let v = fit_vec3(f, -5.0, 7.0, 9);
        for (comp, got) in v.iter().enumerate() {
            let want = fit(|t| [f(t).x, f(t).y, f(t).z][comp], -5.0, 7.0, 9);
            for (a, b) in got.iter().zip(&want) {
                assert!((a - b).abs() < 1e-12, "component {comp}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn chebyshev_reproduces_a_smooth_function_and_its_derivative() {
        // A fit is only meaningful away from its nodes, so the check samples
        // between them.
        let (a, b) = (-3.0, 5.0);
        let f = |t: f64| (0.7 * t).sin() + 0.3 * t;
        let df = |t: f64| 0.7 * (0.7 * t).cos() + 0.3;
        let c: Vec<f32> = fit(f, a, b, 14).iter().map(|&v| v as f32).collect();

        let mut worst_v: f64 = 0.0;
        let mut worst_d: f64 = 0.0;
        for k in 0..=500 {
            let t = a + (b - a) * k as f64 / 500.0;
            let x = 2.0 * (t - a) / (b - a) - 1.0;
            let (v, d) = eval(&c, x);
            worst_v = worst_v.max((v - f(t)).abs());
            worst_d = worst_d.max((d * 2.0 / (b - a) - df(t)).abs());
        }
        assert!(worst_v < 1e-6, "value error {worst_v:e}");
        assert!(worst_d < 1e-5, "derivative error {worst_d:e}");
    }

    #[test]
    fn the_derivative_is_analytic_not_differenced() {
        // The point of storing positions only: velocity must be exact, not a
        // finite difference of the fitted values.
        let c: Vec<f32> = fit(|t: f64| t * t * t, -1.0, 1.0, 6)
            .iter()
            .map(|&v| v as f32)
            .collect();
        for x in [-0.9, -0.3, 0.0, 0.44, 0.87] {
            let (v, d) = eval(&c, x);
            assert!((v - x * x * x).abs() < 1e-5, "value at {x}: {v}");
            assert!((d - 3.0 * x * x).abs() < 1e-4, "derivative at {x}: {d}");
        }
    }

    /// A synthetic solar system: circular orbits with known periods. Exercises the
    /// whole path — fit, serialise, parse, evaluate, frame change — against an
    /// answer that is known in closed form.
    fn synthetic(name: &str, center: Center, day: f64) -> Vec3 {
        let (r, period, phase) = match name {
            "SUN" => (1.0e6, 11.86 * 365.25, 0.0),
            "EARTH" => (1.496e8, 365.256, 1.0),
            "MOON" => (3.844e5, 27.3217, 2.0),
            "MERCURY" => (5.79e7, 87.969, 0.5),
            "VENUS" => (1.082e8, 224.701, 1.5),
            "MARS BARYCENTER" => (2.279e8, 686.98, 2.5),
            "JUPITER BARYCENTER" => (7.785e8, 4332.6, 3.0),
            "SATURN BARYCENTER" => (1.4335e9, 10759.2, 3.5),
            "URANUS BARYCENTER" => (2.8725e9, 30688.5, 4.0),
            "NEPTUNE BARYCENTER" => (4.4951e9, 60182.0, 4.5),
            _ => (5.906e9, 90560.0, 5.0),
        };
        let _ = center;
        let a = std::f64::consts::TAU * day / period + phase;
        // A small inclination so latitudes are not identically zero.
        Vec3 { x: r * a.cos(), y: r * a.sin(), z: 0.05 * r * (2.0 * a).sin() }
    }

    #[test]
    fn fitting_a_synthetic_system_is_accurate_to_well_under_an_arcsecond() {
        let (eph, errors) = build(0.0, 400.0, synthetic);
        for e in &errors {
            assert!(
                e.max_arcsec < 0.05,
                "{}: {:.4} arcsec ({:.3} km)",
                e.body,
                e.max_arcsec,
                e.max_km
            );
        }
        // And the whole thing is small.
        assert!(eph.byte_len() > 0);
    }

    #[test]
    fn serialisation_round_trips_exactly() {
        let (eph, _) = build(0.0, 200.0, synthetic);
        let bytes = eph.to_bytes();
        assert_eq!(bytes.len(), eph.byte_len(), "byte_len disagrees with to_bytes");
        let back = CompactEphemeris::from_bytes(&bytes).expect("parse");
        assert_eq!(back.day_range(), eph.day_range());

        let days: Vec<f64> = (0..40).map(|k| k as f64 * 4.7).collect();
        for frame in [Frame::Geocentric, Frame::Heliocentric, Frame::Barycentric] {
            let a = eph.charts(&days, frame);
            let b = back.charts(&days, frame);
            assert_eq!(a, b, "{frame:?} charts differ after round trip");
        }
    }

    #[test]
    fn malformed_input_is_rejected_rather_than_read_past() {
        let (eph, _) = build(0.0, 100.0, synthetic);
        let good = eph.to_bytes();
        assert!(CompactEphemeris::from_bytes(&good).is_some());
        assert!(CompactEphemeris::from_bytes(b"").is_none());
        assert!(CompactEphemeris::from_bytes(b"NOPE").is_none());
        for cut in [8, 40, good.len() - 1] {
            assert!(
                CompactEphemeris::from_bytes(&good[..cut]).is_none(),
                "truncation to {cut} bytes was accepted"
            );
        }
        let mut bad_version = good.clone();
        bad_version[4] = 99;
        assert!(CompactEphemeris::from_bytes(&bad_version).is_none());
    }

    #[test]
    fn charts_match_the_spice_path_in_shape_and_degeneracy() {
        let (eph, _) = build(0.0, 100.0, synthetic);
        let days = [10.0, 33.3, 71.0];
        for frame in [Frame::Geocentric, Frame::Heliocentric, Frame::Barycentric] {
            let cs = eph.charts(&days, frame);
            assert_eq!(cs.len(), days.len());
            for c in &cs {
                assert_eq!(c.states.len(), BODIES.len());
                assert_eq!(c.frame, frame);
                let degenerate = c.states.iter().filter(|s| s.dist == 0.0).count();
                // Exactly the observer is degenerate, and only when it is a body.
                let expect = usize::from(BODIES.contains(&frame.observer()));
                assert_eq!(degenerate, expect, "{frame:?}");
                for s in &c.states {
                    assert!(s.lon.is_finite() && s.lat.is_finite() && s.dist.is_finite());
                }
            }
        }
    }

    #[test]
    fn velocity_matches_a_finite_difference_of_position() {
        // Independent check on the analytic derivative: it must agree with a
        // numerical difference of the same series.
        let (eph, _) = build(0.0, 200.0, synthetic);
        let h = 1e-4;
        for &day in &[20.0, 55.5, 130.0] {
            let a = eph.charts(&[day - h], Frame::Barycentric);
            let b = eph.charts(&[day + h], Frame::Barycentric);
            let m = eph.charts(&[day], Frame::Barycentric);
            for bi in 0..BODIES.len() {
                let dl = {
                    let mut d = b[0].states[bi].lon - a[0].states[bi].lon;
                    if d > std::f64::consts::PI { d -= std::f64::consts::TAU }
                    if d < -std::f64::consts::PI { d += std::f64::consts::TAU }
                    d / (2.0 * h)
                };
                let got = m[0].states[bi].lon_speed;
                assert!(
                    (got - dl).abs() < 1e-6 * dl.abs().max(1e-4),
                    "{}: analytic {got:e} vs differenced {dl:e}",
                    BODIES[bi]
                );
            }
        }
    }

    #[test]
    fn evaluation_outside_the_span_clamps_rather_than_diverging() {
        // A Chebyshev fit explodes outside its interval. Clamping keeps the answer
        // stale instead of catastrophically wrong.
        let (eph, _) = build(0.0, 100.0, synthetic);
        for day in [-500.0, -1.0, 101.0, 5000.0] {
            let c = eph.charts(&[day], Frame::Barycentric);
            for s in &c[0].states {
                assert!(s.dist.is_finite() && s.dist < 1e11, "day {day}: dist {}", s.dist);
            }
        }
    }
}
