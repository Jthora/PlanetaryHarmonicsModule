//! Lunar cycles, eclipses, stations, and fixed points.
//!
//! # Phases without orbital elements
//!
//! The Moon has three periods that matter and they beat against each other: the
//! synodic (29.53 d, phase), the anomalistic (27.55 d, distance), and the draconic
//! (27.21 d, latitude). Eclipses happen where all three line up. Getting each as a
//! *continuous phase* rather than a raw value matters for a model, because raw
//! distance alone cannot distinguish approach from recession.
//!
//! The usual route is osculating orbital elements. This module uses a cheaper and
//! more robust one: each cycle is recovered as the phase of an analytic signal
//! built from the value and its own time derivative, both of which
//! [`crate::chart`] already supplies exactly from the state vector —
//!
//! ```text
//! anomalistic  =  atan2( ḋ/ω_a , −(d − d̄) )      0 at perigee
//! draconic     =  atan2( β , β̇/ω_d )             0 at the ascending node
//! ```
//!
//! No element conversion, no iteration, and the phase is continuous through the
//! turning points where a value-only feature would fold. The ascending node then
//! falls out as `λ − F`, which is checked against its known 18.6-year regression.
//!
//! # Eclipses as a continuous quantity
//!
//! A boolean "is there an eclipse" throws away almost everything: a near-miss at
//! 1.2° and a quiet quadrature at 90° would look identical. What is given instead
//! is the true angular separation from syzygy, from which a model can learn its own
//! threshold — along with smooth scores for callers that want one.

use crate::chart::{Chart, BODIES};
use crate::chart_features::FeatureSet;
use std::f64::consts::{PI, TAU};

/// Mean lunar distance, km.
const MOON_MEAN_DIST: f64 = 384_400.0;
/// Anomalistic month, days — perigee to perigee.
const ANOMALISTIC_MONTH: f64 = 27.554_549_9;
/// Draconic month, days — node to node.
const DRACONIC_MONTH: f64 = 27.212_220_8;

/// Ecliptic longitude and latitude of the galactic centre at J2000, radians.
///
/// From the IAU galactic pole definition; the centre sits at RA 17h45m37s,
/// dec −28°56′10″, which transforms to these ecliptic coordinates.
pub const GALACTIC_CENTRE: (f64, f64) = (266.840_0_f64 * PI / 180.0, -5.536_4_f64 * PI / 180.0);

fn norm_tau(x: f64) -> f64 {
    let r = x % TAU;
    if r < 0.0 { r + TAU } else { r }
}

/// Angular separation between two directions given as (longitude, latitude).
pub fn separation(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let (s1, c1) = lat1.sin_cos();
    let (s2, c2) = lat2.sin_cos();
    let d = s1 * s2 + c1 * c2 * (lon1 - lon2).cos();
    d.clamp(-1.0, 1.0).acos()
}

/// The Moon's anomalistic phase, radians — 0 at perigee, π at apogee.
pub fn anomalistic_phase(dist: f64, dist_speed: f64) -> f64 {
    let omega = TAU / ANOMALISTIC_MONTH;
    norm_tau((dist_speed / omega).atan2(-(dist - MOON_MEAN_DIST)))
}

/// The Moon's draconic phase (argument of latitude), radians — 0 at the ascending node.
pub fn draconic_phase(lat: f64, lat_speed: f64) -> f64 {
    let omega = TAU / DRACONIC_MONTH;
    norm_tau(lat.atan2(lat_speed / omega))
}

/// Lunar-cycle features. Requires SUN and MOON in the chart; empty otherwise.
pub fn lunar(chart: &Chart, max_harmonic: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    let (Some(moon), Some(sun)) = (chart.body("MOON"), chart.body("SUN")) else {
        return f;
    };
    if moon.dist == 0.0 || sun.dist == 0.0 {
        return f;
    }

    // Synodic: elongation, at every harmonic. Harmonic 2 is the one the tides use.
    let d = norm_tau(moon.lon - sun.lon);
    let (mut c, mut s) = (1.0, 0.0);
    let (ds, dc) = d.sin_cos();
    for h in 1..=max_harmonic {
        let (nc, ns) = (c * dc - s * ds, s * dc + c * ds);
        c = nc;
        s = ns;
        f.push(format!("syn.h{h}.cos"), c);
        f.push(format!("syn.h{h}.sin"), s);
    }
    // Illuminated fraction: the phase as it is actually seen.
    f.push("illumination", (1.0 - d.cos()) / 2.0);

    let anom = anomalistic_phase(moon.dist, moon.dist_speed);
    let drac = draconic_phase(moon.lat, moon.lat_speed);
    for (name, phase) in [("anom", anom), ("drac", drac)] {
        let (mut c, mut s) = (1.0, 0.0);
        let (ps, pc) = phase.sin_cos();
        for h in 1..=max_harmonic.min(6) {
            let (nc, ns) = (c * pc - s * ps, s * pc + c * ps);
            c = nc;
            s = ns;
            f.push(format!("{name}.h{h}.cos"), c);
            f.push(format!("{name}.h{h}.sin"), s);
        }
    }

    // Raw quantities alongside the phases: distance drives the tidal amplitude as
    // 1/d^3, and declination sets how the diurnal tide splits between hemispheres.
    f.push("dist", moon.dist);
    f.push("dist_norm", (moon.dist - MOON_MEAN_DIST) / 21_000.0);
    f.push("tidal_scale", (MOON_MEAN_DIST / moon.dist).powi(3));
    f.push("lat", moon.lat);
    f.push("dec", moon.dec);
    f.push("dec_abs", moon.dec.abs());
    f.push("lon_speed", moon.lon_speed);

    // Ascending node, and where the Moon and Sun sit relative to it. The node's
    // 18.6-year regression is the slowest lunar cycle and modulates the others.
    let node = norm_tau(moon.lon - drac);
    f.push("node.cos", node.cos());
    f.push("node.sin", node.sin());
    let sun_from_node = norm_tau(sun.lon - node);
    f.push("sun_node.cos", sun_from_node.cos());
    f.push("sun_node.sin", sun_from_node.sin());
    // Twice the Sun-node angle: eclipse seasons recur twice per draconic year, so
    // this is the term that is actually near its extremum when eclipses are possible.
    f.push("sun_node2.cos", (2.0 * sun_from_node).cos());
    f.push("sun_node2.sin", (2.0 * sun_from_node).sin());
    f
}

/// Eclipse-geometry features.
///
/// `sep_solar` is the true Moon-Sun angular separation and `sep_lunar` the Moon's
/// separation from the antisolar point; an eclipse of some kind requires one of
/// them below roughly 1.5°. The `*_score` companions are smooth bumps of width
/// 1° for callers that want a single number, but the separations are the honest
/// features — a model can pick its own threshold from them.
pub fn eclipses(chart: &Chart) -> FeatureSet {
    let mut f = FeatureSet::default();
    let (Some(moon), Some(sun)) = (chart.body("MOON"), chart.body("SUN")) else {
        return f;
    };
    if moon.dist == 0.0 || sun.dist == 0.0 {
        return f;
    }
    let sep_solar = separation(moon.lon, moon.lat, sun.lon, sun.lat);
    let sep_lunar = separation(moon.lon, moon.lat, sun.lon + PI, -sun.lat);
    let width = 1.0_f64.to_radians();
    f.push("sep_solar", sep_solar);
    f.push("sep_lunar", sep_lunar);
    f.push("score_solar", (-(sep_solar / width).powi(2)).exp());
    f.push("score_lunar", (-(sep_lunar / width).powi(2)).exp());
    f.push("sep_min", sep_solar.min(sep_lunar));
    // Apparent radii decide whether a central solar eclipse is total or annular.
    // Ratio > 1 means the Moon covers the Sun.
    let moon_radius = (1_737.4_f64 / moon.dist).asin();
    let sun_radius = (696_000.0_f64 / sun.dist).asin();
    f.push("radius_ratio", moon_radius / sun_radius);
    f.push("umbral_margin", moon_radius + sun_radius - sep_solar);
    f
}

/// Aspects from every body to the galactic centre and anticentre.
pub fn fixed_points(chart: &Chart, max_harmonic: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    // The galactic centre is fixed in space, so in an of-date zodiac it would drift
    // with precession -- but body longitudes and this constant are both J2000, so
    // the difference is already precession-free and no correction belongs on either.
    for (point_name, lon) in [("gc", GALACTIC_CENTRE.0), ("gac", GALACTIC_CENTRE.0 + PI)] {
        for (bi, body) in BODIES.iter().enumerate() {
            let s = &chart.states[bi];
            if s.dist == 0.0 {
                continue;
            }
            let d = norm_tau(s.lon - lon);
            let (mut c, mut sn) = (1.0, 0.0);
            let (ds, dc) = d.sin_cos();
            for h in 1..=max_harmonic {
                let (nc, ns) = (c * dc - sn * ds, sn * dc + c * ds);
                c = nc;
                sn = ns;
                f.push(format!("{point_name}.{body}.h{h}.cos"), c);
                f.push(format!("{point_name}.{body}.h{h}.sin"), sn);
            }
        }
    }
    f
}

/// A retrograde station: the instant a body's apparent longitude stops advancing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Station {
    /// Index into [`BODIES`].
    pub body: usize,
    /// Day of the station, interpolated between the bracketing samples.
    pub day: f64,
    /// True if the body is entering retrograde (speed going positive to negative).
    pub retrograde: bool,
}

/// Locate every station in a chart series.
///
/// Stations are sign changes in `lon_speed`, refined by linear interpolation.
/// **The series must be sampled finely enough that no two stations of the same body
/// fall inside one step** — daily sampling is ample, since the shortest retrograde
/// (Mercury's) lasts about three weeks.
pub fn stations(charts: &[Chart]) -> Vec<Station> {
    let mut out = Vec::new();
    for bi in 0..BODIES.len() {
        for w in charts.windows(2) {
            let (a, b) = (&w[0].states[bi], &w[1].states[bi]);
            if a.dist == 0.0 || b.dist == 0.0 {
                continue;
            }
            if a.lon_speed == 0.0 || a.lon_speed.signum() == b.lon_speed.signum() {
                continue;
            }
            let frac = a.lon_speed / (a.lon_speed - b.lon_speed);
            out.push(Station {
                body: bi,
                day: w[0].day + frac * (w[1].day - w[0].day),
                retrograde: a.lon_speed > 0.0,
            });
        }
    }
    out.sort_by(|x, y| x.day.partial_cmp(&y.day).unwrap());
    out
}

/// Per-epoch station-timing features for a chart series.
///
/// For each body: days since the previous station and days until the next. Where
/// the series does not reach a station in one direction the value is the distance
/// to the end of the series, negated in sign meaning — callers wanting to exclude
/// those should trim the series ends rather than trust a saturated value.
///
/// Returns one [`FeatureSet`] per chart, aligned with the input.
pub fn station_timing(charts: &[Chart]) -> Vec<FeatureSet> {
    let all = stations(charts);
    let mut out: Vec<FeatureSet> = charts.iter().map(|_| FeatureSet::default()).collect();
    for (bi, body) in BODIES.iter().enumerate() {
        let times: Vec<f64> = all.iter().filter(|s| s.body == bi).map(|s| s.day).collect();
        for (ci, chart) in charts.iter().enumerate() {
            let d = chart.day;
            let prev = times.iter().rev().find(|&&t| t <= d).copied();
            let next = times.iter().find(|&&t| t >= d).copied();
            let since = prev.map(|t| d - t).unwrap_or(d - charts[0].day);
            let until = next.map(|t| t - d).unwrap_or(charts[charts.len() - 1].day - d);
            out[ci].push(format!("since.{body}"), since);
            out[ci].push(format!("until.{body}"), until);
            out[ci].push(format!("nearest.{body}"), since.min(until));
            // Retrograde motion is the state between stations; the fraction of the
            // way through the current interval places the body within it.
            let span = since + until;
            out[ci].push(
                format!("frac.{body}"),
                if span > 0.0 { since / span } else { 0.0 },
            );
        }
    }
    out
}

/// Lunar, eclipse and fixed-point features, prefixed and merged.
pub fn all(chart: &Chart, max_harmonic: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    for (prefix, set) in [
        ("moon", lunar(chart, max_harmonic)),
        ("ecl", eclipses(chart)),
        ("fix", fixed_points(chart, max_harmonic.min(6))),
    ] {
        for (n, v) in set.names.iter().zip(&set.values) {
            f.push(format!("{prefix}.{n}"), *v);
        }
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chart::{BodyState, Frame};

    /// A synthetic Moon on a circular inclined orbit, so the cycles are known exactly.
    fn moon_chart(day: f64) -> Chart {
        let mut states = vec![BodyState::default(); BODIES.len()];
        let mi = BODIES.iter().position(|b| *b == "MOON").unwrap();
        let si = BODIES.iter().position(|b| *b == "SUN").unwrap();

        let inc = 5.145_f64.to_radians();
        // Argument of latitude advances at the draconic rate; the node regresses.
        let f_ang = TAU * day / DRACONIC_MONTH;
        let lat = (inc.sin() * f_ang.sin()).asin();
        let lat_speed = inc.sin() * f_ang.cos() * (TAU / DRACONIC_MONTH) / lat.cos().max(1e-12);
        let anom = TAU * day / ANOMALISTIC_MONTH;
        let dist = MOON_MEAN_DIST - 21_000.0 * anom.cos();
        let dist_speed = 21_000.0 * anom.sin() * (TAU / ANOMALISTIC_MONTH);

        states[mi] = BodyState {
            lon: norm_tau(TAU * day / 27.321_582),
            lat,
            dist,
            lon_speed: TAU / 27.321_582,
            lat_speed,
            dist_speed,
            ra: 0.0,
            dec: lat,
        };
        states[si] = BodyState {
            lon: norm_tau(TAU * day / 365.256_363),
            lat: 0.0,
            dist: 1.496e8,
            lon_speed: TAU / 365.256_363,
            lat_speed: 0.0,
            dist_speed: 0.0,
            ra: 0.0,
            dec: 0.0,
        };
        Chart { day, frame: Frame::Geocentric, states }
    }

    #[test]
    fn anomalistic_phase_is_zero_at_perigee_and_advances() {
        // Perigee: minimum distance, zero radial speed.
        let p = anomalistic_phase(MOON_MEAN_DIST - 21_000.0, 0.0);
        assert!(p < 1e-12 || (TAU - p) < 1e-12, "perigee phase {p}");
        // Apogee is half a cycle away.
        let a = anomalistic_phase(MOON_MEAN_DIST + 21_000.0, 0.0);
        assert!((a - PI).abs() < 1e-12, "apogee phase {a}");
        // Just after perigee the phase must have increased, not decreased.
        let after = anomalistic_phase(MOON_MEAN_DIST - 20_000.0, 2000.0);
        assert!(after > 0.0 && after < PI / 2.0, "phase just after perigee {after}");
    }

    #[test]
    fn draconic_phase_is_zero_at_the_ascending_node() {
        let n = draconic_phase(0.0, 1.0);
        assert!(n.abs() < 1e-12, "ascending node phase {n}");
        // Maximum latitude is a quarter cycle later.
        let hi = draconic_phase(0.09, 0.0);
        assert!((hi - PI / 2.0).abs() < 1e-12, "max-latitude phase {hi}");
        // Descending node, half a cycle on.
        let d = draconic_phase(0.0, -1.0);
        assert!((d - PI).abs() < 1e-12, "descending node phase {d}");
    }

    #[test]
    fn phases_advance_at_their_own_periods() {
        // The point of the analytic-signal construction: each phase must complete
        // exactly one turn per its own period, independent of the others.
        for (period, extract) in [
            (ANOMALISTIC_MONTH, 0usize),
            (DRACONIC_MONTH, 1usize),
        ] {
            let mut total = 0.0;
            let mut prev: Option<f64> = None;
            let steps = 200;
            for k in 0..=steps {
                let day = 100.0 + period * k as f64 / steps as f64;
                let c = moon_chart(day);
                let m = c.body("MOON").unwrap();
                let ph = if extract == 0 {
                    anomalistic_phase(m.dist, m.dist_speed)
                } else {
                    draconic_phase(m.lat, m.lat_speed)
                };
                if let Some(p) = prev {
                    let mut d = ph - p;
                    if d > PI { d -= TAU }
                    if d < -PI { d += TAU }
                    total += d;
                }
                prev = Some(ph);
            }
            assert!(
                (total - TAU).abs() < 1e-3,
                "period {period}: phase swept {total} rad, expected {TAU}"
            );
        }
    }

    #[test]
    fn the_lunar_node_regresses_once_in_186_years() {
        // The strongest available check on the draconic construction: the node is
        // derived as lambda - F, and its regression rate is not put in anywhere.
        // If it comes out at -19.34 deg/yr the whole chain is right.
        let mut total = 0.0;
        let mut prev: Option<f64> = None;
        let (t0, t1) = (0.0, 2000.0);
        let steps = 4000;
        for k in 0..=steps {
            let day = t0 + (t1 - t0) * k as f64 / steps as f64;
            let c = moon_chart(day);
            let f = lunar(&c, 1);
            let node = f.get("node.sin").unwrap().atan2(f.get("node.cos").unwrap());
            if let Some(p) = prev {
                let mut d = node - p;
                if d > PI { d -= TAU }
                if d < -PI { d += TAU }
                total += d;
            }
            prev = Some(node);
        }
        let deg_per_year = total.to_degrees() / ((t1 - t0) / 365.25);
        assert!(
            (deg_per_year + 19.34).abs() < 0.6,
            "node moves {deg_per_year} deg/yr, expected -19.34"
        );
    }

    #[test]
    fn separation_is_symmetric_and_bounded() {
        let cases = [
            (0.0, 0.0, PI, 0.0, PI),
            (0.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, PI / 2.0, 0.0, PI / 2.0),
            (0.0, PI / 2.0, 1.0, -PI / 2.0, PI),
        ];
        for (l1, b1, l2, b2, want) in cases {
            let s = separation(l1, b1, l2, b2);
            assert!((s - want).abs() < 1e-12, "{s} vs {want}");
            assert!((separation(l2, b2, l1, b1) - s).abs() < 1e-12, "not symmetric");
        }
    }

    #[test]
    fn eclipse_separation_is_small_only_near_syzygy_at_low_latitude() {
        // A new moon far from the node must NOT read as an eclipse -- this is the
        // failure a longitude-only test would let through.
        let mut states = vec![BodyState::default(); BODIES.len()];
        let mi = BODIES.iter().position(|b| *b == "MOON").unwrap();
        let si = BODIES.iter().position(|b| *b == "SUN").unwrap();
        states[si] = BodyState { lon: 1.0, lat: 0.0, dist: 1.496e8, ..Default::default() };

        // Conjunction in longitude, but 5 degrees of latitude away: no eclipse.
        states[mi] = BodyState {
            lon: 1.0,
            lat: 5.0_f64.to_radians(),
            dist: MOON_MEAN_DIST,
            ..Default::default()
        };
        let c = Chart { day: 0.0, frame: Frame::Geocentric, states: states.clone() };
        let f = eclipses(&c);
        assert!(f.get("sep_solar").unwrap() > 4.9_f64.to_radians());
        assert!(f.get("score_solar").unwrap() < 1e-6, "false eclipse off-node");
        assert!(f.get("umbral_margin").unwrap() < 0.0, "off-node margin should be negative");

        // Same conjunction, on the node: eclipse.
        states[mi].lat = 0.0;
        let c = Chart { day: 0.0, frame: Frame::Geocentric, states };
        let f = eclipses(&c);
        assert!(f.get("sep_solar").unwrap() < 1e-9);
        assert!((f.get("score_solar").unwrap() - 1.0).abs() < 1e-12);
        assert!(f.get("umbral_margin").unwrap() > 0.0, "on-node margin should be positive");
        // Moon and Sun subtend nearly the same angle -- the reason total eclipses
        // are barely total. Ratio near 1 confirms the radii are in the right units.
        let ratio = f.get("radius_ratio").unwrap();
        assert!((0.9..1.1).contains(&ratio), "apparent radius ratio {ratio}");
    }

    #[test]
    fn stations_are_found_where_longitude_speed_reverses() {
        // Build a series whose speed reverses at two known days.
        let mi = BODIES.iter().position(|b| *b == "MERCURY BARYCENTER" || *b == "MERCURY")
            .unwrap_or(0);
        let charts: Vec<Chart> = (0..120)
            .map(|k| {
                let day = k as f64;
                let mut states = vec![BodyState::default(); BODIES.len()];
                for s in states.iter_mut() {
                    s.dist = 1.0;
                }
                // Speed crosses zero at day 30 (going direct) and day 90 (going
                // retrograde), both strictly inside the series.
                states[mi].lon_speed = (TAU * (day - 30.0) / 120.0).sin();
                Chart { day, frame: Frame::Geocentric, states }
            })
            .collect();
        let all_stations = stations(&charts);
        let found: Vec<&Station> = all_stations.iter().filter(|s| s.body == mi).collect();
        assert_eq!(found.len(), 2, "found {found:?}");
        assert!((found[0].day - 30.0).abs() < 0.1, "{:?}", found[0]);
        assert!((found[1].day - 90.0).abs() < 0.1, "{:?}", found[1]);
        assert!(!found[0].retrograde, "speed turning positive is the direct station");
        assert!(found[1].retrograde, "speed turning negative enters retrograde");

        let timing = station_timing(&charts);
        let body = BODIES[mi];
        assert!(timing[30].get(&format!("nearest.{body}")).unwrap() < 0.2);
        // Day 45 sits 15 days after the first station and 45 before the second.
        assert!((timing[45].get(&format!("since.{body}")).unwrap() - 15.0).abs() < 0.2);
        assert!((timing[45].get(&format!("until.{body}")).unwrap() - 45.0).abs() < 0.2);
        assert!((timing[45].get(&format!("frac.{body}")).unwrap() - 0.25).abs() < 0.02);
    }

    #[test]
    fn every_feature_is_finite_and_named_once() {
        for k in 0..40 {
            let c = moon_chart(k as f64 * 3.7);
            let f = all(&c, 8);
            for (n, v) in f.names.iter().zip(&f.values) {
                assert!(v.is_finite(), "{n} = {v}");
            }
            let unique: std::collections::HashSet<_> = f.names.iter().collect();
            assert_eq!(unique.len(), f.names.len());
        }
    }
}
