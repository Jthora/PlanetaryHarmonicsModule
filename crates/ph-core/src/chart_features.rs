//! Chart features — everything derivable from [`crate::chart`] primitives.
//!
//! Aspects, declination aspects, chart shape and whole-chart resonance, all as
//! cheap trigonometry on stored primitives rather than materialised columns.
//!
//! # No filtering
//!
//! Nothing here is selected for having a physical story. Every pair, every
//! harmonic order, every frame. What survives is decided downstream by
//! validation, not upstream by plausibility.
//!
//! # Every feature is named
//!
//! Values come with names in the same order. That is not decoration: when a model
//! reports skill, the first question is *which columns*, and an unlabelled vector
//! cannot answer it.

use crate::chart::{Chart, BODIES};
use std::f64::consts::TAU;

/// Named feature values.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FeatureSet {
    pub names: Vec<String>,
    pub values: Vec<f64>,
}

impl FeatureSet {
    pub fn push(&mut self, name: impl Into<String>, value: f64) {
        self.names.push(name.into());
        self.values.push(value);
    }

    pub fn extend(&mut self, other: FeatureSet) {
        self.names.extend(other.names);
        self.values.extend(other.values);
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Look a feature up by name. Linear, for diagnostics rather than hot paths.
    pub fn get(&self, name: &str) -> Option<f64> {
        self.names.iter().position(|n| n == name).map(|i| self.values[i])
    }
}

fn short(body: &str) -> &str {
    body.split(' ').next().unwrap_or(body)
}

/// Aspect features: `cos(nΔλ)` and `sin(nΔλ)` for every body pair and every
/// harmonic order up to `max_harmonic`.
///
/// The classical aspects are the low orders of this basis — conjunction and
/// opposition at n = 1 and 2, trine at 3, square at 4, sextile at 6, the twelvefold
/// division at 12. Carrying n = 1…24 contains all of them **plus every base nobody
/// has tried**, which is the point.
///
/// Keeping the sine term is what lets a model learn a phase offset. Dropping it
/// would assume the response peaks at exact alignment.
pub fn aspects(chart: &Chart, max_harmonic: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    let lons = chart.longitudes();
    for i in 0..BODIES.len() {
        for j in (i + 1)..BODIES.len() {
            let d = lons[j] - lons[i];
            let (a, b) = (short(BODIES[i]), short(BODIES[j]));
            let (s1, c1) = d.sin_cos();
            let (mut c, mut s) = (c1, s1);
            for n in 1..=max_harmonic {
                f.push(format!("asp.{a}.{b}.h{n}.cos"), c);
                f.push(format!("asp.{a}.{b}.h{n}.sin"), s);
                if n < max_harmonic {
                    // Angle-sum recurrence: O(n) multiply-adds, not n trig calls.
                    let (cn, sn) = (c * c1 - s * s1, s * c1 + c * s1);
                    c = cn;
                    s = sn;
                }
            }
        }
    }
    f
}

/// Declination aspects — parallels and contraparallels.
///
/// Two bodies are *parallel* at equal declination and *contraparallel* at equal and
/// opposite. Long used in astrology, almost never harmonically encoded. Emitted as
/// smooth closeness measures rather than orb thresholds, so no arbitrary cut is
/// imposed.
pub fn declination_aspects(chart: &Chart) -> FeatureSet {
    let mut f = FeatureSet::default();
    let decs = chart.declinations();
    for i in 0..BODIES.len() {
        for j in (i + 1)..BODIES.len() {
            let (a, b) = (short(BODIES[i]), short(BODIES[j]));
            f.push(format!("dec.{a}.{b}.parallel"), (decs[i] - decs[j]).cos());
            f.push(format!("dec.{a}.{b}.contra"), (decs[i] + decs[j]).cos());
            f.push(format!("dec.{a}.{b}.diff"), decs[i] - decs[j]);
        }
    }
    f
}

/// Chart shape — circular statistics on the longitude distribution.
///
/// The classical chart shapes are distributional facts wearing names. A *bundle*
/// is high concentration with one large gap; a *splash* is low concentration with
/// no large gap; a *bucket* is a bowl plus an isolated body, which is a large gap
/// next to a small one.
///
/// Emitting the statistics rather than the categories keeps them continuous and
/// avoids inventing thresholds. Cheap, and as far as I know untried — which makes
/// them the family least contaminated by prior expectation.
pub fn shape(chart: &Chart) -> FeatureSet {
    let mut f = FeatureSet::default();
    let mut lons = chart.longitudes();
    let n = lons.len() as f64;

    // Resultant vector: length is concentration, 0 for uniform and 1 for coincident.
    let (mut sx, mut sy) = (0.0f64, 0.0f64);
    for &l in &lons {
        let (s, c) = l.sin_cos();
        sx += c;
        sy += s;
    }
    let (rx, ry) = (sx / n, sy / n);
    let r = rx.hypot(ry);
    f.push("shape.concentration", r);
    f.push("shape.circular_variance", 1.0 - r);
    // Mean direction as its components, so it stays continuous across the wrap.
    f.push("shape.mean_dir.cos", if r > 0.0 { rx / r } else { 0.0 });
    f.push("shape.mean_dir.sin", if r > 0.0 { ry / r } else { 0.0 });

    // Gap structure.
    lons.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut gaps: Vec<f64> = (0..lons.len())
        .map(|i| {
            let next = if i + 1 == lons.len() {
                lons[0] + TAU
            } else {
                lons[i + 1]
            };
            next - lons[i]
        })
        .collect();
    gaps.sort_by(|a, b| b.partial_cmp(a).unwrap());
    f.push("shape.max_gap", gaps[0]);
    f.push("shape.second_gap", gaps[1]);
    // Occupied span: a bowl fills half the circle, a bundle much less.
    f.push("shape.span", TAU - gaps[0]);
    // Clusters, counted as gaps exceeding twice the mean gap.
    let mean_gap = TAU / n;
    f.push(
        "shape.n_clusters",
        gaps.iter().filter(|&&g| g > 2.0 * mean_gap).count() as f64,
    );
    // Entropy of the gap distribution: high for even spacing, low for clumping.
    let ent: f64 = gaps
        .iter()
        .map(|&g| {
            let p = g / TAU;
            if p > 0.0 {
                -p * p.ln()
            } else {
                0.0
            }
        })
        .sum();
    f.push("shape.gap_entropy", ent);
    f
}

/// Whole-chart harmonic resonance — the CosmicCypher `base_N` idea, corrected.
///
/// Project every pair's separation onto a base-N harmonic and average across the
/// chart, giving one number per base: how strongly the whole configuration rings
/// at that division.
///
/// Two fixes against the original. Its `base7` and `base11` computed
/// `sin(θ·π/(1/7))`, a period of 0.29° — noise, so two of twelve bases never
/// worked. And it emitted sine only, which is *zero* at exact aspect; both
/// components are emitted here so the phase is learnable rather than assumed.
pub fn resonance(chart: &Chart, max_base: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    let lons = chart.longitudes();
    let mut pairs = 0.0f64;
    let mut acc = vec![(0.0f64, 0.0f64); max_base + 1];
    for i in 0..lons.len() {
        for j in (i + 1)..lons.len() {
            let d = lons[j] - lons[i];
            pairs += 1.0;
            for (base, slot) in acc.iter_mut().enumerate().skip(1) {
                let (s, c) = (d * base as f64).sin_cos();
                slot.0 += c;
                slot.1 += s;
            }
        }
    }
    for (base, (c, s)) in acc.iter().enumerate().skip(1) {
        f.push(format!("res.base{base}.cos"), c / pairs);
        f.push(format!("res.base{base}.sin"), s / pairs);
        f.push(format!("res.base{base}.mag"), (c / pairs).hypot(s / pairs));
    }
    f
}

/// Motion features — retrograde structure and speed.
pub fn motion(chart: &Chart) -> FeatureSet {
    let mut f = FeatureSet::default();
    let mut retro = 0.0;
    for (i, s) in chart.states.iter().enumerate() {
        let b = short(BODIES[i]);
        f.push(format!("mot.{b}.lon_speed"), s.lon_speed);
        f.push(format!("mot.{b}.lat_speed"), s.lat_speed);
        f.push(format!("mot.{b}.dist"), s.dist);
        f.push(format!("mot.{b}.dist_speed"), s.dist_speed);
        f.push(format!("mot.{b}.retro"), if s.lon_speed < 0.0 { 1.0 } else { 0.0 });
        // Proximity to a station: small |speed| means near a turning point.
        f.push(format!("mot.{b}.station"), (-s.lon_speed.abs() * 20.0).exp());
        if s.lon_speed < 0.0 {
            retro += 1.0;
        }
    }
    f.push("mot.n_retrograde", retro);
    f
}

/// Everything, for one chart in one frame.
///
/// Names are prefixed by frame so several frames can be concatenated without
/// collision.
pub fn all(chart: &Chart, max_harmonic: usize, max_base: usize) -> FeatureSet {
    let tag = match chart.frame {
        crate::chart::Frame::Geocentric => "geo",
        crate::chart::Frame::Heliocentric => "helio",
        crate::chart::Frame::Barycentric => "bary",
    };
    let mut out = FeatureSet::default();
    for part in [
        aspects(chart, max_harmonic),
        declination_aspects(chart),
        shape(chart),
        resonance(chart, max_base),
        motion(chart),
    ] {
        for (n, v) in part.names.into_iter().zip(part.values) {
            out.push(format!("{tag}.{n}"), v);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chart::{BodyState, Frame};

    fn chart_with(lons: &[f64]) -> Chart {
        let mut states = vec![BodyState::default(); BODIES.len()];
        for (i, &l) in lons.iter().enumerate().take(BODIES.len()) {
            states[i].lon = l;
        }
        Chart {
            day: 0.0,
            frame: Frame::Geocentric,
            states,
        }
    }

    #[test]
    fn aspect_count_and_naming() {
        let c = chart_with(&[0.0; 11]);
        let f = aspects(&c, 12);
        let pairs = BODIES.len() * (BODIES.len() - 1) / 2;
        assert_eq!(f.len(), pairs * 12 * 2);
        assert!(f.names.iter().any(|n| n.starts_with("asp.SUN.MOON.h1.")));
    }

    #[test]
    fn aspect_recurrence_matches_direct_evaluation() {
        let mut lons = vec![0.0; BODIES.len()];
        lons[0] = 0.3;
        lons[1] = 2.1;
        let c = chart_with(&lons);
        let f = aspects(&c, 16);
        let d = 2.1 - 0.3;
        for n in 1..=16 {
            let got = f.get(&format!("asp.SUN.EARTH.h{n}.cos")).unwrap();
            assert!((got - (n as f64 * d).cos()).abs() < 1e-12, "order {n}");
        }
    }

    #[test]
    fn conjunction_and_opposition_sit_at_the_expected_orders() {
        // All bodies conjunct: every harmonic cosine is 1.
        let conj = aspects(&chart_with(&[1.0; 11]), 12);
        assert!(conj
            .names
            .iter()
            .zip(&conj.values)
            .filter(|(n, _)| n.ends_with(".cos"))
            .all(|(_, v)| (v - 1.0).abs() < 1e-12));

        // Two bodies opposed: h1 cosine is -1, h2 cosine is +1.
        let mut l = vec![0.0; BODIES.len()];
        l[1] = std::f64::consts::PI;
        let opp = aspects(&chart_with(&l), 4);
        assert!((opp.get("asp.SUN.EARTH.h1.cos").unwrap() + 1.0).abs() < 1e-12);
        assert!((opp.get("asp.SUN.EARTH.h2.cos").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn shape_separates_bundle_from_splash() {
        // Bundle: everything within a narrow span.
        let bundle: Vec<f64> = (0..BODIES.len()).map(|i| i as f64 * 0.02).collect();
        let b = shape(&chart_with(&bundle));
        // Splash: evenly spread around the circle.
        let splash: Vec<f64> = (0..BODIES.len())
            .map(|i| i as f64 * TAU / BODIES.len() as f64)
            .collect();
        let s = shape(&chart_with(&splash));

        assert!(
            b.get("shape.concentration").unwrap() > 0.99,
            "bundle concentration {}",
            b.get("shape.concentration").unwrap()
        );
        assert!(
            s.get("shape.concentration").unwrap() < 0.01,
            "splash concentration {}",
            s.get("shape.concentration").unwrap()
        );
        assert!(b.get("shape.max_gap").unwrap() > s.get("shape.max_gap").unwrap());
        // Even spacing maximises gap entropy.
        assert!(s.get("shape.gap_entropy").unwrap() > b.get("shape.gap_entropy").unwrap());
    }

    #[test]
    fn shape_detects_a_bucket() {
        // Bowl plus one isolated body: one very large gap next to small ones.
        let mut l: Vec<f64> = (0..BODIES.len() - 1).map(|i| i as f64 * 0.15).collect();
        l.push(std::f64::consts::PI);
        let f = shape(&chart_with(&l));
        assert!(f.get("shape.max_gap").unwrap() > 2.0, "handle gap should be large");
        assert!(f.get("shape.n_clusters").unwrap() >= 1.0);
    }

    #[test]
    fn resonance_peaks_at_the_matching_base() {
        // Bodies at multiples of 90 degrees ring at base 4 and its multiples.
        let l: Vec<f64> = (0..BODIES.len())
            .map(|i| (i % 4) as f64 * std::f64::consts::FRAC_PI_2)
            .collect();
        let f = resonance(&chart_with(&l), 12);
        let b4 = f.get("res.base4.mag").unwrap();
        let b3 = f.get("res.base3.mag").unwrap();
        assert!(b4 > 0.9, "base4 magnitude {b4}");
        assert!(b4 > b3, "base4 {b4} should exceed base3 {b3}");
    }

    #[test]
    fn resonance_covers_every_base_including_seven_and_eleven() {
        // The original implementation's base7 and base11 were broken; make sure
        // ours are present and finite.
        let f = resonance(&chart_with(&[0.4; 11]), 12);
        for base in 1..=12 {
            let v = f.get(&format!("res.base{base}.mag")).unwrap();
            assert!(v.is_finite(), "base{base} not finite");
        }
        assert!(f.get("res.base7.mag").unwrap() > 0.9);
        assert!(f.get("res.base11.mag").unwrap() > 0.9);
    }

    #[test]
    fn motion_flags_retrograde_and_stations() {
        let mut c = chart_with(&[0.0; 11]);
        c.states[3].lon_speed = -0.01; // retrograde
        c.states[4].lon_speed = 0.0; // exact station
        c.states[5].lon_speed = 1.0; // fast direct
        let f = motion(&c);
        assert_eq!(f.get("mot.n_retrograde").unwrap(), 1.0);
        assert_eq!(f.get(&format!("mot.{}.retro", short(BODIES[3]))).unwrap(), 1.0);
        assert!((f.get(&format!("mot.{}.station", short(BODIES[4]))).unwrap() - 1.0).abs() < 1e-12);
        assert!(f.get(&format!("mot.{}.station", short(BODIES[5]))).unwrap() < 1e-6);
    }

    #[test]
    fn all_is_frame_prefixed_and_unique() {
        let f = all(&chart_with(&[0.5; 11]), 8, 8);
        assert!(f.names.iter().all(|n| n.starts_with("geo.")));
        let mut seen = std::collections::HashSet::new();
        assert!(f.names.iter().all(|n| seen.insert(n.clone())), "duplicate name");
        assert_eq!(f.names.len(), f.values.len());
        assert!(f.values.iter().all(|v| v.is_finite()));
    }
}

#[cfg(test)]
mod placeholder_tests {
    use super::*;
    use crate::chart::{Chart, Frame};

    /// The names a placeholder yields must be exactly the names a real chart
    /// yields, or `featureNames` in the WASM layer would describe a different
    /// matrix than `chartFeatures` returns and every column would be mislabelled.
    ///
    /// The trap this guards is specific: the derived layers skip bodies at zero
    /// distance, so a naive dummy chart with every body present would include the
    /// frame's observer and shift every column after it.
    #[test]
    fn placeholder_names_match_a_populated_chart_of_the_same_frame() {
        for frame in [Frame::Geocentric, Frame::Heliocentric, Frame::Barycentric] {
            let ph = Chart::placeholder(frame);

            // A "real" chart: same frame, arbitrary but non-degenerate values, and
            // crucially the same body absent.
            let mut real = Chart::placeholder(frame);
            for (i, s) in real.states.iter_mut().enumerate() {
                if s.dist == 0.0 {
                    continue;
                }
                s.lon = 1.7 * i as f64;
                s.lat = -0.03 * i as f64;
                s.dist = 3.3e5 * (i as f64 + 1.0);
                s.lon_speed = -0.02 * i as f64;
            }

            let a = all(&ph, 6, 8);
            let b = all(&real, 6, 8);
            assert_eq!(a.names, b.names, "{frame:?}: placeholder column set differs");
            assert!(!a.names.is_empty());

            // The families differ in how they treat the observer, and that is
            // exactly why a placeholder must reproduce it rather than guess.
            //
            // aspects() and declination_aspects() carry every body unconditionally,
            // so the observer appears there as a constant column -- this is the
            // source of the "constant columns dropped" seen when fitting. The
            // fixed-point family skips zero-distance bodies, so the observer is
            // genuinely absent there. A placeholder that got the zero-distance
            // body wrong would shift every fix.* column.
            //
            // Barycentric is the case that makes this worth asserting rather than
            // assuming: its observer is the solar-system barycentre, which is not
            // a member of BODIES at all, so nothing is skipped and no body is
            // degenerate. Geocentric and heliocentric each lose one.
            let observer = frame.observer();
            let observer_is_a_body = crate::chart::BODIES.contains(&observer);
            assert_eq!(
                observer_is_a_body,
                frame != Frame::Barycentric,
                "{frame:?}: unexpected observer membership"
            );
            if observer_is_a_body {
                assert!(
                    a.names.iter().any(|n| n.contains(&format!(".{observer}."))),
                    "{frame:?}: {observer} should still appear in the unguarded families"
                );
                assert_eq!(
                    ph.states
                        .iter()
                        .filter(|s| s.dist == 0.0)
                        .count(),
                    1,
                    "{frame:?}: exactly one body should be degenerate"
                );
            } else {
                assert!(ph.states.iter().all(|s| s.dist != 0.0),
                        "{frame:?}: no body should be degenerate");
            }
        }
    }

    #[test]
    fn placeholder_names_are_stable_across_calls() {
        // The column order is an API contract once JavaScript caches it.
        let a = all(&Chart::placeholder(Frame::Geocentric), 4, 5);
        let b = all(&Chart::placeholder(Frame::Geocentric), 4, 5);
        assert_eq!(a.names, b.names);
    }
}
