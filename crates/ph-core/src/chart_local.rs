//! Site-local features — the layer that makes the feature vector depend on *where*.
//!
//! Every other feature family is global. At a given hour, `geo.asp.sun_moon.h2.cos`
//! is identical for Tokyo and for Lima. A model fed only those families can learn
//! *when* earthquakes happen but has no channel through which to learn *where*: the
//! only spatial signal available to it is the background rate, which it would then
//! reproduce and nothing more.
//!
//! This module supplies the missing channel. Given a site's latitude and longitude
//! it produces the angles an observer there would measure — sidereal time, the
//! ascendant and midheaven, and each body's hour angle and altitude. These are the
//! same quantities that carry real tidal meaning (the Doodson argument τ *is* the
//! lunar hour angle), which is a reason to expect them to be informative, though
//! nothing here depends on that being true.
//!
//! # Frames of reference
//!
//! Positions from [`crate::chart`] are J2000-referenced. Sidereal time is referred
//! to the equinox *of date*, so mixing the two directly would introduce the full
//! accumulated precession — 0.34° over the catalogue span. [`equatorial_of_date`]
//! precesses the ecliptic longitude and applies the obliquity of date before the
//! conversion, keeping the two consistent.
//!
//! UT1 is approximated by UTC. The difference is bounded by ±0.9 s by construction
//! (leap seconds keep |DUT1| under that), which is 0.004° of Earth rotation.

use crate::chart::{Chart, BODIES};
use crate::chart_features::FeatureSet;
use std::f64::consts::TAU;

/// A point on the Earth's surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Site {
    /// Geodetic latitude, radians, north positive.
    pub lat: f64,
    /// Longitude, radians, east positive.
    pub lon: f64,
}

impl Site {
    pub fn from_degrees(lat_deg: f64, lon_deg: f64) -> Self {
        Site { lat: lat_deg.to_radians(), lon: lon_deg.to_radians() }
    }
}

fn norm_tau(x: f64) -> f64 {
    let r = x % TAU;
    if r < 0.0 { r + TAU } else { r }
}

/// Mean obliquity of the ecliptic at `day`, radians (IAU 1980).
pub fn obliquity_of_date(day: f64) -> f64 {
    let t = (day - 0.5) / 36525.0;
    (23.439_291_1 - 0.013_004_2 * t - 1.64e-7 * t * t + 5.036e-7 * t * t * t).to_radians()
}

/// Greenwich mean sidereal time at `day`, radians.
///
/// `day` is days since 2000-01-01T00:00 UTC, so the J2000 offset of half a day is
/// removed first. At 2000-01-01T12:00 this returns 18h 41m 50.5s, the defined value.
pub fn gmst(day: f64) -> f64 {
    let d = day - 0.5;
    let t = d / 36525.0;
    let deg = 280.460_618_37 + 360.985_647_366_29 * d + 0.000_387_933 * t * t
        - t * t * t / 38_710_000.0;
    norm_tau(deg.to_radians())
}

/// Local mean sidereal time at a site, radians.
pub fn lst(day: f64, site: Site) -> f64 {
    norm_tau(gmst(day) + site.lon)
}

/// Right ascension and declination referred to the equinox **of date**.
///
/// Takes the J2000 ecliptic longitude and latitude, advances the longitude by
/// general precession, and rotates by the obliquity of date.
pub fn equatorial_of_date(j2000_lon: f64, lat: f64, day: f64) -> (f64, f64) {
    let lambda = crate::chart::tropical_lon(j2000_lon, day);
    let eps = obliquity_of_date(day);
    let (sl, cl) = lambda.sin_cos();
    let (sb, cb) = lat.sin_cos();
    let (se, ce) = eps.sin_cos();
    let ra = norm_tau((sl * ce - (sb / cb) * se).atan2(cl));
    let dec = (sb * ce + cb * se * sl).asin();
    (ra, dec)
}

/// Ecliptic longitude of the midheaven, radians.
///
/// The MC is where the meridian cuts the ecliptic: `atan2(sin θ, cos θ · cos ε)`,
/// with θ the local sidereal time. Independent of latitude.
pub fn midheaven(theta: f64, eps: f64) -> f64 {
    let (st, ct) = theta.sin_cos();
    norm_tau(st.atan2(ct * eps.cos()))
}

/// Ecliptic longitude of the ascendant, radians.
///
/// Where the eastern horizon cuts the ecliptic. Unlike the MC this depends on
/// latitude, and it **degenerates near the poles**: above the polar circles the
/// ecliptic can fail to intersect the horizon in the usual way and the ascendant
/// jumps. Callers at extreme latitude should expect discontinuity, not an error.
pub fn ascendant(theta: f64, eps: f64, lat: f64) -> f64 {
    let (st, ct) = theta.sin_cos();
    let (se, ce) = eps.sin_cos();
    norm_tau(ct.atan2(-(st * ce + lat.tan() * se)))
}

/// Altitude and azimuth of a body, radians. Azimuth is measured east of north.
pub fn horizon(ha: f64, dec: f64, lat: f64) -> (f64, f64) {
    let (sh, ch) = ha.sin_cos();
    let (sd, cd) = dec.sin_cos();
    let (sp, cp) = lat.sin_cos();
    let alt = (sp * sd + cp * cd * ch).asin();
    let az = norm_tau((-sh * cd).atan2(sd * cp - cd * sp * ch));
    (alt, az)
}

/// Site-local features for one chart.
///
/// `max_harmonic` controls the aspect harmonics computed between each body and the
/// two angles; the hour-angle and horizon features are unaffected by it.
///
/// Names are unprefixed here — [`all`] adds the `local.` prefix.
pub fn angles(chart: &Chart, site: Site, max_harmonic: usize) -> FeatureSet {
    let mut f = FeatureSet::default();
    let theta = lst(chart.day, site);
    let eps = obliquity_of_date(chart.day);
    let mc = midheaven(theta, eps);
    let asc = ascendant(theta, eps, site.lat);

    f.push("last.cos", theta.cos());
    f.push("last.sin", theta.sin());
    f.push("mc.cos", mc.cos());
    f.push("mc.sin", mc.sin());
    f.push("asc.cos", asc.cos());
    f.push("asc.sin", asc.sin());
    // The ascendant-midheaven separation varies with latitude and season; it is a
    // compact summary of how obliquely the ecliptic meets this site's horizon.
    let sep = norm_tau(asc - mc);
    f.push("asc_mc.cos", sep.cos());
    f.push("asc_mc.sin", sep.sin());

    let mut n_above = 0.0;
    let mut alt_sum = 0.0;
    for (bi, body) in BODIES.iter().enumerate() {
        let s = &chart.states[bi];
        // A body at the observer's own location has no meaningful direction.
        if s.dist == 0.0 {
            continue;
        }
        let (ra, dec) = equatorial_of_date(s.lon, s.lat, chart.day);
        let ha = norm_tau(theta - ra);
        let (alt, az) = horizon(ha, dec, site.lat);

        // Hour angle as cos/sin, and at the semidiurnal harmonic: the tidal
        // potential's dominant terms are diurnal and semidiurnal in exactly this
        // angle, so both are given directly rather than left to be learned.
        f.push(format!("ha.{body}.cos"), ha.cos());
        f.push(format!("ha.{body}.sin"), ha.sin());
        f.push(format!("ha2.{body}.cos"), (2.0 * ha).cos());
        f.push(format!("ha2.{body}.sin"), (2.0 * ha).sin());
        f.push(format!("alt.{body}"), alt.sin());
        f.push(format!("az.{body}.cos"), az.cos());
        f.push(format!("az.{body}.sin"), az.sin());
        f.push(format!("above.{body}"), if alt > 0.0 { 1.0 } else { 0.0 });
        f.push(format!("dec_date.{body}"), dec.sin());

        if alt > 0.0 {
            n_above += 1.0;
        }
        alt_sum += alt.sin();

        // Aspects to the angles, by the same angle-sum recurrence used for
        // body-to-body aspects: one trig call, then n-1 multiply-adds.
        for (angle_name, angle) in [("asc", asc), ("mc", mc)] {
            let d = norm_tau(s.lon - angle);
            let (mut c, mut sn) = (1.0, 0.0);
            let (d_s, d_c) = d.sin_cos();
            for h in 1..=max_harmonic {
                let (nc, ns) = (c * d_c - sn * d_s, sn * d_c + c * d_s);
                c = nc;
                sn = ns;
                f.push(format!("{angle_name}.{body}.h{h}.cos"), c);
                f.push(format!("{angle_name}.{body}.h{h}.sin"), sn);
            }
        }
    }
    f.push("n_above", n_above);
    f.push("alt_sum", alt_sum);
    f
}

/// Site-local features with the `local.` prefix, ready to merge with the global set.
pub fn all(chart: &Chart, site: Site, max_harmonic: usize) -> FeatureSet {
    let inner = angles(chart, site, max_harmonic);
    let mut f = FeatureSet::default();
    for (n, v) in inner.names.iter().zip(&inner.values) {
        f.push(format!("local.{n}"), *v);
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chart::{BodyState, Frame};

    fn chart_with(day: f64, lon: f64, lat: f64) -> Chart {
        let mut states = vec![BodyState::default(); BODIES.len()];
        for s in states.iter_mut() {
            s.dist = 1.0;
            s.lon = lon;
            s.lat = lat;
        }
        Chart { day, frame: Frame::Geocentric, states }
    }

    #[test]
    fn gmst_matches_the_defined_value_at_j2000() {
        // GMST at 2000-01-01T12:00 UTC is 18h 41m 50.548s by definition.
        let h = gmst(0.5) / TAU * 24.0;
        let expect = 18.0 + 41.0 / 60.0 + 50.548 / 3600.0;
        assert!((h - expect).abs() < 1e-5, "gmst {h} h, expected {expect} h");
    }

    #[test]
    fn sidereal_day_is_shorter_than_a_solar_day() {
        // Sidereal time advances by a full turn plus ~1/365.25 per solar day, which
        // is what makes the sidereal day 3m 56s short.
        let a = gmst(1000.0);
        let b = gmst(1001.0);
        let advance = norm_tau(b - a);
        let excess_seconds = advance / TAU * 86400.0;
        assert!(
            (excess_seconds - 236.0).abs() < 1.0,
            "sidereal day runs {excess_seconds} s ahead, expected ~236 s"
        );
    }

    #[test]
    fn obliquity_decreases_and_is_correct_at_j2000() {
        let e0 = obliquity_of_date(0.5).to_degrees();
        assert!((e0 - 23.439_291).abs() < 1e-5, "obliquity {e0}");
        assert!(obliquity_of_date(36525.0) < obliquity_of_date(0.5));
    }

    #[test]
    fn a_body_on_the_meridian_is_at_the_midheaven() {
        // Put a body at the MC longitude and check its hour angle is zero: this
        // ties the MC formula to the hour-angle convention, catching a sign error
        // in either that a self-consistent test would miss.
        let day = 8000.3;
        let site = Site::from_degrees(0.0, 0.0);
        let theta = lst(day, site);
        let eps = obliquity_of_date(day);
        let mc = midheaven(theta, eps);
        // midheaven() returns a tropical longitude; invert to the J2000 value that
        // equatorial_of_date will precess back.
        let j2000_lon = mc - crate::chart::precession_since_j2000(day);
        let (ra, _) = equatorial_of_date(j2000_lon, 0.0, day);
        let ha = norm_tau(theta - ra);
        let off = if ha > std::f64::consts::PI { ha - TAU } else { ha };
        assert!(off.abs() < 1e-9, "hour angle at MC was {off} rad");
    }

    #[test]
    fn the_sun_is_highest_when_its_hour_angle_is_zero() {
        // Altitude must peak at transit for every latitude where the body rises.
        let dec = 0.2;
        for lat_deg in [-60.0, -20.0, 0.0, 35.0, 65.0] {
            let lat = (lat_deg as f64).to_radians();
            let (peak, _) = horizon(0.0, dec, lat);
            for k in 1..24 {
                let ha = k as f64 / 24.0 * TAU;
                let (alt, _) = horizon(ha, dec, lat);
                assert!(alt <= peak + 1e-12, "lat {lat_deg}: alt {alt} > transit {peak}");
            }
        }
    }

    #[test]
    fn altitude_at_transit_follows_the_standard_relation() {
        // At upper transit, altitude = 90 - |lat - dec|.
        let lat = 40.0_f64.to_radians();
        let dec = 10.0_f64.to_radians();
        let (alt, _) = horizon(0.0, dec, lat);
        let expect = std::f64::consts::FRAC_PI_2 - (lat - dec).abs();
        assert!((alt - expect).abs() < 1e-12, "{} vs {expect}", alt);
    }

    #[test]
    fn the_ascendant_is_rising_meaning_its_altitude_is_zero() {
        // The defining property: the ascendant is the ecliptic point on the eastern
        // horizon. Convert it to equatorial and its altitude must vanish.
        let day = 5000.75;
        for lat_deg in [-50.0, -10.0, 25.0, 55.0] {
            let site = Site::from_degrees(lat_deg, 12.0);
            let theta = lst(day, site);
            let eps = obliquity_of_date(day);
            let asc = ascendant(theta, eps, site.lat);
            let j2000_lon = asc - crate::chart::precession_since_j2000(day);
            let (ra, dec) = equatorial_of_date(j2000_lon, 0.0, day);
            let (alt, az) = horizon(norm_tau(theta - ra), dec, site.lat);
            assert!(alt.abs() < 1e-9, "lat {lat_deg}: ascendant altitude {alt}");
            // And it must be rising -- in the eastern half, azimuth in (0, pi).
            assert!(az > 0.0 && az < std::f64::consts::PI, "lat {lat_deg}: azimuth {az}");
        }
    }

    #[test]
    fn features_differ_between_sites_at_the_same_instant() {
        // The entire purpose of the module: two sites at one epoch must not share a
        // feature vector, or the model has no spatial channel.
        let c = chart_with(9000.0, 1.0, 0.05);
        let tokyo = all(&c, Site::from_degrees(35.7, 139.7), 4);
        let lima = all(&c, Site::from_degrees(-12.0, -77.0), 4);
        assert_eq!(tokyo.names, lima.names);
        let differing = tokyo
            .values
            .iter()
            .zip(&lima.values)
            .filter(|(a, b)| (*a - *b).abs() > 1e-6)
            .count();
        assert!(
            differing > tokyo.len() / 2,
            "only {differing} of {} features differ between Tokyo and Lima",
            tokyo.len()
        );
    }

    #[test]
    fn every_feature_is_finite_and_named_once() {
        let c = chart_with(7123.4, 2.5, -0.1);
        for lat in [-89.0, -45.0, 0.0, 45.0, 89.0] {
            let f = all(&c, Site::from_degrees(lat, 30.0), 6);
            for (n, v) in f.names.iter().zip(&f.values) {
                assert!(v.is_finite(), "{n} = {v} at latitude {lat}");
            }
            let unique: std::collections::HashSet<_> = f.names.iter().collect();
            assert_eq!(unique.len(), f.names.len(), "duplicate feature name");
        }
    }

    #[test]
    fn hour_angles_sweep_a_full_turn_across_a_day() {
        // A body's hour angle must advance by one turn per sidereal day. Checking it
        // through the feature surface confirms the wiring, not just the formula.
        let site = Site::from_degrees(20.0, 45.0);
        let mut prev: Option<f64> = None;
        let mut total = 0.0;
        for k in 0..=24 {
            let day = 6000.0 + k as f64 / 24.0;
            let c = chart_with(day, 0.0, 0.0);
            let f = all(&c, site, 1);
            let ha = f.get("local.ha.SUN.sin").unwrap().atan2(f.get("local.ha.SUN.cos").unwrap());
            if let Some(p) = prev {
                let mut d = ha - p;
                if d > std::f64::consts::PI { d -= TAU }
                if d < -std::f64::consts::PI { d += TAU }
                total += d;
            }
            prev = Some(ha);
        }
        assert!((total.abs() - TAU).abs() < 0.05, "swept {total} rad in a day");
    }
}
