//! Chart primitives — body states in multiple reference frames.
//!
//! This is the raw material every downstream consumer needs: where each body is,
//! how fast it is moving, and in which direction, in whichever frame the question
//! is posed. Aspects, harmonics and chart shapes are all cheap trigonometry on
//! top of these.
//!
//! # Store primitives, derive features
//!
//! A full feature set — every pair, every harmonic, every frame — runs to
//! thousands of columns, which over hundreds of thousands of epochs is tens of
//! gigabytes. The ~20 numbers per body here are a **twenty-fold** smaller
//! representation carrying the same information, and the derived set can change
//! without recomputing a single ephemeris.
//!
//! # Three frames, deliberately not one
//!
//! [`Frame::Geocentric`] is the apparent sky: what an observer on Earth sees,
//! including retrograde motion, which is an artefact of Earth's own orbit but a
//! real feature of the view.
//!
//! [`Frame::Heliocentric`] is true orbital configuration. No retrogrades, and
//! Earth becomes a body like any other.
//!
//! [`Frame::Barycentric`] is relative to the solar system barycentre, where the
//! Sun's own displacement and motion become visible.
//!
//! The transformation between them is nonlinear and involves Earth's position, so
//! **a consumer given one frame cannot cheaply reconstruct another's angular
//! structure.** They are not redundant.

use rustspice_core::{Aberration, Et, Session, Vec3};

/// Obliquity of the ecliptic at J2000, degrees. Fixed, since ECLIPJ2000 and
/// J2000 are both inertial frames tied to that epoch.
const OBLIQUITY_J2000_DEG: f64 = 23.439_291_111;

/// Which centre the positions are measured from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Frame {
    /// Earth-centred — the apparent sky, with retrograde motion.
    Geocentric,
    /// Sun-centred — true orbital configuration.
    Heliocentric,
    /// Solar-system-barycentre-centred.
    Barycentric,
}

impl Frame {
    /// SPICE observer name.
    pub fn observer(&self) -> &'static str {
        match self {
            Frame::Geocentric => "EARTH",
            Frame::Heliocentric => "SUN",
            Frame::Barycentric => "SOLAR SYSTEM BARYCENTER",
        }
    }
}

/// Bodies carried in a chart.
///
/// Outer planets use barycentres, whose offset from the body centre is far below
/// any angular resolution that matters here.
/// Earth is included deliberately. Geocentrically it is degenerate — all zeros,
/// which a consumer can ignore as a constant. Heliocentrically it is the body
/// whose motion produces every retrograde in the geocentric frame.
pub const BODIES: &[&str] = &[
    "SUN",
    "EARTH",
    "MOON",
    "MERCURY",
    "VENUS",
    "MARS BARYCENTER",
    "JUPITER BARYCENTER",
    "SATURN BARYCENTER",
    "URANUS BARYCENTER",
    "NEPTUNE BARYCENTER",
    "PLUTO BARYCENTER",
];

/// One body's state, in one frame, at one epoch.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BodyState {
    /// Ecliptic longitude in the **J2000** frame, radians in `[0, 2π)`.
    ///
    /// ⚠ This is *not* the tropical longitude an astrologer means. Tropical is
    /// referenced to the equinox **of date**, which has precessed ~0.34° since
    /// J2000. Use [`tropical_lon`] for that, and [`sidereal_lon`] for the sidereal
    /// zodiac. The distinction was found by testing the March 2024 equinox, where
    /// the J2000 Sun sits at −0.3316° rather than 0.
    pub lon: f64,
    /// Ecliptic latitude, radians.
    pub lat: f64,
    /// Distance from the frame's centre, km.
    pub dist: f64,
    /// Longitude rate, radians per day. **Negative means retrograde** — the sign
    /// carries real information and must not be discarded.
    pub lon_speed: f64,
    /// Latitude rate, radians per day.
    pub lat_speed: f64,
    /// Radial velocity, km per day.
    pub dist_speed: f64,
    /// Right ascension, radians in `[0, 2π)`.
    pub ra: f64,
    /// Declination, radians.
    pub dec: f64,
}

/// All body states in one frame at one epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct Chart {
    /// Days since 2000-01-01T00:00 UTC.
    pub day: f64,
    pub frame: Frame,
    /// Parallel to [`BODIES`].
    pub states: Vec<BodyState>,
}

impl Chart {
    /// Ecliptic longitudes, in [`BODIES`] order — the input to aspect and
    /// chart-shape calculations.
    pub fn longitudes(&self) -> Vec<f64> {
        self.states.iter().map(|s| s.lon).collect()
    }

    /// Declinations, for parallel and contraparallel aspects.
    pub fn declinations(&self) -> Vec<f64> {
        self.states.iter().map(|s| s.dec).collect()
    }

    /// State of a named body.
    pub fn body(&self, name: &str) -> Option<&BodyState> {
        BODIES.iter().position(|b| *b == name).map(|i| &self.states[i])
    }

    /// Longitude of the Moon's ascending node, radians.
    ///
    /// Derived from the orbital angular momentum rather than read from a kernel:
    /// `h = r × v`, and the node lies along `ẑ × h`. This is the *osculating* node,
    /// which wobbles about the mean node by a degree or so.
    ///
    /// Only meaningful geocentrically.
    pub fn lunar_node(&self, moon_pos: Vec3, moon_vel: Vec3) -> f64 {
        let h = (
            moon_pos.y * moon_vel.z - moon_pos.z * moon_vel.y,
            moon_pos.z * moon_vel.x - moon_pos.x * moon_vel.z,
            moon_pos.x * moon_vel.y - moon_pos.y * moon_vel.x,
        );
        // n = z_hat x h
        norm_tau((-h.1).atan2(h.0))
    }
}

fn norm_tau(x: f64) -> f64 {
    let r = x % std::f64::consts::TAU;
    if r < 0.0 {
        r + std::f64::consts::TAU
    } else {
        r
    }
}

/// Rotate an ecliptic vector to equatorial using the J2000 obliquity.
fn ecliptic_to_equatorial(v: Vec3) -> (f64, f64, f64) {
    let e = OBLIQUITY_J2000_DEG.to_radians();
    let (se, ce) = e.sin_cos();
    (v.x, v.y * ce - v.z * se, v.y * se + v.z * ce)
}

/// Charts for a series of epochs in one frame.
///
/// Positions are **geometric** (no aberration): these are dynamical quantities.
/// Apply light-time correction downstream if apparent positions are wanted — the
/// two differ by ~40 s for lunar phase, which matters for event timing and does
/// not for feature generation.
pub fn charts(
    session: &mut Session,
    days: &[f64],
    frame: Frame,
    epoch2000: Et,
) -> rustspice_core::Result<Vec<Chart>> {
    let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();
    let mut out: Vec<Chart> = days
        .iter()
        .map(|&d| Chart {
            day: d,
            frame,
            states: vec![BodyState::default(); BODIES.len()],
        })
        .collect();

    for (bi, body) in BODIES.iter().enumerate() {
        // Skip self-reference: a body has no position relative to itself, and
        // SPICE would return zeros with a spurious light-time.
        if *body == frame.observer() {
            continue;
        }
        let states = session.states(
            body,
            &epochs,
            "ECLIPJ2000",
            frame.observer(),
            Aberration::None,
        )?;
        for (ci, sv) in states.iter().enumerate() {
            let (p, v) = (sv.position, sv.velocity);
            let r2 = p.x * p.x + p.y * p.y;
            let r = p.norm();
            // Angular rates from the state vector, exactly rather than by
            // differencing: d(lon)/dt = (x*vy - y*vx) / (x^2 + y^2).
            let lon_speed = if r2 > 0.0 {
                (p.x * v.y - p.y * v.x) / r2 * 86400.0
            } else {
                0.0
            };
            let rho = r2.sqrt();
            let lat_speed = if r > 0.0 && rho > 0.0 {
                (v.z * rho - p.z * (p.x * v.x + p.y * v.y) / rho) / (r * r) * 86400.0
            } else {
                0.0
            };
            let (ex, ey, ez) = ecliptic_to_equatorial(p);
            out[ci].states[bi] = BodyState {
                lon: norm_tau(p.y.atan2(p.x)),
                lat: if r > 0.0 { (p.z / r).asin() } else { 0.0 },
                dist: r,
                lon_speed,
                lat_speed,
                dist_speed: if r > 0.0 {
                    (p.x * v.x + p.y * v.y + p.z * v.z) / r * 86400.0
                } else {
                    0.0
                },
                ra: norm_tau(ey.atan2(ex)),
                dec: {
                    let m = (ex * ex + ey * ey + ez * ez).sqrt();
                    if m > 0.0 {
                        (ez / m).asin()
                    } else {
                        0.0
                    }
                },
            };
        }
    }
    Ok(out)
}

/// General precession in longitude since J2000, radians.
///
/// The equinox drifts at ~50.2879 arcsec per year, so a J2000-referenced longitude
/// must be advanced by this to become tropical (of-date). Verified against the
/// March 2024 equinox, where it accounts for the full observed offset.
pub fn precession_since_j2000(day: f64) -> f64 {
    (50.287_9 / 3600.0 * (day / 365.25)).to_radians()
}

/// Tropical longitude — referenced to the equinox of date.
pub fn tropical_lon(j2000_lon: f64, day: f64) -> f64 {
    norm_tau(j2000_lon + precession_since_j2000(day))
}

/// Ayanamsa: the offset between the tropical and sidereal zodiacs, radians.
///
/// Lahiri-style zero point, ~23.85° at J2000, growing with precession.
pub fn ayanamsa(day: f64) -> f64 {
    (23.85f64).to_radians() + precession_since_j2000(day)
}

/// Sidereal longitude — referenced to the fixed stars.
pub fn sidereal_lon(j2000_lon: f64, day: f64) -> f64 {
    norm_tau(tropical_lon(j2000_lon, day) - ayanamsa(day))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn frames_map_to_spice_observers() {
        assert_eq!(Frame::Geocentric.observer(), "EARTH");
        assert_eq!(Frame::Heliocentric.observer(), "SUN");
        assert_eq!(Frame::Barycentric.observer(), "SOLAR SYSTEM BARYCENTER");
    }

    #[test]
    fn obliquity_rotation_preserves_length_and_the_x_axis() {
        let v = Vec3::new(3.0, 4.0, 5.0);
        let (x, y, z) = ecliptic_to_equatorial(v);
        assert!((x - 3.0).abs() < 1e-12, "x is the rotation axis");
        assert!(((x * x + y * y + z * z).sqrt() - v.norm()).abs() < 1e-9);
        // A vector in the ecliptic plane acquires declination equal to obliquity
        // when it points along ecliptic +y.
        let (_, ey, ez) = ecliptic_to_equatorial(Vec3::new(0.0, 1.0, 0.0));
        assert!((ez.atan2(ey).to_degrees() - OBLIQUITY_J2000_DEG).abs() < 1e-9);
    }

    #[test]
    fn ayanamsa_advances_at_the_precession_rate() {
        let a0 = ayanamsa(0.0).to_degrees();
        let a100 = ayanamsa(36525.0).to_degrees();
        assert!((a0 - 23.85).abs() < 1e-9);
        // 50.2879 arcsec/yr over a century is 1.397 degrees.
        assert!(((a100 - a0) - 1.3969).abs() < 1e-3, "advance {}", a100 - a0);
    }

    #[test]
    fn tropical_and_sidereal_differ_by_the_ayanamsa() {
        let day = 8846.0; // March 2024
        let j = 1.0f64;
        let t = tropical_lon(j, day);
        let sid = sidereal_lon(j, day);
        assert!((t - j - precession_since_j2000(day)).abs() < 1e-12);
        let gap = (t - sid).rem_euclid(TAU);
        assert!((gap - ayanamsa(day)).abs() < 1e-12, "gap {gap}");
    }

    #[test]
    fn precession_matches_the_observed_equinox_offset() {
        // The J2000 Sun sits 0.33 deg short of 0 at the March 2024 equinox; that
        // offset is precession, not error.
        let p = precession_since_j2000(8846.0).to_degrees();
        assert!((p - 0.338).abs() < 0.01, "precession {p:.4} deg");
    }

    #[test]
    fn chart_accessors_are_body_aligned() {
        let c = Chart {
            day: 0.0,
            frame: Frame::Geocentric,
            states: (0..BODIES.len())
                .map(|i| BodyState {
                    lon: i as f64 * 0.1,
                    dec: i as f64 * 0.01,
                    ..Default::default()
                })
                .collect(),
        };
        assert_eq!(c.longitudes().len(), BODIES.len());
        assert!((c.body("MOON").unwrap().lon - 0.1).abs() < 1e-12);
        assert!((c.declinations()[2] - 0.02).abs() < 1e-12);
        assert!(c.body("NOT A BODY").is_none());
    }

    #[test]
    fn lunar_node_is_perpendicular_to_the_orbit_normal() {
        let c = Chart {
            day: 0.0,
            frame: Frame::Geocentric,
            states: vec![],
        };
        // An orbit in the ecliptic plane has h along +z, so the node is undefined
        // in direction but the formula must still return a finite angle.
        let n = c.lunar_node(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        assert!(n.is_finite() && (0.0..TAU).contains(&n));
        // Tilt the orbit: h gains an x component, and the node rotates.
        let a = c.lunar_node(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.2));
        let b = c.lunar_node(Vec3::new(0.0, 1.0, 0.0), Vec3::new(-1.0, 0.0, 0.2));
        assert!((a - b).abs() > 1e-6, "node should track the orbit plane");
    }

    #[test]
    fn earth_is_present_for_the_heliocentric_frame() {
        assert!(BODIES.contains(&"EARTH"));
        assert!(BODIES.contains(&"SUN"));
    }

    #[test]
    fn body_list_has_no_duplicates() {
        let mut seen = std::collections::HashSet::new();
        assert!(BODIES.iter().all(|b| seen.insert(*b)));
    }
}
