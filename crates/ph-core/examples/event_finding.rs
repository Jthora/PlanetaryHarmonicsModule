//! Angle-domain event finding against published new moon times.
//!
//!     cargo run --release --example event_finding
//!
//! New moon is the instant of zero geocentric ecliptic longitude difference
//! between Moon and Sun — a real angular event with published times, so it
//! validates [`ph_core::events`] against something external rather than against
//! more of our own code.
//!
//! It is also the shape of what Star Seer needs: exact times of angular
//! configurations, to arbitrary precision, over long spans.
//!
//! # Aberration: apparent, not geometric
//!
//! Published new moon times use **apparent** geocentric longitudes, so this uses
//! `LightTimeStellar`. That is the opposite of what tidal work needs
//! (`docs/10-rustspice-requirements.md` §2 argues for geometric positions, since
//! tidal force acts on the instantaneous configuration).
//!
//! The difference is not academic here. Solar aberration is ~20.5 arcsec and
//! elongation advances at ~30.5 arcsec per minute, so using geometric positions
//! shifts every new moon **~40 s late** — a systematic bias larger than the
//! minute-rounding of the published values.

use ph_core::events;
use rustspice_core::{Aberration, Et, KernelSet};
use std::cell::Cell;
use std::f64::consts::TAU;

/// Published new moons, 2024 (UTC). Astronomical Almanac values.
const KNOWN: &[(&str, &str)] = &[
    ("2024-01-11", "11:57"),
    ("2024-02-09", "22:59"),
    ("2024-03-10", "09:00"),
    ("2024-04-08", "18:21"),
    ("2024-05-08", "03:22"),
    ("2024-06-06", "12:38"),
];

const SYNODIC_DAYS: f64 = 29.530588;

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let t0 = spice.parse_time("2024-01-01T00:00:00")?;
    let t1 = spice.parse_time("2024-07-01T00:00:00")?;

    // Elongation: Moon's ecliptic longitude minus the Sun's, unwrapped by adding
    // the mean motion back so the angle is monotonic and the linear predictor works.
    let calls = Cell::new(0u32);
    let rate = TAU / SYNODIC_DAYS;
    let spice_cell = std::cell::RefCell::new(spice);
    let elong = |days: f64| -> f64 {
        calls.set(calls.get() + 1);
        let mut s = spice_cell.borrow_mut();
        let et = Et(t0.0 + days * 86400.0);
        let m = s
            .position("MOON", et, "ECLIPJ2000", "EARTH", Aberration::LightTimeStellar)
            .unwrap();
        let u = s
            .position("SUN", et, "ECLIPJ2000", "EARTH", Aberration::LightTimeStellar)
            .unwrap();
        let d = m.y.atan2(m.x) - u.y.atan2(u.x);
        // Re-add the mean advance so theta grows monotonically.
        let wrapped = d.rem_euclid(TAU);
        wrapped + TAU * (days / SYNODIC_DAYS - (wrapped / TAU)).round()
    };

    let span = (t1.0 - t0.0) / 86400.0;
    // One millisecond of elongation, in radians.
    let tol = rate * (0.001 / 86400.0);
    let found = events::crossings(&elong, rate, 0.0, 0.0, span, tol);
    let used = calls.get();

    println!("span {span:.0} days, tolerance {:.1e} rad (~1 ms)\n", tol);
    println!(
        "{:<22} {:<22} {:>9} {:>6}",
        "computed (UTC)", "published", "diff (s)", "iters"
    );
    println!("{}", "-".repeat(62));

    let mut worst: f64 = 0.0;
    let mut spice = spice_cell.into_inner();
    for (c, (d, hm)) in found.iter().zip(KNOWN) {
        let et = Et(t0.0 + c.time * 86400.0);
        let got = spice.format_utc(et, "ISOC", 3)?;
        let want = spice.parse_time(&format!("{d}T{hm}:00"))?;
        let diff = et.0 - want.0;
        worst = worst.max(diff.abs());
        println!("{:<22} {:<13} {:>9.1} {:>15}", got, format!("{d} {hm}"), diff, c.iterations);
    }

    println!(
        "\nlargest disagreement: {:.1} s   (published times are given to the minute)",
        worst
    );
    println!("ephemeris evaluations: {used} for {} events", found.len());
    println!(
        "  = {:.1} per event; a 1 ms scan of this span would need {:.0}",
        used as f64 / found.len() as f64,
        span * 86400.0 * 1000.0
    );
    Ok(())
}
