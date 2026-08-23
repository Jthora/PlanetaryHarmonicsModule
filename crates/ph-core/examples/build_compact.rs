//! Fit a compact chart ephemeris from SPICE and check it against the source.
//!
//!     cargo run --release --example build_compact -- 1900 2100 kernels/chart.phce
//!
//! Reports the fitting error per body and then validates the finished file the way
//! a consumer will use it: by comparing whole charts against the SPICE path at
//! epochs that were never fitting nodes.

use ph_core::chart::{self, Frame, BODIES};
use ph_core::compact::{self, Center};
use rustspice_core::{Aberration, Et, KernelSet, Vec3};
use std::cell::RefCell;

fn days_for_year(y: i64) -> f64 {
    // Days from 2000-01-01 to y-01-01.
    let mut d = 0i64;
    if y >= 2000 {
        for yy in 2000..y {
            d += if (yy % 4 == 0 && yy % 100 != 0) || yy % 400 == 0 { 366 } else { 365 };
        }
    } else {
        for yy in y..2000 {
            d -= if (yy % 4 == 0 && yy % 100 != 0) || yy % 400 == 0 { 366 } else { 365 };
        }
    }
    d as f64
}

fn main() -> rustspice_core::Result<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let y0: i64 = a.first().map(|s| s.parse().unwrap()).unwrap_or(1900);
    let y1: i64 = a.get(1).map(|s| s.parse().unwrap()).unwrap_or(2100);
    let out = a.get(2).cloned().unwrap_or_else(|| "kernels/chart.phce".into());

    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let spice = RefCell::new(ks.open()?);
    let epoch2000 = spice.borrow_mut().parse_time("2000-01-01T00:00:00")?;

    let (d0, d1) = (days_for_year(y0), days_for_year(y1));
    println!("fitting {y0}-{y1}  (day {d0:.0} to {d1:.0}, {:.0} years)", (d1 - d0) / 365.25);

    let sample = |name: &str, center: Center, day: f64| -> Vec3 {
        let observer = match center {
            Center::Ssb => "SOLAR SYSTEM BARYCENTER",
            Center::Earth => "EARTH",
        };
        let et = Et(epoch2000.0 + day * 86400.0);
        spice
            .borrow_mut()
            .position(name, et, "ECLIPJ2000", observer, Aberration::None)
            .unwrap_or(Vec3 { x: 0.0, y: 0.0, z: 0.0 })
    };

    let t0 = std::time::Instant::now();
    let (eph, errors) = compact::build(d0, d1, sample);
    let bytes = eph.to_bytes();
    println!("built in {:.0}s\n", t0.elapsed().as_secs_f64());

    println!("{:<22} {:>10} {:>12}", "body", "max km", "max arcsec");
    for e in &errors {
        println!("{:<22} {:>10.3} {:>12.4}", e.body, e.max_km, e.max_arcsec);
    }
    let worst = errors.iter().fold(0.0f64, |m, e| m.max(e.max_arcsec));
    println!("\nworst angular error: {worst:.4} arcsec");

    std::fs::write(&out, &bytes).expect("write");
    let spice_size = std::fs::metadata("kernels/de440s.bsp").map(|m| m.len()).unwrap_or(0);
    println!(
        "wrote {out}: {:.2} MB   (de440s.bsp is {:.2} MB, {:.1}x larger)",
        bytes.len() as f64 / 1e6,
        spice_size as f64 / 1e6,
        spice_size as f64 / bytes.len() as f64
    );

    // Validate as a consumer would: whole charts, at epochs that were never nodes.
    println!("\nvalidating whole charts against SPICE at non-node epochs:");
    let days: Vec<f64> = (0..400).map(|k| d0 + (d1 - d0) * (k as f64 + 0.317) / 400.0).collect();
    for frame in [Frame::Geocentric, Frame::Heliocentric, Frame::Barycentric] {
        let mut s = spice.borrow_mut();
        let truth = chart::charts(&mut s, &days, frame, epoch2000)?;
        drop(s);
        let got = eph.charts(&days, frame);

        let mut worst_lon = 0.0f64;
        let mut worst_body = "";
        let mut worst_speed = 0.0f64;
        for (t, g) in truth.iter().zip(&got) {
            for (bi, (a, b)) in t.states.iter().zip(&g.states).enumerate() {
                if a.dist == 0.0 {
                    continue;
                }
                let mut d = (a.lon - b.lon).abs();
                if d > std::f64::consts::PI {
                    d = std::f64::consts::TAU - d;
                }
                let arcsec = d.to_degrees() * 3600.0;
                if arcsec > worst_lon {
                    worst_lon = arcsec;
                    worst_body = BODIES[bi];
                }
                let rel = (a.lon_speed - b.lon_speed).abs() / a.lon_speed.abs().max(1e-9);
                worst_speed = worst_speed.max(rel);
            }
        }
        println!(
            "  {frame:<14?} worst longitude {worst_lon:.4} arcsec ({worst_body}), \
             worst speed {:.2e} relative",
            worst_speed
        );
    }
    Ok(())
}
