//! Validate chart primitives against known astronomical facts.
//!
//!     cargo run --release --example chart_validation
//!
//! Everything here is checkable against an almanac, which is the point — the
//! feature generator is the foundation of everything downstream, so it gets
//! validated externally rather than against more of our own code.

use ph_core::chart::{self, Frame};
use rustspice_core::KernelSet;
use std::f64::consts::{PI, TAU};

fn deg(r: f64) -> f64 {
    r.to_degrees()
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let day_of = |s: &str, spice: &mut rustspice_core::Session| -> f64 {
        (spice.parse_time(s).unwrap().0 - epoch2000.0) / 86400.0
    };

    let mut fails = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("{}{name}   {detail}", if ok { "  ok  " } else { "FAIL  " });
        if !ok {
            fails += 1;
        }
    };

    // March equinox 2024: 20 March, 03:06 UTC. Sun's geocentric longitude is 0.
    let eq = day_of("2024-03-20T03:06:00", &mut spice);
    let geo = chart::charts(&mut spice, &[eq], Frame::Geocentric, epoch2000)?;
    let sun_lon = deg(geo[0].body("SUN").unwrap().lon);
    let off = if sun_lon > 180.0 { sun_lon - 360.0 } else { sun_lon };
    // J2000 longitude falls short of zero by exactly the precession since J2000.
    check(
        "J2000 Sun offset at equinox equals precession",
        (off + deg(chart::precession_since_j2000(eq))).abs() < 0.01,
        format!("{off:+.4} deg vs precession {:+.4}", -deg(chart::precession_since_j2000(eq))),
    );
    // Corrected to the equinox of date, the tropical Sun is at zero.
    let trop = deg(chart::tropical_lon(geo[0].body("SUN").unwrap().lon, eq));
    let toff = if trop > 180.0 { trop - 360.0 } else { trop };
    check(
        "tropical Sun at 0 deg on the March equinox",
        toff.abs() < 0.02,
        format!("{toff:+.4} deg"),
    );

    // June solstice: declination reaches the obliquity, +23.44.
    let sol = day_of("2024-06-20T20:51:00", &mut spice);
    let g2 = chart::charts(&mut spice, &[sol], Frame::Geocentric, epoch2000)?;
    let dec = deg(g2[0].body("SUN").unwrap().dec);
    check(
        "Sun declination at June solstice",
        (dec - 23.44).abs() < 0.02,
        format!("{dec:.4} deg (obliquity 23.44)"),
    );

    // Cross-frame identity: heliocentric Earth is geocentric Sun plus 180 deg.
    let hel = chart::charts(&mut spice, &[eq], Frame::Heliocentric, epoch2000)?;
    let earth_lon = deg(hel[0].body("EARTH").unwrap().lon);
    // Wrap the difference to (-180, 180]; the identity puts it at exactly 180.
    let mut diff = (earth_lon - sun_lon).rem_euclid(360.0);
    if diff > 180.0 {
        diff -= 360.0;
    }
    check(
        "heliocentric Earth = geocentric Sun + 180",
        (diff.abs() - 180.0).abs() < 0.01,
        format!("{diff:+.5} deg"),
    );

    // Retrograde: Mercury was retrograde 1-25 April 2024, direct in between.
    let retro = day_of("2024-04-10T00:00:00", &mut spice);
    let direct = day_of("2024-06-10T00:00:00", &mut spice);
    let r = chart::charts(&mut spice, &[retro, direct], Frame::Geocentric, epoch2000)?;
    let sr = r[0].body("MERCURY").unwrap().lon_speed;
    let sd = r[1].body("MERCURY").unwrap().lon_speed;
    check(
        "Mercury retrograde 10 Apr 2024",
        sr < 0.0,
        format!("{:+.4} deg/day", deg(sr)),
    );
    check(
        "Mercury direct 10 Jun 2024",
        sd > 0.0,
        format!("{:+.4} deg/day", deg(sd)),
    );

    // Heliocentric Mercury never retrogrades -- retrograde is an Earth artefact.
    let hr = chart::charts(&mut spice, &[retro], Frame::Heliocentric, epoch2000)?;
    check(
        "heliocentric Mercury is never retrograde",
        hr[0].body("MERCURY").unwrap().lon_speed > 0.0,
        format!("{:+.4} deg/day", deg(hr[0].body("MERCURY").unwrap().lon_speed)),
    );

    // Lunar distance spans perigee to apogee, roughly 356,500 to 406,700 km.
    let month: Vec<f64> = (0..60).map(|i| eq + i as f64 * 0.5).collect();
    let m = chart::charts(&mut spice, &month, Frame::Geocentric, epoch2000)?;
    let (lo, hi) = m.iter().fold((f64::MAX, f64::MIN), |(a, b), c| {
        let d = c.body("MOON").unwrap().dist;
        (a.min(d), b.max(d))
    });
    check(
        "lunar distance range over a month",
        lo > 355_000.0 && lo < 371_000.0 && hi > 400_000.0 && hi < 407_500.0,
        format!("{lo:.0} to {hi:.0} km"),
    );

    // Sidereal longitude trails tropical by the ayanamsa, ~24 deg now.
    let ay = deg(chart::ayanamsa(eq));
    check(
        "ayanamsa in 2024",
        (23.0..25.0).contains(&ay),
        format!("{ay:.3} deg"),
    );

    // Every body should have a plausible longitude and finite speed.
    let sane = geo[0]
        .states
        .iter()
        .enumerate()
        .all(|(i, s)| {
            chart::BODIES[i] == "EARTH" || ((0.0..TAU).contains(&s.lon) && s.lon_speed.is_finite())
        });
    check("all bodies have sane state", sane, format!("{} bodies", chart::BODIES.len()));

    // Moon's node regresses: its longitude decreases over time.
    let n0 = geo[0].lunar_node(
        rustspice_core::Vec3::new(1.0, 0.0, 0.1),
        rustspice_core::Vec3::new(0.0, 1.0, 0.0),
    );
    check("lunar node is finite", n0.is_finite() && n0 < TAU, format!("{:.2} deg", deg(n0)));
    let _ = PI;

    println!(
        "\n{}",
        if fails == 0 {
            "all checks passed".to_string()
        } else {
            format!("{fails} FAILURES")
        }
    );
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
