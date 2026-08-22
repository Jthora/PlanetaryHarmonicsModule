//! End-to-end feature generation: ephemeris to the full feature vector.
//!
//!     cargo run --release --example feature_generation
//!
//! Reports the real column count across three frames, and times generation so the
//! cost of a full catalogue-scale run is known rather than guessed.

use ph_core::{chart, chart_cycles, chart_features, chart_local};
use rustspice_core::KernelSet;
use std::time::Instant;

const MAX_HARMONIC: usize = 24;
const MAX_BASE: usize = 24;
const LOCAL_HARMONIC: usize = 12;

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;

    // A day of hourly epochs, as a sample of the real grid.
    let days: Vec<f64> = (0..24).map(|i| 8766.0 + i as f64 / 24.0).collect();

    let frames = [
        chart::Frame::Geocentric,
        chart::Frame::Heliocentric,
        chart::Frame::Barycentric,
    ];

    let t0 = Instant::now();
    let mut charts = Vec::new();
    for f in frames {
        charts.push(chart::charts(&mut spice, &days, f, epoch2000)?);
    }
    let ephem = t0.elapsed();

    let t1 = Instant::now();
    let mut rows: Vec<chart_features::FeatureSet> = Vec::with_capacity(days.len());
    for i in 0..days.len() {
        let mut row = chart_features::FeatureSet::default();
        for frame_charts in &charts {
            row.extend(chart_features::all(&frame_charts[i], MAX_HARMONIC, MAX_BASE));
        }
        // Lunar cycles, eclipse geometry and fixed points are geocentric quantities.
        let cyc = chart_cycles::all(&charts[0][i], MAX_HARMONIC);
        for (n, v) in cyc.names.iter().zip(&cyc.values) {
            row.push(format!("geo.{n}"), *v);
        }
        rows.push(row);
    }
    // Station timing needs the whole series at once.
    let timing = chart_cycles::station_timing(&charts[0]);
    for (row, t) in rows.iter_mut().zip(timing) {
        for (n, v) in t.names.iter().zip(&t.values) {
            row.push(format!("geo.stn.{n}"), *v);
        }
    }
    let derive = t1.elapsed();

    // Site-local features, for a handful of cells spread over the globe.
    let sites = [
        ("tokyo", chart_local::Site::from_degrees(35.7, 139.7)),
        ("lima", chart_local::Site::from_degrees(-12.0, -77.0)),
        ("reykjavik", chart_local::Site::from_degrees(64.1, -21.9)),
    ];
    let t2 = Instant::now();
    let mut local_rows = Vec::new();
    for (_, site) in &sites {
        let mut per_site = Vec::with_capacity(days.len());
        for c in &charts[0] {
            per_site.push(chart_local::all(c, *site, LOCAL_HARMONIC));
        }
        local_rows.push(per_site);
    }
    let local = t2.elapsed();

    let n_local = local_rows[0][0].len();
    let n = rows[0].len();
    println!("bodies {}   frames {}   harmonics {MAX_HARMONIC}   bases {MAX_BASE}", chart::BODIES.len(), frames.len());
    println!("\nglobal features per epoch:     {n}");
    println!("site-local features per site: {n_local}");
    println!("total per (epoch, site):      {}", n + n_local);

    // Breakdown by family.
    let mut families: std::collections::BTreeMap<&str, usize> = Default::default();
    for name in &rows[0].names {
        let fam = name.splitn(3, '.').nth(1).unwrap_or("?");
        *families.entry(fam).or_insert(0) += 1;
    }
    for (fam, count) in &families {
        println!("  {fam:<6} {count:>6}");
    }

    println!("\ntiming for {} epochs:", days.len());
    println!("  ephemeris  {:>8.1} ms   ({:.2} ms/epoch)", ephem.as_secs_f64() * 1e3, ephem.as_secs_f64() * 1e3 / days.len() as f64);
    println!("  derivation {:>8.1} ms   ({:.2} ms/epoch)", derive.as_secs_f64() * 1e3, derive.as_secs_f64() * 1e3 / days.len() as f64);
    println!("  site-local {:>8.1} ms   ({:.2} ms/epoch/site)", local.as_secs_f64() * 1e3, local.as_secs_f64() * 1e3 / (days.len() * sites.len()) as f64);

    // Extrapolate to the real grid.
    let grid = 430_000.0;
    let per = (ephem + derive).as_secs_f64() / days.len() as f64;
    println!("\nfull hourly grid, 1976-2024 ({grid:.0} epochs):");
    println!("  generation      {:.1} minutes", per * grid / 60.0);
    println!("  materialised    {:.1} GB   (f64, global features only)", n as f64 * grid * 8.0 / 1e9);
    println!(
        "  primitives only {:.2} GB   ({} numbers/epoch across {} frames)",
        (chart::BODIES.len() * 8 * frames.len()) as f64 * grid * 8.0 / 1e9,
        chart::BODIES.len() * 8 * frames.len(),
        frames.len()
    );

    // Real-ephemeris check on the site-local layer: at Greenwich apparent noon the
    // Sun is on the meridian, so the midheaven must sit within a couple of degrees
    // of the Sun's tropical longitude (the equation of time is the difference).
    let noon = chart::charts(&mut spice, &[8766.5], chart::Frame::Geocentric, epoch2000)?;
    let greenwich = chart_local::Site::from_degrees(0.0, 0.0);
    let sun = noon[0].body("SUN").unwrap();
    let sun_lon = chart::tropical_lon(sun.lon, noon[0].day).to_degrees();
    let theta = chart_local::lst(noon[0].day, greenwich);
    let mc = chart_local::midheaven(theta, chart_local::obliquity_of_date(noon[0].day)).to_degrees();
    println!("\nGreenwich apparent noon: Sun {sun_lon:.3} deg, MC {mc:.3} deg, offset {:.3} deg", sun_lon - mc);

    // Spot-check a few values so the run proves the pipeline, not just its size.
    let r = &rows[12];
    for name in [
        "geo.shape.concentration",
        "geo.res.base12.mag",
        "helio.mot.n_retrograde",
        "geo.mot.n_retrograde",
    ] {
        println!("  {name} = {:.5}", r.get(name).unwrap());
    }
    println!();
    for (name, row) in sites.iter().map(|s| s.0).zip(&local_rows) {
        println!(
            "  {name:<10} n_above {:.0}   alt_sum {:+.4}   MC {:>7.2} deg",
            row[12].get("local.n_above").unwrap(),
            row[12].get("local.alt_sum").unwrap(),
            row[12].get("local.mc.sin").unwrap().atan2(row[12].get("local.mc.cos").unwrap()).to_degrees(),
        );
    }
    Ok(())
}
