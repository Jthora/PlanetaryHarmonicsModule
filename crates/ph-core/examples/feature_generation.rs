//! End-to-end feature generation: ephemeris to the full feature vector.
//!
//!     cargo run --release --example feature_generation
//!
//! Reports the real column count across three frames, and times generation so the
//! cost of a full catalogue-scale run is known rather than guessed.

use ph_core::{chart, chart_features};
use rustspice_core::KernelSet;
use std::time::Instant;

const MAX_HARMONIC: usize = 24;
const MAX_BASE: usize = 24;

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
        rows.push(row);
    }
    let derive = t1.elapsed();

    let n = rows[0].len();
    println!("bodies {}   frames {}   harmonics {MAX_HARMONIC}   bases {MAX_BASE}", chart::BODIES.len(), frames.len());
    println!("\nfeatures per epoch: {n}");

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

    // Extrapolate to the real grid.
    let grid = 430_000.0;
    let per = (ephem + derive).as_secs_f64() / days.len() as f64;
    println!("\nfull hourly grid, 1976-2024 ({grid:.0} epochs):");
    println!("  generation      {:.1} minutes", per * grid / 60.0);
    println!("  materialised    {:.1} GB   (f64, all features)", n as f64 * grid * 8.0 / 1e9);
    println!(
        "  primitives only {:.2} GB   ({} numbers/epoch across {} frames)",
        (chart::BODIES.len() * 8 * frames.len()) as f64 * grid * 8.0 / 1e9,
        chart::BODIES.len() * 8 * frames.len(),
        frames.len()
    );

    // Spot-check a few values so the run proves the pipeline, not just its size.
    let r = &rows[12];
    for name in [
        "geo.shape.concentration",
        "geo.res.base12.mag",
        "helio.mot.n_retrograde",
        "geo.mot.n_retrograde",
    ] {
        println!("\n  {name} = {:.5}", r.get(name).unwrap());
    }
    Ok(())
}
