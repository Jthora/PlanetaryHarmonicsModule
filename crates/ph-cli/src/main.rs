//! `ph-features` — export tidal features as CSV, with provenance.
//!
//! The interface between PlanetaryHarmonics and a downstream forecasting stack
//! (`docs/14`, `docs/15`). Rust computes the physics; Python consumes a table.
//!
//! Every export carries a provenance header. A feature whose reference frame or
//! epoch system is ambiguous is worthless in a statistical test, and across a repo
//! boundary that ambiguity is easy to introduce — so the header is not optional
//! and must not be stripped.

use ph_core::{doodson, fault, field::TidalField, love::Elastic};
use rustspice_core::{Et, KernelSet};
use std::fmt::Write as _;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const KERNELS: &[&str] = &["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"];

struct Args {
    kernels: String,
    lat: f64,
    lon: f64,
    strike: f64,
    dip: f64,
    rake: f64,
    mu: f64,
    start: String,
    days: f64,
    step: f64,
    constituents: Vec<String>,
    out: Option<String>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            kernels: "kernels".into(),
            lat: 35.635,
            lon: -120.150,
            strike: 137.0,
            dip: 90.0,
            rake: 180.0,
            mu: 0.4,
            start: "2001-01-01T00:00:00".into(),
            days: 365.0,
            step: 0.02,
            constituents: ["M2", "O1", "N2", "Q1", "Mf", "Msf", "Mm", "Ssa", "Sa"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            out: None,
        }
    }
}

fn usage() -> ! {
    eprintln!(
        "ph-features {VERSION} -- tidal feature export

  --kernels DIR       kernel directory            (default: kernels)
  --lat DEG           site latitude               (default: 35.635)
  --lon DEG           site longitude              (default: -120.150)
  --strike DEG        fault strike                (default: 137)
  --dip DEG           fault dip                   (default: 90)
  --rake DEG          fault rake                  (default: 180)
  --mu F              effective friction          (default: 0.4)
  --start ISO         first epoch                 (default: 2001-01-01T00:00:00)
  --days F            span in days                (default: 365)
  --step F            sample step in days         (default: 0.02)
  --constituents CSV  comma-separated names, or 'none'
  --out PATH          output file                 (default: stdout)

Defaults describe the deep San Andreas at Parkfield."
    );
    std::process::exit(2)
}

fn parse_args() -> Args {
    let mut a = Args::default();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let need = |i: usize| -> String {
            argv.get(i + 1).cloned().unwrap_or_else(|| usage())
        };
        let num = |i: usize| -> f64 { need(i).parse().unwrap_or_else(|_| usage()) };
        match argv[i].as_str() {
            "--kernels" => a.kernels = need(i),
            "--lat" => a.lat = num(i),
            "--lon" => a.lon = num(i),
            "--strike" => a.strike = num(i),
            "--dip" => a.dip = num(i),
            "--rake" => a.rake = num(i),
            "--mu" => a.mu = num(i),
            "--start" => a.start = need(i),
            "--days" => a.days = num(i),
            "--step" => a.step = num(i),
            "--out" => a.out = Some(need(i)),
            "--constituents" => {
                let v = need(i);
                a.constituents = if v == "none" {
                    Vec::new()
                } else {
                    v.split(',').map(|s| s.trim().to_string()).collect()
                };
            }
            "-h" | "--help" => usage(),
            _ => usage(),
        }
        i += 2;
    }
    a
}

fn main() -> rustspice_core::Result<()> {
    let a = parse_args();

    for name in &a.constituents {
        if doodson::constituent(name).is_none() {
            eprintln!("unknown constituent: {name}");
            std::process::exit(2);
        }
    }

    let mut ks = KernelSet::new();
    for k in KERNELS {
        ks.add_file(format!("{}/{k}", a.kernels))?;
    }
    let mut spice = ks.open()?;

    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let t_start = spice.parse_time(&a.start)?;
    let d0 = (t_start.0 - epoch2000.0) / 86400.0;
    let n = (a.days / a.step).ceil() as usize + 1;

    let days: Vec<f64> = (0..n).map(|i| d0 + i as f64 * a.step).collect();
    let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();

    let plane = fault::FaultPlane::new(a.strike, a.dip, a.rake);
    let elastic = Elastic::EARTH;
    let earth = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let tensors = earth.tensors(&mut spice, &epochs)?;

    let local: Vec<_> = tensors
        .iter()
        .map(|t| fault::to_local_ned(t, a.lat, a.lon))
        .collect();
    let cfs: Vec<f64> = local
        .iter()
        .map(|t| elastic.stress(fault::coulomb(t, &plane, a.mu)))
        .collect();

    let mut out = String::with_capacity(n * 120);
    writeln!(out, "# PlanetaryHarmonicsModule feature export").unwrap();
    writeln!(out, "# generator = ph-features {VERSION}").unwrap();
    writeln!(out, "# frame = IAU_EARTH (body-fixed), local North-East-Down at site").unwrap();
    writeln!(out, "# epoch_system = days since 2000-01-01T00:00 UTC").unwrap();
    writeln!(out, "# aberration = NONE (geometric; tidal force acts on instantaneous geometry)").unwrap();
    writeln!(out, "# kernels = {}", KERNELS.join(", ")).unwrap();
    writeln!(out, "# site_lat_deg = {}", a.lat).unwrap();
    writeln!(out, "# site_lon_deg = {}", a.lon).unwrap();
    writeln!(out, "# fault_strike_dip_rake_deg = {},{},{}", a.strike, a.dip, a.rake).unwrap();
    writeln!(out, "# effective_friction_mu = {}", a.mu).unwrap();
    writeln!(out, "# elastic = degree-2 Love h2={} l2={} shear={:.1e} Pa", elastic.h2, elastic.l2, elastic.shear_modulus).unwrap();
    writeln!(out, "# tier_A = t_days, cfs_pa, dcfs_dt_pa_per_day, tensor components").unwrap();
    writeln!(out, "# tier_B = phase_* (analytic Doodson arguments)").unwrap();
    writeln!(out, "# units = tensor s^-2, stress Pa, phase radians in [0,2pi)").unwrap();
    writeln!(out, "# caveat = stress calibration is a degree-2 scalar approximation, good to ~2x").unwrap();
    writeln!(out, "# caveat = dcfs_dt is a central difference on the sample grid, not an analytic derivative").unwrap();

    let mut header = String::from("t_days,cfs_pa,dcfs_dt_pa_per_day,t_nn,t_ee,t_dd,t_ne,t_nd,t_ed");
    for c in &a.constituents {
        write!(header, ",phase_{c}").unwrap();
    }
    writeln!(out, "{header}").unwrap();

    let cons: Vec<_> = a
        .constituents
        .iter()
        .map(|n| doodson::constituent(n).unwrap())
        .collect();

    for i in 0..n {
        // Central difference, one-sided at the ends.
        let d = if i == 0 {
            (cfs[1] - cfs[0]) / a.step
        } else if i == n - 1 {
            (cfs[n - 1] - cfs[n - 2]) / a.step
        } else {
            (cfs[i + 1] - cfs[i - 1]) / (2.0 * a.step)
        };
        let m = &local[i].m;
        write!(
            out,
            "{:.6},{:.4},{:.4},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
            days[i], cfs[i], d, m[0][0], m[1][1], m[2][2], m[0][1], m[0][2], m[1][2]
        )
        .unwrap();
        for c in &cons {
            write!(out, ",{:.6}", c.phase_at(days[i])).unwrap();
        }
        out.push('\n');
    }

    match a.out {
        Some(path) => {
            std::fs::write(&path, out).expect("write output");
            eprintln!("wrote {n} rows to {path}");
        }
        None => print!("{out}"),
    }
    Ok(())
}
