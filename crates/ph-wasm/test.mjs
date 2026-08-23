// End-to-end check of the WASM surface, run under Node.
//   node crates/ph-wasm/test.mjs
import { readFileSync } from "fs";
import { Harmonics } from "./pkg/planetary_harmonics_wasm.js";

let failures = 0;
const check = (name, cond, detail = "") => {
  console.log(`${cond ? "  ok  " : "FAIL  "}${name}${detail ? "   " + detail : ""}`);
  if (!cond) failures++;
};

const h = new Harmonics();

// --- analytic: no kernels needed ---
const names = h.constituentNames();
check("constituent names", names.includes("M2") && names.length === 13, `${names.length} names`);
check("M2 period", Math.abs(h.constituentPeriod("M2") - 0.5175) < 1e-3);
check("Sa period", Math.abs(h.constituentPeriod("Sa") - 365.26) < 0.01);

// Phases must be uniform for times uniform in the argument.
const days = Float64Array.from({ length: 100000 }, (_, i) => i * 0.037);
const ph = h.constituentPhases("M2", days, -120.15);
const bins = new Array(12).fill(0);
for (const p of ph) bins[Math.floor((p / (2 * Math.PI)) * 12) % 12]++;
const worst = Math.max(...bins.map((b) => Math.abs(b - days.length / 12) / (days.length / 12)));
check("phase uniformity", worst < 0.05, `worst bin ${(worst * 100).toFixed(2)}%`);

// Commensurabilities.
const k = h.enumerateCommensurabilities(2, 2);
const rows = k.length / 2;
check("d'Alembert holds", [...Array(rows)].every((_, i) => k[2 * i] + k[2 * i + 1] === 0), `${rows} combos`);
// Jupiter-Saturn great conjunction from mean motions alone.
const per = h.commensurabilityPeriods(Int32Array.from([1, -1]), Float64Array.from([0.001450, 0.0005839]));
check("Jupiter-Saturn 19.86 yr", Math.abs(per[0] / 365.25 - 19.86) < 0.1, `${(per[0] / 365.25).toFixed(2)} yr`);

// Love-number stress conversion.
const pa = h.stressFromTensor(Float64Array.from([8.632e-14]));
check("M2 stress ~595 Pa", Math.abs(pa[0] - 595) < 10, `${pa[0].toFixed(0)} Pa`);

// --- ephemeris-backed ---
for (const f of ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"]) {
  h.loadKernel(f, readFileSync(`kernels/${f}`));
}
const t0 = h.parseTime("2024-01-01T00:00:00");
const t1 = h.parseTime("2024-07-01T00:00:00");
check("parseTime", Math.abs(t0 - 8766) < 1, `day ${t0.toFixed(1)}`);

// New moons: Moon-Sun elongation reaching zero, to 1 ms.
const nm = h.aspectTimes("MOON", "SUN", "EARTH", 0, t0, t1, 29.530588, 0.001);
const published = ["2024-01-11 11:57", "2024-02-09 22:59", "2024-03-10 09:00",
                   "2024-04-08 18:21", "2024-05-08 03:22", "2024-06-06 12:38"];
check("new moon count", nm.length === 6, `${nm.length} found`);
let maxDiff = 0;
for (let i = 0; i < Math.min(nm.length, published.length); i++) {
  const want = h.parseTime(published[i].replace(" ", "T") + ":00");
  maxDiff = Math.max(maxDiff, Math.abs(nm[i] - want) * 86400);
}
check("new moons match published", maxDiff < 35, `worst ${maxDiff.toFixed(1)} s`);

// Tidal tensors, batched.
const tt = h.tidalTensors(["MOON", "SUN"], Float64Array.from([8766, 8767]), "J2000", "EARTH");
check("tensor shape", tt.length === 12, `${tt.length} values`);
const trace = tt[0] + tt[1] + tt[2];
check("tensor is trace-free", Math.abs(trace) < 1e-20, `trace ${trace.toExponential(1)}`);

// --- charts and derived features ---
const bodies = h.bodyNames();
check("body names", bodies.length === 11 && bodies[0] === "SUN", `${bodies.length} bodies`);

// Chart primitives: eight numbers per body per epoch.
const cdays = Float64Array.from([8766, 8766.5, 8767]);
const prim = h.charts(cdays, "geocentric");
check("chart shape", prim.length === cdays.length * bodies.length * 8, `${prim.length} values`);

// The geocentric Sun on 2024-01-01 sits near 280 degrees tropical. The chart is
// J2000-referenced, so precession since J2000 must be added to compare.
const sunIdx = bodies.indexOf("SUN");
const sunLon = prim[sunIdx * 8];
const precession = ((50.2879 / 3600) * (8766 / 365.25) * Math.PI) / 180;
const tropical = (((sunLon + precession) * 180) / Math.PI + 360) % 360;
check("geocentric Sun longitude", Math.abs(tropical - 280.1) < 0.5, `${tropical.toFixed(2)} deg`);

// Earth is the observer geocentrically, so it is degenerate by construction.
const earthIdx = bodies.indexOf("EARTH");
check("observer is degenerate", prim[earthIdx * 8 + 2] === 0, `dist ${prim[earthIdx * 8 + 2]}`);

// Heliocentric Earth is the geocentric Sun turned around.
const helio = h.charts(Float64Array.from([8766]), "heliocentric");
const eLon = (helio[earthIdx * 8] * 180) / Math.PI;
const sLon = (prim[sunIdx * 8] * 180) / Math.PI;
check("helio Earth = geo Sun + 180", Math.abs(((eLon - sLon - 180 + 540) % 360) - 180) < 1e-4);

// Feature names need no kernels and define the column order.
const spec = [["geocentric", "heliocentric"], 8, 6, 6, true, 35.7, 139.7, 4];
const fnames = h.featureNames(...spec);
check("feature names", fnames.length > 1000, `${fnames.length} columns`);
check("names are unique", new Set(fnames).size === fnames.length);
check("names are stable", h.featureNames(...spec).join() === fnames.join());

// THE contract: the matrix width must equal the name count, or every column a
// consumer reads is mislabelled and nothing downstream can notice.
const fdays = Float64Array.from([8766, 8766.25, 8766.5, 8766.75]);
const mat = h.features(fdays, ...spec);
check("matrix width matches names",
      mat.length === fdays.length * fnames.length,
      `${mat.length} = ${fdays.length} x ${mat.length / fdays.length}`);
check("all features finite", mat.every(Number.isFinite));

// Site-local angles must actually vary with the site, or the spatial channel is
// dead and a consumer would silently get the same numbers everywhere.
const lima = h.features(fdays, ["geocentric", "heliocentric"], 8, 6, 6, true, -12.0, -77.0, 4);
let differing = 0;
for (let i = 0; i < fnames.length; i++) if (Math.abs(mat[i] - lima[i]) > 1e-6) differing++;
check("site changes the vector", differing > 50, `${differing} of ${fnames.length} columns differ`);

// A spec asking for site angles without a geocentric chart must fail loudly
// rather than return meaningless numbers.
let threw = false;
try { h.featureNames(["heliocentric"], 4, 4, 4, true, 0, 0, 4); } catch { threw = true; }
check("site without geocentric is rejected", threw);

let threwFrame = false;
try { h.featureNames(["ecliptic"], 4, 4, 4, false, 0, 0, 4); } catch { threwFrame = true; }
check("unknown frame is rejected", threwFrame);

// Lunar phase from the feature vector against the new moon found independently
// above: at a new moon the synodic cosine is +1.
const nmDay = nm[0];
const nmSpec = [["geocentric"], 2, 2, 2, false, 0, 0, 2];
const nmNames = h.featureNames(...nmSpec);
const synIdx = nmNames.indexOf("geo.moon.syn.h1.cos");
check("synodic feature present", synIdx >= 0);
const nmRow = h.features(Float64Array.from([nmDay]), ...nmSpec);
check("new moon has syn.h1.cos = +1", Math.abs(nmRow[synIdx] - 1) < 1e-3,
      `${nmRow[synIdx].toFixed(5)}`);

// --- compact chart ephemeris ---
// A second instance, fed only the 6 MB fit and the 5 KB leapsecond kernel, must
// agree with the 33 MB SPICE kernel it was derived from. Two instances rather
// than one, because the compact path takes precedence once loaded and the
// comparison would otherwise be against itself.
const c = new Harmonics();
c.loadKernel("naif0012.tls", readFileSync("kernels/naif0012.tls"));
const span = c.loadCompact(readFileSync("kernels/chart.phce"));
check("compact loaded", c.hasCompact() && span.length === 2,
      `days ${span[0].toFixed(0)} to ${span[1].toFixed(0)}`);

const compactBytes = readFileSync("kernels/chart.phce").length;
const spiceBytes = readFileSync("kernels/de440s.bsp").length;
check("compact is much smaller", compactBytes * 4 < spiceBytes,
      `${(compactBytes / 1e6).toFixed(2)} MB vs ${(spiceBytes / 1e6).toFixed(2)} MB`);

// Whole charts, at epochs chosen to fall between fitting nodes.
const cmpDays = Float64Array.from(
  { length: 60 }, (_, i) => -30000 + i * 1013.37 + 0.4271);
let worstArcsec = 0, worstBody = "";
for (const frame of ["geocentric", "heliocentric", "barycentric"]) {
  const ref = h.charts(cmpDays, frame);
  const got = c.charts(cmpDays, frame);
  check(`${frame} shape matches`, ref.length === got.length);
  for (let e = 0; e < cmpDays.length; e++) {
    for (let b = 0; b < bodies.length; b++) {
      const o = (e * bodies.length + b) * 8;
      if (ref[o + 2] === 0) continue; // observer is degenerate in both
      let d = Math.abs(ref[o] - got[o]);
      if (d > Math.PI) d = 2 * Math.PI - d;
      const arcsec = (d * 180 / Math.PI) * 3600;
      if (arcsec > worstArcsec) { worstArcsec = arcsec; worstBody = bodies[b]; }
    }
  }
}
check("compact matches SPICE below 0.1 arcsec", worstArcsec < 0.1,
      `worst ${worstArcsec.toFixed(4)} arcsec (${worstBody})`);

// Features from the compact path must match the SPICE path too, not just the
// primitives they are built from.
const fSpec = [["geocentric"], 6, 4, 4, false, 0, 0, 4];
const fd = Float64Array.from([8766.37, 9000.11]);
const fNames = h.featureNames(...fSpec);
const refF = h.features(fd, ...fSpec);
const gotF = c.features(fd, ...fSpec);
// Relative, not absolute. A few features carry raw kilometres -- Pluto's distance
// is 5.9e9 km, where f32 quantises to ~500 km -- so an absolute tolerance would
// be comparing the storage format against itself. Angular features are unaffected.
let worstRel = 0, worstName = "";
for (let i = 0; i < refF.length; i++) {
  const rel = Math.abs(refF[i] - gotF[i]) / Math.max(Math.abs(refF[i]), 1);
  if (rel > worstRel) { worstRel = rel; worstName = fNames[i % fNames.length]; }
}
// 1e-3 rather than something tighter: rate features are formed by cancellation
// -- dist_speed is (p.v)/r -- and pass through zero, so f32 noise in the inputs
// is amplified without bound in relative terms near the crossing. A genuine
// error here (wrong body, wrong frame, shifted column) would be order one, not
// order 1e-5.
check("features agree across sources", worstRel < 1e-3,
      `worst ${worstRel.toExponential(2)} on ${worstName}`);

// Outside the fitted span the answer must be refused, not extrapolated -- a
// Chebyshev fit diverges violently past its interval.
let outOfRange = false;
try { c.charts(Float64Array.from([500000]), "geocentric"); } catch { outOfRange = true; }
check("out-of-span day is rejected", outOfRange);

let badCompact = false;
try { new Harmonics().loadCompact(Uint8Array.from([1, 2, 3, 4, 5])); } catch { badCompact = true; }
check("malformed compact file is rejected", badCompact);

console.log(failures === 0 ? "\nall checks passed" : `\n${failures} FAILURES`);
process.exit(failures === 0 ? 0 : 1);
