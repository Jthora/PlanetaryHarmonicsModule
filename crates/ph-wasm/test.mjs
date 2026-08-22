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

console.log(failures === 0 ? "\nall checks passed" : `\n${failures} FAILURES`);
process.exit(failures === 0 ? 0 : 1);
