# EarthquakeForecastModule — Handoff

**Copy this file into the new repository as `HANDOFF.md`.** It is written to be
read cold, without access to the PlanetaryHarmonics research log.

Source of record for everything below: the `docs/` tree of
[PlanetaryHarmonicsModule](https://github.com/Jthora/PlanetaryHarmonicsModule).

---

## 1. What this repository is

A probabilistic earthquake forecasting system that tests whether celestial, tidal,
gravimetric and geophysical features improve earthquake-rate estimates over
established baselines.

**What it is not:** an earthquake prediction system. In seismology, *prediction*
denotes deterministic time/place/magnitude claims and is largely discredited.
*Forecasting* means probabilistic rate estimation — what CSEP evaluates. Use
"forecast" everywhere, in code, docs and outputs.

**It does not depend on AstrologyCore.** The dependency chain forks:

```text
RustSPICE
  └─> PlanetaryHarmonicsModule
        ├─> AstrologyCore ──> Cosmic Cypher, Star Seer, Resonant Finder
        └─> EarthquakeForecastModule          ← you are here
```

This is deliberate. Depending on the interpretive layer would inherit its baggage
and reviewers would be right to discount the work. **The dependency graph is part
of the argument — do not add that edge.**

---

## 2. The scientific state of play

Established from primary sources. Numbers matter; do not round them into folklore.

### Effect sizes are real, small, and setting-dependent

| Setting | Effect | Source |
|---|---|---|
| Shallow thrust / subduction | **rate ×3**, p < 10⁻⁴ | Cochran, Vidale & Tanaka 2004, *Science* 306 |
| Global (442k events) | ~99% confidence, uplift phase | Métivier et al. 2009, *EPSL* 278 |
| Southern California | **null** | *GJI* 205, 681 |
| Tectonic tremor | very strong, 12.4 h and 24–25 h pulsing | Rubinstein et al. 2008, *Science* 319 |

Neither "tides cause earthquakes" nor "tides do nothing" is supportable.

### Sample size — you need ~10⁴ events, not 10⁵–10⁶

Beeler & Lockner 2003 (*JGR* 108(B8), **free from USGS**), equation 18:

```
N ≥ ln(P_rw) / (Δτ_u / (2 a σ_n))²
```

Worked examples: **6.2×10³–5.5×10⁴** events for Δτ = 0.01 MPa; **>13,000** for daily
Earth tides. N scales as the **inverse square** of normalised stress amplitude — a
10× amplitude increase cuts the requirement 100×.

*(An earlier draft of our notes recorded 10⁵–10⁶ from an unverified search summary.
It was wrong by two orders of magnitude. Regional subsets are viable.)*

### The response is band-limited — the key operating hypothesis ★

Two timescales bound the responsive band:

```
t_n   nucleation duration      damps response above 1/t_n
      Beeler & Lockner extrapolate t_n ≥ 1 year for the San Andreas

T_a = 2π Aσ₀ / τ̇               critical period; response peaks here
      Ader et al. 2014, GJI 198 — ~20–200 yr for ordinary crust

responsive band:   t_n < T < T_a   ≈   1 year to 200 years
```

Which predicts:

| Constituent | Period | Predicted |
|---|---|---|
| Semidiurnal, diurnal (M2, S2, O1, K1) | 12–26 h | **damped** |
| Fortnightly, monthly (Mf, Mm) | 14–28 d | **damped** |
| Annual (Sa) | 365 d | **band edge** |
| Lunar nodal | 18.61 yr | **inside** |
| LOD / decadal | decades | **inside** |

**The short-period constituents most studies use are the ones physics predicts
should not work.** This independently explains weak semidiurnal correlation,
working annual/monsoon signals, the plausibility of Bendick & Bilham's decadal LOD
result, and why tremor responds at tidal periods when earthquakes do not.

⚠ `t_n ≥ 1 yr` is a **lab-to-field extrapolation**, the paper's own inference, not a
measurement. Carry it as a hypothesis with wide uncertainty. It is testable by
measuring the response spectrum — which is the point.

### Sensitivity: two separable parameters

Heimisson & Avouac 2020 (*GRL* 47) and Ader et al. 2014:

```
amplitude    R̃/r = Δτ / (aσ̄)         ← Aσ₀ ALONE sets sensitivity
peak period  T_a = 2π Aσ₀ / τ̇         ← τ̇ only sets WHERE the peak sits
```

Stressing rate does **not** set sensitivity. Because the two observables are
separable, the spectrum yields **Aσ₀ and τ̇ as independent measurements** — and τ̇,
accumulating stress, is what forecasting wants.

Confirmed in data: *EPSL* (2025) finds San Andreas LFE modulation amplitude
controlled by background effective stress, frequency-dependence by frictional
properties and nucleation time.

### Non-linearity

```
R(t) = R₀ · exp( S_T(t) / (Aσ₀) )        M = ⟨exp(S_T/Aσ₀)⟩ ≥ 1,  R₀ = r/M
```

- Ordinary crust: tidal stress ~10⁻³–10⁻⁴ MPa, Aσ₀ ~ 0.01–0.1 MPa → ratio 10⁻²–10⁻¹.
  **Linearisation fine.**
- Parkfield tremor: Aσ₀ = 6×10⁻⁴ MPa → ratio **0.2–2. Strongly non-linear**; carry
  the full exponential and M > 1.

### The invariant that constrains every claim

**⟨R⟩ = r exactly** (Heimisson & Avouac eq. 6). Oscillatory stress does **not change
the mean rate** — only the timing. Tides redistribute *when* events occur; they do
not create them.

Any code path, model, or output implying otherwise is a bug.

---

## 3. What PlanetaryHarmonics provides

Rust crate `ph-core`, tested, in the upstream repo. Consume it; do not reimplement.

| Module | Provides |
|---|---|
| `tidal` | Degree-2 tidal tensor `T_ij = (GM/d³)(3n̂ᵢn̂ⱼ − δᵢⱼ)`, linear superposition, symmetric eigendecomposition, principal axis (the "concentration node") |
| `harmonics` | Fourier angular encoding by angle-sum recurrence, least-squares amplitude/phase decomposition |
| `stats` | Generalised Schuster test, Schuster periodogram, geometric trial-period spacing, peak finding, time-shifted null distributions |
| `ephemeris` | Batched geometric states via RustSPICE, tidal body list, GM lookup |
| `catalog` | Minimal `Event`/`Catalog` types — deliberately thin, an interface not a substrate |

**Validated.** Phase 1 recovered five known deep moonquake periodicities from the
Apollo catalogue to better than 0.2% (13.609 d, 27.19 d, 27.567 d, 29.58 d,
206.19 d). The instrument works on data where tidal forcing is dominant.

**Access:** `pyo3`/`maturin` Python bindings are planned upstream. Until they land,
exchange via CSV/Parquet from a `ph-core` CLI. Do not fork the Rust.

**Every feature carries provenance** — reference frame, epoch system, kernel
versions, feature tier, units. A feature whose frame is ambiguous is worthless in a
statistical test, and across a repo boundary that ambiguity is easy to introduce.
Do not strip it.

---

## 4. What this repository must build

1. **USGS ComCat ingestion** + magnitude-of-completeness analysis per region and
   epoch
2. **GCMT focal mechanisms** → fault geometry (strike/dip/rake)
3. **Coulomb projection** — `ΔCFS = Δτ + μ′Δσₙ` and `dΔCFS/dt`.
   ⚠ *Check upstream first — this is physics and may land in `ph-core`.*
4. **ETAS baseline**, fitted and **frozen**
5. **Residual model** — `λ = λ_ETAS · exp(f_θ(features))`
6. **β(x,t) sensitivity field** — the primary research target
7. **Hydrological covariates** (GRACE/GLDAS) — see §6
8. **CSEP evaluation** via `pyCSEP`

---

## 5. Methodology — non-negotiable

**It is a point process, not classification.** Target the conditional intensity
`λ(x,y,t,m)`, trained on point-process log-likelihood. The integral term makes the
*absence* of earthquakes informative. Classification framing gives 99.9% accuracy
and learns nothing.

**Fit ETAS first, freeze it, learn only the residual.**
`λ = λ_ETAS · exp(f_θ(features))` with a prior pulling `f_θ → 0`. The model is then
structurally incapable of taking credit for clustering ETAS already explains, and
`exp(f_θ)` is directly interpretable: "features modulate rate by ±2.8%" is a real,
falsifiable claim.

**Time-shifted null, not random shuffling.** Shuffling destroys aftershock
autocorrelation and gives falsely tight nulls. Shift the whole catalogue by large
offsets (±1 to ±50 yr), preserving both the catalogue's clustering and the forcing's
structure, breaking only their alignment.

> ⚠ **Two caveats found the hard way, both in `ph-core::stats`:**
>
> 1. **The null is degenerate at a single exact frequency.** A global shift rotates
>    the phase cluster without changing its concentration, so `D²ₙ` is *exactly
>    invariant* and the null has **zero power**. Its power comes from the forcing
>    being **quasi-periodic** — real tidal phase is frequency-modulated by the nodal
>    and perigee cycles. **Compute phase from the full forcing, never a single
>    idealised constituent.**
> 2. **The empirical p-value floor is `1/(n+1)`.** p < 0.05 needs ≥20 shifts;
>    p < 0.005 needs ≥200. Budget shift count for the significance you intend to
>    claim.

**Mc filtering is mandatory.** Detection improved dramatically over the catalogue
span; uncorrected, that trend projects onto any long-period feature and
manufactures signal — which is fatal given §2's band prediction points at exactly
those periods.

**FDR on any feature scan.** Benjamini–Hochberg, or pre-register a small set.

**Report information gain per earthquake** in bits/event via the CSEP T-test. Run
N/S/M/L tests. Do not invent metrics.

---

## 6. Validity gates — run these BEFORE believing any positive result

These need no external data and are cheap. They are gates, not analyses.

**M2 vs S2 — cultural-noise discriminator.** S2 is *exactly* 12.000 h, locked to the
day–night cycle, and catalogue completeness varies with time of day (cultural noise
raises detection thresholds). S2 is also contaminated by the solar *thermal* tide,
which is not a body tide. M2 is 12.42 h, so it precesses through local solar time
over a lunar month and decorrelates from time-of-day artifacts.

- Signal at S2 but not M2 → **artifact**
- Signal at M2 ≥ S2, in the ratio the tidal potential predicts → real

Likewise **K1 (23.93 h) and S1 (24.00 h) are near-degenerate with the diurnal
detection artifact and are effectively unusable.** O1 (25.82 h) and P1 (24.07 h) are
safe.

**Explicit alias analysis.** Enumerate catalogue periodicities (daily, weekly,
seasonal maintenance, network upgrades), compute beats against the constituent
list, blacklist collisions.

**Two free amplitude knobs.** Lunar distance varies 5.5% over the anomalistic month
and tides go as 1/d³ → **18% amplitude modulation** at 27.55 d with known phase. The
18.61-yr nodal cycle modulates diurnal amplitudes ~±11%, giving **seven envelope
cycles** over a 130-year catalogue. Response tracking either at predicted depth and
phase is a fingerprint no instrumental effect plausibly mimics.

**Same-band constituent ratios.** Within a band (M2/S2/N2; O1/K1/P1) the transfer
function is effectively constant, so response ratios should equal the known
amplitude ratios. Isolates the amplitude law from the frequency law.

**Hydrological loading is a confounder *and* an instrument.** Seasonal groundwater,
snow and atmospheric loading are annual — and so are Sa and Ssa. Any annual
"celestial" correlation is confounded with hydrology by default. But hydrological
amplitude is independently measurable from GRACE/GRACE-FO and GLDAS, so model it
jointly: it becomes a **second probe of the transfer function** at annual period,
where tidal amplitude is weak. Precedent: *Science Advances* sciadv.ady6350.

---

## 7. Data — all free, no credentials

| Source | Contents |
|---|---|
| USGS ComCat | Global earthquake catalogue, public API |
| GCMT | Focal mechanisms |
| IRIS / EarthScope | Waveforms and catalogues |
| IERS EOP | Polar motion, length of day |
| GRACE / GRACE-FO, GLDAS | Hydrological loading |
| FES2014 / TPXO9 / GOT | Ocean tide models |
| NAIF | DE440/DE441 SPICE kernels |

**Literature access without institutional credentials:** author self-archived pages
first (highest yield), then **USGS Publications Warehouse** — USGS-authored work is
public domain and covers much of this field (Beeler, Cochran, Hardebeck) — then
arXiv, **ESS Open Archive** (AGU preprints: GRL, JGR), EarthArXiv, and the Unpaywall
API by DOI.

---

## 8. Prior art — use it, do not rebuild

- **RECAST** — neural temporal point process, GRU encoder-decoder; matches or beats
  ETAS on Southern California given enough events. `github.com/keliankaz/recast`
- **EarthquakeNPP** — benchmark suite, arXiv:2410.08226
- **pyCSEP** — the evaluation framework, `cseptesting.org`
- **Pulsar timing array methods** — red noise as a Fourier-domain Gaussian process,
  analytic marginalisation over Fourier amplitudes, empirical false-alarm
  estimation. The same statistical problem shape (small periodic signal, strongly
  red noise, rigorous false-alarm control) and far more rigorous than the Schuster
  test, which assumes independence. Tooling: `enterprise`.

---

## 9. Terminology

| Avoid | Use |
|---|---|
| prediction | probabilistic forecast, conditional rate estimate |
| nexus event | coherence maximum, commensurability alignment |
| resonance point | phase coherence peak |
| astrology feature | derived planetary feature |
| base-N harmonics | harmonic order N of relative longitude |

Label every feature by tier: **A** established physics, **B** established method /
novel application, **C** exploratory. Report Tier C separately and never as
mechanism.

---

## 10. Mistakes already made — do not repeat

1. **Event count off by 100×** — recorded 10⁵–10⁶ from a search summary; the paper
   says ~10⁴. *Verify against primary sources, not abstracts.*
2. **Wrong response model** — assumed a simple 1/ω low-pass giving a monotonic
   constituent ranking. The literature gives a two-regime, band-limited response.
   *An abstract-based reconstruction was directionally right and structurally
   wrong.*
3. **Conflated two parameters** — claimed sensitivity ∝ τ̇/(Aσ₀). Amplitude depends
   on Aσ₀ alone; τ̇ only sets the peak period. The corrected form is more useful.
4. **Claimed novelty that existed since 1995** — the "generalised Doodson expansion
   with planetary arguments" is HW95: 12,935 waves, 1,483 of them planetary.
   *Search for prior art before claiming a contribution.*
5. **Reasoned from a null with no power** — see §5's degeneracy caveat. Found only
   by writing a test that failed.

---

## 11. Open questions

- Does the 1 yr – 200 yr band prediction survive measurement at global scale? This
  is the project's central question and is **not** answerable from the literature.
- Free access to HW95 / ETERNA catalogue data files. The Kudryavtsev & Cionco
  (2025) release (arXiv:2508.18111) ships in HW95 format and is the best free route.
- Exact transfer function near `T ≈ t_a` — Heimisson & Avouac state their model is
  **invalid** there, which is precisely where the spectral peak sits. Candidate for
  symbolic regression.
- What gets published if the answer is null? **Decide now:** a rigorous upper bound
  on celestial–terrestrial coupling, plus the moonquake and tremor validations.
  Pre-committing removes the incentive to keep scanning until something crosses
  p < 0.05.
