# 07 — Research Log

Dated record of literature passes, findings, and decisions. Append-only.

---

## 2026-08-21 — First literature pass

**Scope:** tidal triggering of earthquakes; angular feature encoding; ML
forecasting architecture; validation datasets.

### Headline findings

**1. The effect is real, small, and setting-dependent.**
Range spans null results (Southern California, general catalogues) to a
**factor-of-3 rate change** (Cochran et al. 2004, shallow thrust faults under high
ocean-tidal loading, p < 10⁻⁴). Métivier et al. 2009 found ~99% confidence on
442,412 NEIC events, with the anomaly *larger for smaller and shallower* events.

Neither "tides cause earthquakes" nor "tides have no effect" is supportable. The
defensible position is that the effect exists, is a few percent globally, and
grows large where tidal stress amplitude is large.

**2. Sample-size requirement is a hard design constraint.**
Beeler & Lockner (2003): **10⁵–10⁶ earthquakes** needed for a statistically robust
correlation. A 10× stress amplitude increase cuts the requirement 100×.

*Decision:* full USGS catalogue pull is a hard requirement. Power analysis before
any modelling. Region selection driven by tidal stress amplitude, which is a
physical prediction rather than cherry-picking.

**3. The nucleation-time argument predicts long periods should win.**
Beeler & Lockner: correlation is facilitated when forcing period **exceeds the
characteristic nucleation time** of frictional instability. Tidal semidiurnal
periods sit near or below it — hence weak correlation.

This predicts that **long-period constituents (Mf 13.66 d, Mm 27.55 d, Ssa 182.6 d,
Sa 365.26 d, nodal 18.61 yr) should correlate better than semidiurnal ones.**
Sharp, testable, and it is exactly the regime where orbital/planetary cycles live.

*This is the most useful physical hook found in the entire pass.* See
[08-hypotheses.md](08-hypotheses.md) §1 for what it becomes.

**4. Phase is diagnostic, and can be counter-intuitive.**
Métivier: events cluster at the **uplift** phase. Scholz et al. 2019: mid-ocean
ridge events peak at **low tide** — opposite the naive expectation, explained via
magma chamber inflation.

*Decision:* angular encodings must retain the sine component so phase is learned,
never assumed. Settles the encoding question in [02-angular-encoding.md](02-angular-encoding.md).

**5. Tidal sensitivity is a stress-state proxy — the strongest reframe available.**
Beaucé et al. 2023: tidal sensitivity of seismicity rose **~1.5 years before** the
M7.1 Ridgecrest mainshock (>150,000 event catalogue). Ide et al. 2016: very large
events cluster near maximum tidal amplitude, modulating the **size distribution**
rather than the rate.

*Decision:* primary ML target changed from rate modulation to the **time-varying
sensitivity field β(x,t)**. See [04-ml-architecture.md](04-ml-architecture.md) §4.

**6. Tremor is an order of magnitude easier than earthquakes.**
Rubinstein et al. 2008: clear tremor pulsing at 12.4 h and 24–25 h; tides
"influence the genesis of tremor much more effectively than they do the genesis of
normal earthquakes." Ide et al.: tremor rate rises **exponentially** with tidal
stress.

**7. Deep moonquakes are tidally dominated, not tidally nudged.**
Apollo PSE catalogue. Known monthly / 206-day / 6-year periodicities.

*Decision:* adopt a graded validation ladder — **moonquakes → tremor →
earthquakes**. Each phase has a known-answer test before the next. See
[05-research-frontier.md](05-research-frontier.md).

### Methodological finding

**Holt & Newman (2025), arXiv:2508.07064**, *Tidal Triggering of M7+ Earthquakes by
Jupiter*: **28,782 chi-squared tests**, **1,071 significant** reported. At α = 0.05
chance predicts ≈ 1,439. They found *fewer significant results than the null
expects*.

(Tests are heavily correlated — overlapping 5–45 day intervals — so the effective
independent count is well below 28,782 and 1,439 is not the exact expectation. The
conclusion is unchanged: uncorrected and unadjusted for correlation, the result is
indistinguishable from noise.)

*Decision:* FDR correction mandatory on any feature scan; time-shifted null
distribution rather than random shuffling. Documented as the project's canonical
cautionary example in [00-framing.md](00-framing.md).

### Physics check performed

Tidal effect scales as M/d³. Computed:

| Body | M/d³ (SI) | vs. Moon |
|---|---|---|
| Moon | 1.30 × 10⁻³ | 1 |
| Sun | 5.9 × 10⁻⁴ | 0.46 |
| Venus (closest) | 8.9 × 10⁻⁸ | 7 × 10⁻⁵ |
| Jupiter (closest) | 9.4 × 10⁻⁹ | 7 × 10⁻⁶ |

Matches published figures (Venus ~10⁴× weaker than Sun+Moon combined).

*Conclusion:* direct planetary tidal triggering is below detection and must be
documented as such. Planetary influence enters legitimately via **orbital
perturbation of the lunisolar tide** and possibly via **principal-axis deflection**
(see [03-tidal-tensor.md](03-tidal-tensor.md) §5).

*Also recorded:* a uniform gravitational field produces no stress. Barycentres,
Lagrange points, and inter-planetary "gravitational nulls" are not stress
mechanisms. Written into doc 03 explicitly to prevent the reasoning error.

### Tooling identified

| Tool | Purpose |
|---|---|
| `pyCSEP` (cseptesting.org) | N/S/M/L tests, information gain per earthquake |
| RECAST (`github.com/keliankaz/recast`) | Neural temporal point process, GRU encoder-decoder |
| EarthquakeNPP (arXiv:2410.08226) | Benchmark suite for neural point processes |
| SPOTL (Agnew) | Ocean tidal loading |
| GCMT | Focal mechanisms → fault geometry |
| IERS | Polar motion, LOD series |
| Apollo PSE (Nakamura) | Deep moonquake catalogue |

### Confounder identified — high priority

**Seasonal hydrological loading is annual. So are Sa and Ssa.** Any annual-period
"celestial" correlation is confounded with groundwater, snow load, and atmospheric
pressure by default. Must enter as an explicit covariate, not be ignored.

Currently the most dangerous unaddressed confounder in the project.

### Open threads for next pass

1. Rate-and-state transfer function — analytic form of amplitude and phase
   response to periodic forcing (Ader & Avouac; Dieterich 1994; Beeler & Lockner)
2. Cartwright–Tayler–Edden tidal potential catalogue — amplitudes per constituent
3. Farrell (1972) load Love numbers and Green's functions for ocean loading
4. Schuster test null distribution generalised to harmonic order n
5. ETAS space–time–magnitude likelihood, exact form
6. Deep moonquake catalogue structure and prior periodicity analyses

Carried into [09-deep-dive-agenda.md](09-deep-dive-agenda.md).
