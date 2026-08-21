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

---

## 2026-08-21 — Second pass: rate-and-state response theory

**Scope:** agenda items 1 (transfer function), 2 (spectral prediction test),
PTA statistical methods. Triggered by hypothesis 9.

### Finding 1 — the spectral prediction in hypothesis 9 was wrong ✗

The first-pass guess was a simple low-pass with response ∝ 1/ω across the tidal
band, giving a monotonic ranking of constituents by period (nodal > Sa > Ssa > Mm >
Mf > semidiurnal).

**The literature does not support this.** The actual result (Ader, Lapusta, Avouac
& Ampuero 2014, *GJI* 198, 385–413; Heimisson & Avouac 2020, *GRL*) is a
**two-regime response** separated by a critical period.

Define the Dieterich characteristic time and critical period:

```
t_a = a σ / τ̇          T_c ≈ 2π t_a
```

- **T ≪ t_a** — seismicity rate **tracks stress**. Response amplitude ∝ A,
  approximately flat in period.
- **T ≫ t_a** — seismicity rate **tracks stressing rate**. Response amplitude
  ∝ A/T, *falling* with increasing period.

So long-period constituents are **not** automatically favoured. The nodal-term
ranking prediction is dead.

**What survives, and is sharper:** the response has a **knee, and plausibly a peak,
near T ≈ t_a**. A 2012 AGU abstract from the same group describes an amplitude
maximum at a characteristic period *decreasing at both smaller and larger periods*
— i.e. **band-pass**.

*Reconstruction (unverified — see caveat):* band-pass behaviour requires **two**
timescales, not one. The Dieterich relaxation time t_a governs the long-period
roll-off; a separate **nucleation duration** governs the short-period cutoff
(Beeler & Lockner's original argument). Response peaks between them.

If correct, this is *better* for the inverse problem than the original guess: the
spectrum would identify **two** fault properties rather than one, and a peak is far
more diagnostic than a slope — an artifact can produce a slope, but is unlikely to
produce a peak at a physically predicted location.

**Caveat, important:** Heimisson & Avouac state their analytical model is
**not valid for periods similar to t_a** — precisely the region of interest. The
band-pass reconstruction above is inference from abstracts, not from the papers
themselves. Both primary sources are paywalled and were not obtained this pass.
**Do not build on this until the actual equations are in hand.**

### Finding 2 — hypothesis 1 is already an active research programme

Frequency-dependent tidal sensitivity as a probe of fault properties is
**published work**, not a novel idea:

- **Beeler et al. (2018)**, *JGR* 123 — *Constraints on friction, dilatancy,
  diffusivity, and effective stress from low-frequency earthquake rates on the deep
  San Andreas Fault*
- **EPSL (2025)**, *Probing lower-crustal fault properties with frequency-dependent
  tidal tremor triggering* — San Andreas LFEs. Finds that **modulation amplitude is
  controlled mainly by background effective stress**, while the **diurnal vs.
  semidiurnal variation is controlled by frictional properties and nucleation
  time.**

**This substantially validates hypothesis 9** — β does track effective normal
stress, as published work confirms. It also means the core idea is not ours to
claim.

*Strategic consequence.* Novelty must move to where the existing work is not:

1. **Long-period constituents.** Published work is diurnal/semidiurnal. Mf, Mm,
   Ssa, Sa, and the 18.61 yr nodal term are largely unexploited.
2. **Global scale.** Existing studies are single-fault (San Andreas, Cascadia). A
   global β field is genuinely new.
3. **The multi-basis / generalised Doodson extension.** Still unclaimed.
4. **Other bodies.** Heimisson & Avouac explicitly frame their model as applicable
   to "other solid-surface bodies" — supporting the moonquake plan.

This is a better position than it sounds: the mechanism is now *established
literature we can cite* rather than something we must prove.

### Finding 3 — pulsar timing array methods transfer cleanly ✓

Hypothesis 10 confirmed as viable. PTA analysis solves the same statistical
problem shape — small periodic signal, strongly red noise, rigorous false-alarm
control:

- Red noise modelled as a **Gaussian process in the Fourier domain** with
  harmonically related sinusoids
- **Analytic marginalisation** over individual Fourier component amplitudes while
  fitting the power-law amplitude and spectral index of the underlying process
- **Hierarchical Bayesian** treatment with hyperparameter marginalisation
- False-alarm null estimated **empirically**, because correlations make analytic
  nulls unreliable

Mature tooling exists (`enterprise`, and the EPTA/NANOGrav/MeerKAT analysis
stacks).

**Direct mapping to our problem:** catalogue red noise (aftershock clustering, Mc
drift, hydrological loading) becomes a fitted power-law Gaussian process; tidal
harmonics are the deterministic signal on top. Strictly more rigorous than the
Schuster test, which assumes independent samples and has no red-noise treatment.

Pleasingly, PTA's empirical false-alarm estimation is the same instinct as our
time-shifted null (doc 04 §6a) — independent convergence on the same safeguard.

### Finding 4 — the hydrological confounder is also the control

*Tidal and hydrological seismicity modulations reveal pore fluid diffusion during
earthquake nucleation* (**Science Advances**, sciadv.ady6350) treats tidal and
seasonal/hydrological modulation **jointly** to infer pore fluid diffusion.

Reframes doc 07's "most dangerous confounder": hydrological loading is a
**second periodic forcing with independently known amplitude and phase**. Rather
than merely controlling for it, use it as a **calibration signal** — it probes the
same transfer function at annual period, where tidal amplitude is weak.

Two forcings at different periods constrain the transfer function far better than
one. The confounder becomes an instrument.

### Not obtained this pass

- Ader et al. 2014 full text (cert failure on the Caltech mirror; paywalled at OUP)
- Heimisson & Avouac 2020 full text
- Science Advances paper full text (403)
- CTE / HW95 constituent catalogues

All four remain on the critical path. Acquiring the actual equations is the
blocking task for the next pass.

---

## 2026-08-21 — Third pass: primary sources obtained

**Access note.** No institutional access. Both critical-path papers were obtained
free from the author's public page at `web.gps.caltech.edu/~avouac/publications/`
(self-archived copies). **Author pages are the highest-yield OSINT route** — see
[11-osint-access.md](11-osint-access.md).

### The actual equations

**Heimisson & Avouac (2020), GRL 47** — general oscillatory stress:

```
Short period (T ≪ t_a):     R(t) = R₀ · exp( S_T(t) / (Aσ₀) )        (eq. 4, 8)
                            R₀ = r / M,   M = ⟨exp(S_T/Aσ₀)⟩ ≥ 1     (eq. 5, 7)

Long period  (T ≫ t_a):     R(t)/r ≈ 1 / (1 − t_a Ṡ_T(t)/(Aσ₀))     (eq. 9)

Characteristic time:        t_a = Aσ₀ / ṡ₀
```

with `S(t) = τ(t) − μσ(t)` the *modified* Coulomb stress, and μ = τ₀/σ₀ − α where
α is the Linker–Dieterich constant (0–0.25). **μ here is not a friction coefficient
in the traditional sense** — worth noting, since it is easy to misread.

**Ader et al. (2014), GJI 198** — finite rate-and-state fault:

```
Critical period:            T_a = 2π t_a = 2π aσ̄ / τ̇_a               (eq. 8)
Small-perturbation amplitude: R̃/r = Δτ / (aσ̄)                        (eq. B7)
```

### Finding 1 — band-pass CONFIRMED, and it comes from fault finiteness ✓

The second-pass reconstruction was right. Ader et al. state explicitly there is
"a critical period T_a, **at which the amplitude of the seismicity response
peaks**."

Structure of the finite-fault response:

| Regime | Amplitude | Phase |
|---|---|---|
| T < T_a | **increases with period** | ≈ 0, drifting to about **−π/4** as T → T_a |
| T ≥ T_a | falls; Coulomb-failure behaviour, rate ∝ stress rate | — |

The rising limb below T_a is a **finite-fault effect not present in the
spring–slider (SRM) model**, and so is the gradual −π/4 phase drift. Both are
therefore **discriminators between fault models**, exactly as hypothesis 2 wanted.

Also: finite-fault response amplitude is **always much larger than SRM predicts,
sometimes by more than an order of magnitude**. Consequence stated by the authors:
**aσ̄ values inferred from observations using SRM are systematically
underestimated.** Heimisson & Avouac note a second underestimation bias — Aσ₀ is
underestimated if stress is approximated by a single harmonic. Two independent
biases, both in the same direction.

### Finding 2 — hypothesis 9 was structurally wrong; the corrected version is better ✓

The synthesis claimed β ∝ τ̇/(Aσ₀). **That conflates two separate roles.**

```
Response AMPLITUDE      ←  Aσ₀ alone        R̃/r = Δτ/(aσ̄)
Peak/critical PERIOD    ←  T_a = 2π Aσ₀/τ̇   (sets which periods are favoured)
```

Stressing rate τ̇ does **not** set sensitivity. It only sets where the spectral peak
sits.

**This is better than the original claim**, because the two observables are
separable:

> Measure response **amplitude** → recover **Aσ₀** (effective normal stress × a).
> Measure the **peak period** T_a → recover **Aσ₀/τ̇**, hence **τ̇** independently.

**The harmonic spectrum yields effective normal stress and stressing rate as two
independent measurements.** And τ̇ — accumulating stress — is precisely the quantity
earthquake forecasting wants.

This also matches EPSL (2025) exactly: modulation *amplitude* controlled by
background effective stress; *frequency dependence* controlled by frictional
properties and nucleation time. Two observables, two parameters, confirmed in
published data.

### Finding 3 — numbers, and a resurrection of the long-period ranking

`T_a = 2π Aσ₀/τ̇`. Plugging in literature values:

| Setting | Aσ₀ | τ̇ | T_a |
|---|---|---|---|
| Ordinary crust (aftershock studies) | 0.01–0.1 MPa | ~0.003 MPa/yr | **~20–200 yr** |
| Parkfield tremor (Thomas et al. 2012) | 6 × 10⁻⁴ MPa | ~0.003 MPa/yr | ~1.3 yr |
| Tremor during slow slip (τ̇ ~100× higher) | 6 × 10⁻⁴ MPa | elevated | **~days** |
| Ader et al. simulation (σ̄ = 5 MPa) | — | — | ≈ 0.03 yr (~11 d) |

Two consequences:

**1. Tremor's T_a sits inside the tidal band.** That explains strong, clearly
frequency-dependent tidal modulation of tremor — and the fortnightly modulation —
without needing any special pleading. The mechanism for doc 05's "tremor is the
high-SNR testbed" is now quantitative.

**2. For ordinary crust, T_a ~ decades — so every tidal constituent, out to the
18.61 yr nodal term, lies in the T < T_a regime where amplitude RISES with
period.** The long-period ranking (nodal > Sa > Ssa > Mm > Mf > diurnal >
semidiurnal) is **plausible again** — but as the *rising limb toward a peak at
T_a*, not as the 1/ω low-pass falloff of the first pass. Same ordering, completely
different mechanism, and now with a theoretical amplitude curve to fit rather than
a bare slope.

### Finding 4 — two incidental confirmations

- **Ader et al. use Schuster spectra** (their Fig. 9) as the analysis tool. Our
  generalised-Schuster approach (doc 02) is the field's own method.
- **Harmonic forcing measurably modifies the magnitude distribution** in their
  simulations — whole magnitude bands "completely disappeared" under perturbation.
  Independent theoretical support for the Ide et al. (2016) b-value target
  (doc 04 §4).

### Non-linearity caveat

Response is `exp(S_T/Aσ₀)`, not linear. Linearisation holds when |S_T|/Aσ₀ ≪ 1:

- Ordinary crust: tidal stress ~10⁻³–10⁻⁴ MPa, Aσ₀ ~ 0.01–0.1 MPa → ratio 10⁻²–10⁻¹.
  **Linear is fine.**
- Parkfield tremor: Aσ₀ = 6×10⁻⁴ MPa → ratio **0.2–2. Strongly non-linear.**

So the tremor testbed needs the full exponential form, and M > 1 must be carried.
Note also **⟨R⟩ = r exactly** (eq. 6): oscillatory stress does not change the mean
rate, only the timing. Tides redistribute *when*, they do not create events. This
belongs in doc 00's claim discipline.

### RustSPICE update

Pulled to `c7f180b` (tag `v0.1.0`). Substantially restructured — see revised
[10-rustspice-requirements.md](10-rustspice-requirements.md).
