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

---

## 2026-08-21 — Fourth pass: USGS public-domain set and tidal catalogues

**Access.** USGS Publications Warehouse and the USGS Rock Physics Lab page
(`earthquake.usgs.gov/research/eqproc/rockphysics/pubs.php`) delivered the entire
Beeler set free — USGS-authored work is public domain. Confirms doc 11's
prediction that this was the highest-value untried route.

Obtained: Beeler & Lockner (2003), Beeler et al. (2013) LFE rheology, Beeler et al.
(2016) effective stress, Lockner & Beeler (1999) premonitory slip.

### Finding 1 — the event-count requirement was wrong by two orders of magnitude ✓✓

**Previously recorded: 10⁵–10⁶ events. Actual: ~10⁴.**

That figure came from a search summary in the first pass and was never verified
against the paper. The paper itself, equation 18:

```
N ≥ ln(P_rw) / ( Δτ_u / (2 a σ_n) )²
```

**N scales as the inverse square of normalised stress amplitude.** Worked examples
from the paper:

| Stress amplitude | Events needed |
|---|---|
| Δτ = 0.01 MPa (a = 0.0045, σ_n = 18 MPa/km, 5–15 km) | **6.2 × 10³ – 5.5 × 10⁴** |
| Δτ = 0.1 MPa (10× larger) | **60 – 541** |

And from the abstract: daily Earth tides require **">13,000 earthquakes to
detect."**

**Consequence — the project is far more feasible than doc 01 recorded.** USGS
ComCat holds millions of events; 13,000 is trivially available, and many *regional*
subsets clear the bar comfortably. We are not power-limited at the global scale,
and the doc 04 §7 power analysis should be re-run with the correct law. The
"hard requirement" framing for the full catalogue pull stands, but for coverage and
confounder control rather than raw statistical power.

**Doc 01 corrected.**

### Finding 2 — the second timescale, and a startling band prediction ★★

Beeler & Lockner identify the **nucleation duration t_n**, distinct from Dieterich's
t_a. Their result:

> Fault strength's second-order dependence on sliding rate "determines the duration
> of nucleation and **damps the response to stress change at frequencies greater
> than 1/t_n**."

And, critically:

> "The experiments suggest that the **minimum typical duration of earthquake
> nucleation on the San Andreas fault system is 1 year.**"

This is the missing short-period cutoff — the second timescale the band-pass
reconstruction required (third pass). Combining with Ader's critical period:

```
responsive band:    t_n  <  T  <  T_a
                   ~1 yr  <  T  <  ~20–200 yr        (ordinary crust)
```

**Everything faster than about a year is damped by nucleation. Everything slower
than T_a rolls off.**

If this holds, the implications are severe and specific:

| Constituent | Period | Predicted |
|---|---|---|
| Semidiurnal (M2, S2, N2) | 12–12.4 h | **damped** |
| Diurnal (O1, K1, P1) | 24–26 h | **damped** |
| Fortnightly, monthly (Mf, Mm) | 13.7–27.6 d | **damped** |
| Semiannual (Ssa) | 182.6 d | marginal |
| **Annual (Sa)** | 365.26 d | **at the band edge** |
| **Lunar nodal** | 18.61 yr | **inside the band** |
| **LOD / decadal** | decades | **inside the band** |

**This is the strongest theoretical result the project has found.** It says the
responsive band for ordinary earthquakes is *years to decades* — orbital and
long-period timescales — and that the short-period tidal constituents everyone
studies are the ones physics predicts should NOT work.

It independently explains: why semidiurnal tidal correlation is weak everywhere;
why annual/monsoon signals correlate (Ader's Himalayan motivation); why Bendick &
Bilham's decadal LOD correlation is not absurd; and why tremor — with t_n and T_a
both far shorter — responds at tidal periods while earthquakes do not.

**It also justifies the project's long-period focus on physical grounds.** Not
because long cycles are meaningful, but because the fault itself low-passes
everything faster than a year.

⚠ **Caveat.** t_n ≥ 1 yr is an *extrapolation from laboratory experiments to the
San Andreas*, not a field measurement. It is the paper's own inference and should
be carried as a hypothesis with a large uncertainty, not a constant. But it is
published, mechanistic, and testable — precisely by measuring the response spectrum
(hypothesis 12).

### Finding 3 — HW95 already contains the planetary expansion ✗

**Doc 02's central novelty claim was wrong.**

The **HW95 catalogue** (Hartmann & Wenzel 1995, *GRL* 22) is a harmonic development
of the tide-generating potential containing **12,935 waves, of which 1,483 are due
to direct planetary effects** — Venus, Jupiter, Mars, Mercury, Saturn. Based on
DE200, spanning 1850–2150. Scope: Moon to degree 6, Sun to degree 3, planets to
degree 2.

So the "generalised Doodson expansion including planetary arguments" that doc 02
proposed as our novel contribution **was published in 1995.**

**This is a gift, not a setback.** We do not have to derive it — 1,483
rigorously-computed planetary tidal waves with correct amplitudes, produced by
geodesists, are exactly the multi-body harmonic feature basis the project wants.
Using an authoritative published basis is a far stronger position than deriving our
own.

**Corrected novelty claim.** Not the expansion. Rather:
1. Using HW95/KSM03 planetary waves as an **ML feature basis tested against
   seismicity** — not done
2. **Broadband response spectroscopy** across six decades (hypothesis 12)
3. A **global β field** — published work is single-fault
4. The **moonquake → tremor → earthquake validation ladder**

**Doc 02 corrected.**

### Finding 4 — the modern successor catalogue, and it is in our band

**Kudryavtsev & Cionco (2025)**, *New and updated long-periodic terms in harmonic
development of the Earth tide-generating potential*, arXiv:2508.18111, accepted in
*Geodesy and Geodynamics*. **Free on arXiv.**

- **38 terms with period longer than ~18 years**, amplitude ≥ 10⁻⁸ m²s⁻²
- A relatively large term at **~7.4 kyr**
- **Several NEW waves near the 18.61-year lunar nodal period**
- Built on **DE441**, planetary terms included, 30,000+ year spectral analysis
- Updates the **KSM03** expansion (Kudryavtsev 2004)
- Released in both **HW95 and KSM03 format**

This lands exactly in the band Finding 2 predicts is responsive. New waves near
18.61 yr are directly relevant. Highest-priority acquisition.

### Finding 5 — LFE effective stress is kilopascal-scale

Beeler et al. (2013), *Inferring fault rheology from low-frequency earthquakes on
the San Andreas*: observed tidal modulation "restricts ambient stress to be at most
a **few kilopascal**."

A few kPa = 10⁻³ MPa, consistent with Thomas et al.'s Aσ₀ = 6×10⁻⁴ MPa for
Parkfield tremor. Confirms the low-effective-stress mechanism for LFE/tremor tidal
sensitivity and supports tremor as the high-SNR testbed.

### Catalogue access status

**ETERNA 3.30** bundles seven tidal potential catalogues — Doodson (1921),
Cartwright et al. (1971, 1973), Büllesfeld (1985), Tamura (1987), Xi (1989),
Roosbeek (1996), Hartmann & Wenzel (1995). HW95 gives better than 1 nGal accuracy.

⚠ ETERNA 3.30 was distributed on CD-ROM for US$300 in 1997. **Modern free access
route not yet established.** Options to try next pass: ETERNA 3.4 / ETERNA-x
redistributions, IGETS / ICET (International Center for Earth Tides), PyETERNA and
other open reimplementations, and the KSM03/HW95-format catalogue released with
Kudryavtsev & Cionco (2025) — which may be the cleanest free path.

**RATGP95** (Roosbeek 1996) has a free PDF at Oxford Academic.

### Open

- ETERNA / HW95 catalogue data files — free source not yet found
- Kudryavtsev & Cionco (2025) catalogue data
- Ader PhD thesis (Caltech THESIS, open) — likely fuller derivations

---

## 2026-08-21 — Phase 1 validation: PASSED ✓

First end-to-end run. `cargo run --release --example moonquake_validation`.

**Data.** Apollo PSE expanded event catalogue, `levent.1008weber.csv`, from the PDS
Geosciences Node — public domain, no credentials. Fetch script at
`scripts/fetch-apollo.sh`.

**Ingestion finding — use `T2`, not `T1`.** The catalogue carries two
classification columns: `T1` original, `T2` revised (Nakamura 2005). They disagree
for 1,471 events, most of them originally logged as meteoroid impacts (`M`) and
later reclassified as deep moonquakes (`A`). Counting `T1` finds **1,359** deep
moonquakes; `T2` finds **7,082** — the figure the literature reports. Parsing to
usable times yields **6,954** (128 rows have unparseable time fields and are
skipped). Nest A1 has **424** events.

**Result — five of six known periodicities recovered, all nests (N = 6,954):**

| Period | Literature | Recovered | Error | Power |
|---|---|---|---|---|
| Half-month | 13.60 d | **13.609 d** | +0.07% | 869.4 |
| Draconic month | 27.212 d | **27.19 d** | −0.09% | 81.2 |
| Anomalistic month | 27.5546 d | **27.567 d** | +0.05% | 235.5 |
| Synodic month | 29.5306 d | **29.58 d** | +0.16% | 83.5 |
| Solar perturbation | ~206 d | **206.19 d** | +0.09% | 111.7 |
| Draconic/anomalistic beat | ~2190 d | no peak | — | — |

**Nest A1 (N = 424):** the draconic month is the **strongest** peak (power 131.5),
matching the literature's characterisation of A1. The 6-year beat appears at
2083 d (−4.87%, power 43.8) — acceptable given the record spans only ~1.4 cycles.

**Independent check of the "known" periods.** Both long periods derive from lunar
month beats, computed here from first principles:

```
draconic/anomalistic beat      = 2188.6 d    (literature ~2190, the 6-year cycle)
(anomalistic/synodic beat) / 2 =  205.9 d    (literature ~206)
```

So the 206-day period is **half the perigee–syzygy ("full moon") cycle**, and the
6-year period is the draconic/anomalistic beat. Neither is folklore.

### Unexplained peak — flag for the null

All-nests shows the **second strongest peak of the whole spectrum at 1886.6 d
(power 625.0)**, which matches no known deep moonquake period.

Most likely a **catalogue artifact**: the record spans ~2,920 days, Apollo stations
came online between 1969 and 1972 and ran to 1977, so overall event rate varies on
multi-year timescales and the observational envelope has its own spectral
structure. A 901 d peak (power 128.6) is similarly suspect.

**This is exactly what the time-shifted null exists to catch**, and it is a useful
early demonstration that raw periodogram power is not evidence. Next step for
Phase 1 is to run the null against these peaks and confirm the known periods
survive while 1886 d does not.

Smaller peaks at 6.805 d (= 13.609/2) and, in A1, 9.11 d (≈ 27.19/3) are harmonics
of the recovered periods, as expected.

### Assessment

**The instrument works.** The pipeline recovers five known tidal periodicities from
a real catalogue to better than 0.2%, on data where tidal forcing is the dominant
mechanism. Per `docs/12-build-plan.md` this was the gate on everything downstream.

Not yet done in Phase 1: the time-shift null on these peaks, and the ephemeris path
— computing actual tidal phase from Earth–Moon–Sun geometry rather than folding on
trial periods. The latter needs lunar kernels and is what carries forward to
Phases 2 and 3.

---

## 2026-08-21 — A3/A4: tidal phase implemented, and the null kills a 10⁻⁸⁹ result

**Code:** `ph-core::phase`, `examples/moonquake_tidal_phase.rs`.

### What was built

Tidal phase by the standard construction (Tanaka et al. 2002 and successors):
sample the forcing finely, locate successive **maxima**, assign each event a phase
by interpolating between the maxima bracketing it. Peak times refined by parabolic
interpolation for sub-sample accuracy.

Forcing scalar: largest eigenvalue of the lunar tidal tensor from real Earth + Sun
geometry in `MOON_PA`. Because the Moon is tidally locked, the tensor's
*orientation* barely moves and its *magnitude* varies as GM/d³ — so maxima are
lunar perigees and phase 0 means perigee.

**Self-check passed:** mean interval between forcing maxima came out **27.539 d**
against the anomalistic month's 27.555 d, confirming the scalar tracks lunar
distance as intended.

### The result

6,954 deep moonquakes, all phased, none dropped.

```text
analytic Schuster, order 1:   D²/N = 204.9   p = 1.07e-89   phase −92.1°
                     order 2: D²/N =  39.6   p = 6.52e-18   phase −91.7°

time-shift null, 308 offsets: observed D² = 1.42e6
                              null max    = 1.13e7      empirical p = 0.699
```

**The analytic test says p ≈ 10⁻⁸⁹. The time-shift null says p = 0.70.**

Nearly a third of shifted realisations cluster at least as strongly as the true
alignment, and the strongest cluster eight times harder. The analytic p-value
overstates significance by roughly **88 orders of magnitude**.

### Why — the degeneracy is broader than documented

The second-pass caveat said the time-shift null is degenerate against a *single
exact frequency*. That was too narrow. **Quasi-periodicity in the forcing is not
sufficient.** What matters is whether the **catalogue** shares the forcing's period.

Deep moonquakes are themselves locked near the anomalistic month. So is the
forcing. A global shift then rotates the phase cluster to a different phase while
leaving its concentration intact — exactly the invariance, arriving through the
catalogue rather than through the forcing.

Consequence: for this dataset, **"do events cluster in tidal phase?" is not a
falsifiable question.** Both series are periodic at the same rate, so clustering is
guaranteed at *some* phase.

### What this does and does not mean

**It does not mean deep moonquakes are not tidally driven.** The literature is
unambiguous that they are, and our own periodogram recovers five tidal periods to
better than 0.21%.

**It means this particular test cannot demonstrate it** — and, more importantly,
that reporting the analytic p-value alone would have been a false positive of
spectacular size. The null did its job. This is the first time in the project that
the machinery has caught something the naive analysis would have got badly wrong,
and it caught it on the dataset we chose *because* we knew the answer.

**The falsifiable question is which phase**, and whether it is consistent across
independent sub-populations in a way a failure mechanism predicts. Which is
precisely why Weber, Bills & Johnson (2009) work **per moonquake nest** rather than
on the pooled catalogue — a design choice whose necessity is now measured rather
than inherited.

### Consequences for the plan

- **A4 does not pass as specified.** The pooled-catalogue null cannot resolve the
  1886 d and 901 d artifacts either, for the same reason. A4 needs restating as a
  per-nest test.
- **A5/A6 rise in priority.** Per-nest Coulomb projection against Weber's cluster
  constraints is now the only route to a falsifiable Phase 1 claim, not merely a
  second validation.
- **Carry this to Phase 3.** Earthquake catalogues have strong internal
  periodicity too — aftershock sequences, seasonal detection cycles. The same trap
  is waiting, and the analytic Schuster p-value will be just as wrong there.
- `nakamura_2005_dm_locations.csv` (106 nests, already downloaded) becomes the
  next thing to use.

---

## 2026-08-22 — A5/A6: per-nest Coulomb search, and a second null artifact

**Code:** `ph-core::fault`, `examples/moonquake_nest_coulomb.rs`. 43 tests.

### What was built

Fault-plane resolution in the Aki & Richards convention: normal and slip vectors
from strike/dip/rake, traction decomposition, `ΔCFS = τ + μ′σₙ`, and rotation of a
body-fixed tensor into a local North-East-Down frame. Verified against invariants —
unit and orthogonal basis vectors, `|t|² = σₙ² + |τ|²`, rotation-invariant trace and
eigenvalues.

**Two things worth recording.**

*Depth does not matter.* The degree-2 tide-generating potential goes as
`r²·P₂(cos ψ)`, so its second derivatives are **constant throughout the body**. Deep
moonquake nests at 700–1200 km see the same tidal tensor as the surface. Only
latitude and longitude enter, and only to set the local frame. (Elastic response
does vary with depth, but that is the Love-number scale factor, which does not
change timing.)

*Coulomb stress is linear in the tensor*, `ΔCFS = T_ij C_ij` with
`C_ij = ½(uᵢnⱼ + uⱼnᵢ) + μ nᵢnⱼ`. So the variance of ΔCFS over any time series is a
quadratic form `cᵀΣc` on the **covariance of the tensor components**, precomputed
once. A grid search over thousands of planes and thousands of epochs collapses from
a full pass per plane to ~72 flops per plane — exactly, not approximately. This is
what made a 2000-trial null over 74 nests run in a minute.

### The result — and a null that lied

Weber's criterion: search fault orientations for the one minimising
`score = std(ΔCFS at events) / std(ΔCFS over the span)`. Search space 1,728 planes ×
5 friction values. Nests with ≥20 events: **74**.

| Null construction | Nests below p = 0.05 |
|---|---|
| **Uniform** random event times | **73 / 74** |
| **Shift** — slide the sequence, preserve spacing | **17 / 74** |

**The uniform null was wrong and would have produced a spectacular false
positive.** Drawing random times destroys the catalogue's temporal clustering, so
the test silently became "are events clustered in time at all?" — trivially yes,
since events at nearby epochs see nearly identical tensors and therefore low
variance, with no tidal alignment required.

The shift null preserves relative spacing and breaks only the alignment with the
forcing, which is the question actually intended.

**This is the second time a null has manufactured significance in this project**,
after the pooled-phase test (A3/A4). Both had the same root cause: *a null that
fails to preserve the catalogue's own structure tests a different hypothesis than
the one stated.*

### Honest final numbers

```text
shift null, 2000 trials:   17/74 nests nominally p < 0.05   (3.7 expected)
Benjamini-Hochberg, FDR 0.05:   0/74 nests survive
ensemble excess:                17 vs 3.7, sd 1.87  ->  7.1 sigma
```

**No individual nest is defensible after multiple-comparison correction. The
population-level excess is 7.1σ.**

That is a real but modest result, and it matches Weber, Bills & Johnson (2009)
directly: they report that for some clusters the constant-stress fit is good while
for others it is "not strongly dependent on plane orientation." A heterogeneous
population with a minority of well-constrained nests is precisely what 17/74 looks
like.

Note also that the p-value floor bit in practice. BH over 74 tests needs
p ≤ 6.8×10⁻⁴ for the strongest nest to survive; at 200 trials the floor was
5×10⁻³, so **no nest could have passed regardless of signal strength.** Raising to
2000 trials moved the floor to 5×10⁻⁴, below the threshold — after which the answer
was still zero, but now for a real reason rather than an arithmetic one.

### Caveats

1. The shift null may itself be partially degenerate, since catalogue and forcing
   share the anomalistic month (A3/A4). 17/74 could still be inflated.
2. No Love numbers, so this is stress **shape**, not magnitude. Timing and the
   orientation search are unaffected; absolute ΔCFS in Pa is not available and must
   not be reported from this code.
3. Grid resolution is 15°/15°/30°. A finer grid would lower every score, observed
   and null alike, so the comparison holds — but per-nest orientations should not
   be quoted at better than grid resolution.

### Plan status

**A5 complete. A6 complete, with a negative headline and a positive ensemble.**
Phase 1 has now produced its real deliverable, which is not "tides drive
moonquakes" — that was known — but a **validated instrument plus two documented
ways to fool it**. Both traps are waiting in Phase 3, where earthquake catalogues
carry aftershock clustering and seasonal detection cycles of their own.

---

## 2026-08-22 — C1/D1: Parkfield LFEs, and the M2/S2 gate fires

**Code:** `ph-core::parkfield`, `examples/parkfield_constituents.rs`. 48 tests.

### The dataset

Shelly (2017, updated 2024), from USGS ScienceBase — **public domain**.
**1,528,117 LFEs, 88 families, 23.1 years.** Every family holds ≥5,333 events;
median 16,058, max 44,156. Compare the largest Apollo moonquake nest at 85.

LFE **families** are the direct analogue of moonquake **nests**. Most families
individually clear Beeler & Lockner's ~10⁴ event requirement.

### The measurement

Schuster power evaluated at the **exact** periods of named constituents rather
than by blind peak search, which turns doc 08 §13b's validity gate into a direct
measurement. Pooled catalogue:

| Constituent | Period (d) | Power | vs S1 |
|---|---|---|---|
| M2 principal lunar | 0.5175 | 1,831 | 0.11 |
| **S2 principal solar** | 0.5000 | **7,729** | 0.48 |
| N2 larger elliptic | 0.5274 | 21 | 0.00 |
| **K1 luni-solar diurnal** | 0.9973 | **18,861** | **1.16** |
| **S1 solar diurnal** | 1.0000 | **16,245** | 1.00 |
| O1 principal lunar diurnal | 1.0758 | 8,052 | 0.50 |
| Mf lunar fortnightly | 13.661 | 334 | 0.02 |
| Msf lunar synodic fortnightly | 14.765 | 663 | 0.04 |
| Mm lunar monthly | 27.555 | 631 | 0.04 |
| Ssa solar semiannual | 182.62 | 431 | 0.03 |
| Sa solar annual | 365.26 | 4,654 | 0.29 |

Null expectation for `D²/N` is **1**. The largest family alone reproduces every
ratio to within a few percent.

### The gate fires

**S1 is exactly 24.000 h and has essentially no body-tide amplitude.** Its power of
16,245 is therefore a *direct measurement of the detection artifact floor* — LFE
detection is template matching on continuous data, and its sensitivity tracks the
day-night cultural noise cycle.

Three readings follow immediately:

1. **K1 (23.93 h) at 1.16× S1 is unusable**, exactly as doc 08 §13b predicted.
   Degenerate with the diurnal artifact.
2. **S2 (7,729) exceeds M2 (1,831) by 4.2×.** For a genuine body tide the ordering
   must be the other way — the M2 tidal potential is roughly 2.2× S2's. Seeing S2
   dominate is the signature of a **thermal or cultural semidiurnal cycle**, not a
   body tide. **The gate returns "suspect artifact."**
3. **Sa (4,654) is large but confounded** with seasonal variation in detection.

### What this does and does not mean

**It does not contradict the literature.** Thomas et al. (2012) find strong tidal
modulation of Parkfield LFEs — but they resolve *tidal stress on the fault* and
test its phase. They do not fold event times on trial periods.

**The lesson is sharper than "the catalogue is noisy":**

> **Raw period folding on a detection-limited catalogue measures the detector, not
> the Earth.**

Our own doc 02 argued this from theory. It is now measured, on a catalogue chosen
because the physical effect there is known to be strong.

**This is the third trap the machinery has caught**, after the pooled-phase null
(A3/A4) and the uniform-time null (A5/A6). All three share a shape: *a statistic
that looks decisive while silently answering a different question.*

### Consequence for the band prediction

**This run does not test it.** The long-period constituents (Mf 334, Msf 663,
Mm 631, Ssa 431) are all far below the diurnal artifacts — but that comparison is
meaningless while the diurnal band is artifact-dominated. Testing the 1 yr–200 yr
band prediction requires the artifact removed first.

### Next (C2)

Do it properly, as Thomas et al. do: compute the tidal tensor at the Parkfield
family locations, resolve ΔCFS onto San Andreas geometry (`ph-core::fault` already
does this), and test phase against a structure-preserving null. The full
exponential response `R₀exp(S_T/Aσ₀)` with `M > 1` is required here — Thomas et al.
infer `Aσ₀ = 6×10⁻⁴ MPa`, putting Parkfield firmly in the non-linear regime.

Also worth adding: an explicit **alias analysis** (doc 08 §13e). With a 24 h
detection cycle this strong, beats against it will appear elsewhere in the
spectrum and need blacklisting.

---

## 2026-08-22 — C2: first result that survives correction

**Code:** `examples/parkfield_coulomb_phase.rs`.

### The null, fixed before the run

Every trap so far came from choosing a null after seeing the data, so this one was
specified first and justified on structure.

D1 established that the catalogue carries a large detection artifact locked to
**solar time** (S1 power 16,245 against a null expectation of 1). So:

> **Null: shift every event time by a whole number of solar days.**

- The artifact is locked to local solar time → a whole-day shift leaves it
  **exactly invariant**.
- The lunar tide precesses ~50 min per solar day → a whole-day shift **does**
  decorrelate ΔCFS phase.

It holds the confound fixed while sliding only the quantity under test. 400 shifts
drawn from ±(30–4000) days, excluding near-multiples of the 27.55 d anomalistic and
29.53 d synodic months.

### Setup

ΔCFS from the real Earth tidal tensor (Moon, Sun, planets) in `IAU_EARTH`, rotated
to local NED at 35.635 N, −120.150 E, resolved onto deep San Andreas geometry —
strike 137°, dip 90°, rake 180°, μ = 0.4. Sampled at 0.02 d over 23 years plus
padding: **841,272 points**.

**Self-check passed:** mean interval between ΔCFS maxima came out **0.5175 d** —
exactly M2. The forcing is what we think it is.

### Result

| | |
|---|---|
| Families tested (largest 12) | 12 |
| Nominally p < 0.05 | **10 / 12** (0.6 expected) |
| **Surviving Benjamini-Hochberg, FDR 0.05** | **9 / 12** |
| `D²/N` range | 106 – 1702 (null expectation 1) |
| Preferred phases | **all within 71.4°**, spanning 0.7° to 72.1° |

**The null has demonstrable power.** Two families are *not* significant, with
p = 0.14 and p = 0.22, and observed `D²` sits near — not far above — the null
maximum for most families. That spread is the signature of a working null,
in direct contrast to A5/A6's 73/74 pinned at the floor.

Every family clusters **shortly after the ΔCFS maximum**, consistent with Thomas
et al. (2012) for these same LFEs.

Geometry is robust: `D²/N` moves only 404 → 415 across strike ±10° and dip 80–90°.
A thrust geometry gives both lower power (356) and a very different phase
(−108.8°), so the statistic does discriminate between orientations.

### Caveats — none of these are small

1. **The 12 families are co-located within ~30 km**, so they see nearly identical
   forcing. Phase coherence across them is a **coherence check, not independent
   confirmation**. It shows the pipeline returns a consistent answer; it is not 12
   independent tests of phase.
2. **Absolute phase is geometry-dependent.** Strike 127 → 137 → 147 shifts the
   preferred phase 57° → 41° → 25°. The *existence* of clustering is robust; the
   *value* should not be quoted without stating the assumed geometry.
3. **Effect sizes are modest.** For most families observed `D²` is at roughly the
   98th–99th percentile of the null, not orders above it.
4. **No Love numbers**, so this is stress *shape*. Timing is unaffected; absolute
   ΔCFS in Pa is unavailable.
5. **This does not test the band prediction.** M2 is a 12.42 h constituent, deep
   inside the range doc 07 predicts should be *damped* for ordinary crust. That it
   works at Parkfield is consistent with tremor's `T_a` sitting in the tidal band —
   but the 1 yr–200 yr prediction remains **entirely untested**.

### Where this leaves the project

Phase 1 delivered a validated instrument and three documented traps. **C2 is the
first time the instrument has been pointed at contested data with a
pre-specified null and returned a result that survives correction.**

That is the methodological milestone the validation ladder existed to reach: not
"we found a correlation," but "we found one *after* three separate opportunities to
fool ourselves, using a null fixed in advance."

---

## 2026-08-22 — C3: the amplitude law, and a fourth trap

**Code:** `ph_core::phase::cycle_amplitude_at`, `examples/parkfield_amplitude_law.rs`.
50 tests.

### Why amplitude rather than constituent

D1 showed the diurnal and semidiurnal bands are dominated by a detection artifact
locked to solar time, so a constituent-by-constituent spectrum runs straight back
into it.

But **the artifact does not care how strong the tide is.** Binning events by the
peak-to-trough amplitude of the ΔCFS cycle they fall in gives a test that detection
bias cannot produce, because it is amplitude-independent by construction.

This is a direct test of the amplitude law `R̃/r = Δτ/(aσ̄)` (Ader et al. 2014, eq.
B7). Predictions fixed before the run: flat → artifact; `D²/N ∝ amplitude²` →
linear response; steeper → the non-linear exponential regime.

### Result

All 1,528,117 events, six equal-count amplitude bins:

| Bin | Events | Amplitude (rel.) | `D²/N` | ratio to bin 1 |
|---|---|---|---|---|
| 1 | 254,686 | 1.00 | 216 | 1.00 |
| 2 | 254,686 | 1.48 | 2,301 | 10.67 |
| 3 | 254,686 | 1.85 | 3,513 | 16.30 |
| 4 | 254,686 | 2.20 | 3,883 | 18.01 |
| 5 | 254,686 | 2.53 | 7,639 | 35.44 |
| 6 | 254,687 | 3.00 | 15,288 | 70.91 |

**Monotonic across every bin. A 3× amplitude range produces a 71× change in phase
concentration.** Log-log slope **3.56**, against 2 for a linear response.

### The fourth trap — per-bin nulls are compromised

Bin 6 has the **highest** `D²/N` (15,288) and yet **p = 0.45**.

That is the tell. **Binning by amplitude selects on the forcing itself**, and
high-amplitude cycles recur at the ~14.77 d spring–neap period, so every bin
inherits temporal structure derived from the very signal under test. A whole-day
shift does not preserve that structure, so the per-bin null answers a different
question — again.

**The per-bin p-values in the table above are not trustworthy and should not be
quoted.**

### Nulling the claim actually being made

The claim is the *trend*, so the trend is what must be nulled: shift the event
times, **re-bin at the shifted times**, refit the slope.

```text
observed slope 3.56    null median 0.66    null max 3.53    p = 0.0050
```

**The trend survives.** Null slopes centre near 0.66 — no trend, as expected when
alignment is broken — while the observed 3.56 exceeds all 200 draws.

⚠ **But only just.** Null max 3.53 against observed 3.56 is a thin margin, and
p = 0.0050 is the 200-trial floor rather than a measurement. The true value is
plausibly somewhat higher. More trials would sharpen it; the null slope
distribution has a heavy tail and deserves more sampling before this is published.

### Reading

Response scales with tidal forcing amplitude, **faster than linearly**. That is
consistent with the exponential form `R = R₀exp(S_T/Aσ₀)` being in its non-linear
regime — which is exactly what Thomas et al. (2012) imply for Parkfield with
`Aσ₀ = 6×10⁻⁴ MPa`, giving `S_T/Aσ₀ ≈ 0.2–2`.

Two independent lines now agree that Parkfield LFEs respond to tidal Coulomb
stress: C2's phase clustering (9/12 families surviving FDR) and C3's amplitude
scaling. They rest on different statistics and different nulls.

### Trap count: four

1. Pooled-phase null degenerate when catalogue shares the forcing's period
2. Uniform-time null tests temporal clustering, not tidal alignment
3. Raw period folding measures the detector
4. **Per-bin nulls compromised when the binning variable derives from the forcing**

All four share one shape: *a statistic that looks decisive while silently answering
a different question.* Each was caught only because the result was checked against
what the null was actually testing, rather than against intuition.

### Still untested

**The band prediction.** Every result so far sits at M2 (12.42 h) — deep inside the
range doc 07 predicts should be *damped* for ordinary crust. Parkfield works
because tremor's `T_a` sits in the tidal band. The 1 yr–200 yr claim for ordinary
crust remains untouched, and needs a terrestrial earthquake catalogue (Phase 3).

---

## 2026-08-22 — C3b confirms; C4 fails, and the failure is instructive

### C3b — the amplitude law holds at 2,000 trials

```text
observed slope 3.56   null median 0.43   null max 4.09   p = 0.0095
```

Against the 200-trial run (p = 0.0050, null max 3.53). As flagged, **that p was a
floor, not a measurement** — with ten times the sampling, 18 of 2,000 null draws
now exceed the observed slope and the null maximum overtakes it.

**The result survives at p ≈ 0.01** rather than 0.005. Weaker than it first looked,
still significant, and now properly estimated. The prediction that the floor was
hiding the true value was correct.

Also required an optimisation worth noting: sorting 1.5M events per trial would
have dominated runtime, so amplitude bin edges are now fixed once from the observed
distribution and reused. Binning is O(n) and — more importantly — **identical
between observed and null**, which is what a like-for-like comparison needs.

### C4 — frequency-resolved response: method failed

Built `ph_core::demod` (complex demodulation, 5 tests, verified to separate M2 from
a 4× larger S2 when the window exceeds their 14.77 d beat). Then attempted response
per constituent at Parkfield.

**Two nulls tried, both invalid.**

**Trap 5 — the time-shift null, recurring.** Demodulation isolates one constituent,
which makes the band a **near-pure tone** — and against a pure tone a time shift
merely *rotates* the phase cluster without diluting it. `D²` is near-invariant and
the null has no power.

This is trap 1 exactly. `ph_core::stats` documents it in a module-level warning.
**The demodulation step recreated the degenerate case, and I walked into it
anyway.** Worth recording plainly: documenting a trap is not the same as being
immune to it, because the trap reappears wearing different clothes.

**Trap 6 — the sham-frequency null.** Replacing it with "run the same procedure at
frequencies carrying no tide" seemed sound. The result:

```text
sham floor (10 quiet bands):  median D²/N = 1,439,329   max 1,520,587
M2:  17,986      O1: 8,357      Mf: 4,195      Sa: 1,488,398
```

With N = 1,528,117, a sham median of 1.44M means **D² ≈ N²** — events land at
essentially *one* phase at frequencies where there is no tide. Every real
constituent sits **below** the floor.

The diagnosis: at a frequency with no genuine power, `z̄` is dominated by whatever
leaks in from the nearest strong constituent at ω′. Then `z̄ ∝ e^{i(ω′−ω)t}`, so
`arg z̄ = (ω′−ω)t`, and the reported band phase `ωt + arg z̄` collapses to `ω′t` —
**the phase of the leaking constituent, not of the target.** A tide-free frequency
does not give a neutral baseline; it gives a relabelled copy of the dominant band.

### What C4 needs instead

**Do not derive constituent phase by demodulating the composite ΔCFS series.**
Compute it from the **astronomical argument** — the exact Doodson combination of
fundamental arguments for that constituent. That phase is uniform over long spans
by construction, so events uniform in time are uniform in phase, and both the
statistic and its null behave.

This is what the tidal literature does, and it is why the field works in terms of
constituent arguments rather than filtered series.

**It also makes D2 a blocker.** Analytic constituent arguments mean the HW95/KSM03
expansion, which doc 07's fourth pass left as "free access not yet established."
D2 moves from parallel task to **prerequisite for C4**.

### Status

- **C4 is not done.** No frequency-resolved response exists yet, and no transfer
  function. The numbers above are artifacts and must not be quoted.
- `ph_core::demod` is correct and tested for what it does — isolating a constituent
  from a series. The error was using its *phase output* as a statistic against a
  null, not the demodulation itself. It stays.
- **C2 and C3 are unaffected.** Both used the full composite ΔCFS, which is
  genuinely quasi-periodic, so their nulls retain power. C3b independently
  confirmed at 2,000 trials.

### Trap count: six

1. Pooled-phase null degenerate when catalogue shares the forcing's period
2. Uniform-time null tests temporal clustering, not tidal alignment
3. Raw period folding measures the detector
4. Per-bin nulls compromised when the binning variable derives from the forcing
5. **Time-shift null degenerate against a demodulated single constituent** — trap 1
   recurring, walked into despite being documented
6. **Sham-frequency null invalid**: demodulating where there is no tide returns the
   leaking constituent's phase, not a neutral baseline

Six for six, all the same shape: *a statistic that looks decisive while silently
answering a different question.*
