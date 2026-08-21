# 08 — Generated Hypotheses

Ideas produced by the first literature pass. These are candidate research
contributions, not established results. Each states the idea, why the literature
supports it, and what would have to be true for it to fail.

---

## 1. The harmonic spectrum as an inverse method for fault friction ★

**The idea.**

Beeler & Lockner showed that a fault's correlation with periodic stressing depends
on whether the forcing period exceeds the **characteristic nucleation time**. That
is not just a caveat — it is a statement that a fault acts as a **low-pass filter**
on periodic stress.

So a rate-and-state fault has a **transfer function**: amplitude response R(ω) and
phase lag φ(ω) to sinusoidal forcing. Which means the harmonic spectrum we are
already computing (doc 02) should show:

- Strong response at low harmonic orders / long periods
- A **knee** at the nucleation frequency
- Suppressed response above it

And crucially — **the location of the knee estimates the nucleation time, hence
in-situ frictional properties (aσ, D_c).**

**Why this is the most valuable idea here.**

It converts the multi-basis harmonic framework from *exploratory feature
engineering* into a **physically-grounded inverse problem**. We would no longer be
asking "does base N correlate?" but "what does the shape of the response spectrum
tell us about the fault?"

That is a real geophysical measurement, publishable independently of any
forecasting result, and it makes every harmonic order physically meaningful.

**Failure condition.** If the observed spectrum is flat, or the knee is not
resolvable given catalogue size, the inversion is unconstrained. Power analysis
(doc 04 §7) determines this before we invest.

**Deep dive needed:** derive R(ω) and φ(ω) analytically. See doc 09 §1.

---

## 2. Phase-versus-frequency as an anti-artifact test ★

**The idea.**

Rate-and-state predicts not just *whether* seismicity responds, but a **specific
phase lag that varies with frequency**. That is a two-dimensional prediction:
a curve φ(ω), not a single number.

A spurious correlation — multiple testing, catalogue bias, hydrological aliasing —
has **no reason whatsoever to reproduce the predicted lag-versus-frequency
relation.** It might accidentally produce a significant amplitude at one order. It
will not accidentally trace a theoretical phase curve across many orders.

**Why this matters.**

It is a far stronger test than amplitude significance, and it is essentially free
once the Fourier encoding is in place — the phase is already in the coefficients.

This directly answers the Holt & Newman failure mode: they tested amplitude-like
significance across thousands of intervals. Testing *structure across orders*
rather than *significance at any order* is qualitatively harder to fake.

**Recommend making this the project's primary evidentiary standard.**

**Failure condition.** Requires enough events per harmonic order to estimate phase
with useful uncertainty. Phase estimates are noisier than amplitude estimates.

---

## 3. Use the tidal potential catalogue as the prior over harmonic orders ★

**The idea.**

We have been treating "which base numbers matter" as an open search. But the
**Cartwright–Tayler–Edden** tidal potential catalogue already tells us the true
amplitude of the astronomical forcing at each frequency — several hundred
constituents with measured coefficients.

So: **weight the prior over harmonic orders by the actual tidal potential
amplitude at that frequency.**

**Why this is powerful.**

1. It collapses the search space from "all bases" to "physically forced
   frequencies," which massively reduces the multiple-testing burden — the single
   biggest statistical threat to the project.
2. It gives every feature a physical amplitude, so effect sizes become
   *predictions* rather than free parameters. A constituent with 10× the forcing
   amplitude should show a larger response; that ratio is testable.
3. It preserves the exploratory goal: if an **unforced** base lights up beyond its
   prior, that is precisely the anomaly worth investigating — and now it stands out
   against a principled background instead of being one of 28,782 tests.

This reframes "not just base 12" from a scan into a **structured comparison
against known forcing.**

**Failure condition.** If real response is dominated by non-tidal periodicities
(hydrological, catalogue artifacts), the prior misleads. Mitigate by fitting the
unforced orders too and comparing.

---

## 4. The sensitivity field β(x,t) as a tomographic inverse problem

**The idea.**

Beaucé's result is regional and retrospective. Generalised, **β(x,t)** — the
coupling between tidal stress and seismicity rate — is a continuous field to be
estimated from sparse point observations (earthquakes). That is structurally a
**tomography / spatial inverse problem**: Gaussian process or spline field with a
smoothing prior, estimated inside the point-process likelihood.

**The valuable part: resolution.**

Inverse problems come with **resolution kernels**. We can compute *where we have
enough events to constrain β and where we do not*, and publish that map alongside
the estimate.

An honest map of "here is where we can say something, and here is where we cannot"
is more scientifically valuable — and far more credible — than a global field with
uniform false confidence. It also directly addresses the objection that celestial
correlation studies over-claim spatial coverage.

**Failure condition.** β may vary faster than the event rate can track, making the
field under-determined everywhere. Resolution analysis reveals this honestly
rather than hiding it.

---

## 5. Model the detection function instead of truncating at Mc

**The idea.**

Standard practice truncates the catalogue at the magnitude of completeness,
discarding most events. Instead, put a **detection probability function
P_det(m, x, t)** directly into the point-process likelihood.

This is routine in astronomy (selection functions) and under-used in seismology.

**Why it matters here.**

Beeler & Lockner say we need 10⁵–10⁶ events. Mc truncation throws away the
majority of the catalogue — exactly the small, shallow events where Métivier found
the anomaly is **largest**. Modelling detection instead of truncating could
recover an order of magnitude more data in the regime where the signal is
strongest.

It also converts the Mc-drift problem from a confounder into a modelled quantity.
Detection improved over the catalogue span; unmodelled, that trend aliases into
long-period features and manufactures signal (doc 04 §6b). Modelled, it is handled.

**Failure condition.** A misspecified detection function introduces bias worse than
truncation. Validate against Mc-truncated results on subsets where both are viable.

---

## 6. Deep moonquakes yield lunar interior structure as a by-product

Beyond validating the pipeline (doc 05 §1), deep moonquake tidal analysis is an
established route to constraining **lunar Love numbers and deep interior
structure**.

So Phase 1 of the validation ladder is not merely a test harness — it can produce
an independent scientific result. That matters strategically: it establishes
methodological credibility in a domain where tidal forcing is uncontroversial,
*before* any terrestrial claim is made.

---

## 7. Pipeline self-test with synthetic null features

**The idea.**

Before trusting any positive result, construct a **synthetic feature with identical
spectral and autocorrelation properties to a real celestial feature but no physical
link to seismicity**, and run the complete pipeline on it.

The pipeline must return null. If it does not, the pipeline is broken — and we
learn that cheaply, before publishing.

Extend to a small suite: shuffled-phase features, features from a fictitious body
on a plausible orbit, features from a real body evaluated at wrong epochs.

This is cheap, decisive, and almost never done in this literature. It should be a
CI test, not a one-off.

---

## 8. Ocean tidal loading is the compute bottleneck, and probably the signal

Cochran et al.'s factor-of-3 result was driven by **ocean tidal loading**, not
solid Earth tide. Loading is a **spatial convolution** of ocean tide height with
Farrell Green's functions — far more expensive than the five-number degree-2 solid
tide, and unavoidable if we want the settings where the effect is largest.

Design implication: the architecture in doc 06 optimises the *cheap* part. The
expensive part is loading. It needs its own precompute strategy — likely
constituent-wise loading coefficients per location, computed once and stored, so
runtime cost reduces to harmonic synthesis.

Flagged as the most likely source of unpleasant surprise in the implementation
schedule.

---

## Priority for next pass

| # | Idea | Value | Cost | Priority |
|---|---|---|---|---|
| 1 | Rate-and-state spectrum inversion | Very high | Medium | **1** |
| 2 | Phase-vs-frequency test | Very high | Low | **1** |
| 3 | CTE amplitude prior | High | Low | **2** |
| 7 | Synthetic null self-test | High | Low | **2** |
| 5 | Detection function | High | Medium | 3 |
| 4 | β tomography + resolution | High | High | 3 |
| 8 | Ocean loading precompute | Required | High | 3 |
| 6 | Moonquake interior inversion | Medium | Medium | 4 |
| 9 | β = τ̇/(aσ) synthesis | **Very high** | Medium | **1** — partly confirmed |
| 10 | PTA statistical machinery | High | Low | **1** — confirmed viable |
| 11 | Hydrological forcing as second probe | **Very high** | Low | **1** |

Ideas 1, 2 and 3 together would convert this project from a correlation search
into a physically-grounded inverse problem with a built-in artifact test. That is
the strategic goal of the next research pass.

Idea 9 subsumes and sharpens them: if β ∝ τ̇/(aσ) holds, the project acquires an
**external validation channel** (β vs. GPS strain rate) that does not route through
earthquake forecasting at all — and a falsifiable spectral-slope prediction that
distinguishes real response from artifact. Verify or kill it first.

---

## 9. Synthesis: β measures τ̇/(aσ), and that unifies the literature ★★

*Added during the brainstorm pass. Speculative — the arithmetic below is
back-of-envelope and agenda item 1 exists to verify or kill it.*

**The chain.**

Dieterich's characteristic time is

```
t_a = a σ / τ̇
```

with σ the **effective** normal stress and τ̇ the local stressing rate. If the
seismicity response to periodic forcing is a low-pass filter with corner at
ω ≈ 1/t_a, then tidal sensitivity is governed by t_a — and therefore

```
β  ∝  1 / t_a  =  τ̇ / (a σ)
```

**Tidal sensitivity is a measurement of local stressing rate divided by effective
normal stress.**

**Why this is interesting: it explains four separate observations at once.**

| Observation | Explanation via t_a |
|---|---|
| Earthquakes correlate weakly with tides (Beeler & Lockner) | Typical crust: a≈0.01, σ≈150 MPa, τ̇≈0.003 MPa/yr → t_a ~ centuries. Tidal periods sit far above the corner; response strongly suppressed. |
| Subduction thrust faults respond strongly (Cochran) | High pore pressure → small effective σ → shorter t_a → less suppression. |
| Tremor responds far more than earthquakes (Rubinstein) | Near-lithostatic pore pressure *and* transiently enormous τ̇ during slow-slip episodes → t_a collapses → corner frequency enters the tidal band. |
| Sensitivity rose before Ridgecrest (Beaucé) | β rising = τ̇ rising and/or σ_eff falling = approaching failure. Exactly what a precursor should be. |

**Two sharp, falsifiable predictions.**

1. ~~**Spectral slope.**~~ **SUPERSEDED — see revision note below.** The original
   1/ω monotonic-ranking prediction was tested against the literature in the second
   research pass and does not hold.

2. **Spatial correlation.** β should correlate with independent estimates of pore
   pressure and geodetic strain rate. Both are measurable (GPS/InSAR strain,
   seismic velocity ratios, hydrological models) — so β can be validated against
   data that has nothing to do with celestial mechanics.

Prediction 2 is the important one. It provides an **external validation channel**
for the whole framework that does not route through earthquake forecasting at all.
If β correlates with GPS strain rate, the quantity is real regardless of whether
it improves forecasts.

**Consequence for the "all bases" question.**

This gives a *physical* answer to which harmonic orders should matter: those where
the tidal potential has genuine amplitude (CTE/HW95, hypothesis 3) **and** the
period exceeds the nucleation time (this hypothesis). That selection likely favours
the long-period constituents.

Worth stating carefully: those timescales — fortnightly, monthly, annual, 18.6-year
— overlap the ones traditional astrology emphasises. That is a **coincidence of
timescale arising from the same orbital mechanics**, not evidence for interpretive
claims. It explains why both systems attend to the same periods without implying
anything about meaning.

**Consequence for Resonant Finder / Star Seer.**

If the fault is a filter, then "resonance points" should be defined on the
**forced response**, not the forcing. Response = forcing convolved with the
transfer function, which shifts and smooths peaks relative to exact geometric
alignment. **Search the response, not the stimulus.** This is a concrete
architectural correction for those applications.

**Failure conditions.** The linear small-perturbation treatment may not apply; the
Beeler–Lockner second response mode may dominate; t_a may be unresolvable given
catalogue sizes; spectral structure may be masked by catalogue artifacts (Mc drift,
hydrological aliasing). The synthetic-null self-test (hypothesis 7) addresses the
last.

---

### Revision note — 2026-08-21, second research pass

**Prediction 1 was wrong and is withdrawn.** Ader et al. (2014) and Heimisson &
Avouac (2020) give a **two-regime** response, not a simple low-pass:

```
t_a = a σ / τ̇                    T_c ≈ 2π t_a

T ≪ t_a :  seismicity rate tracks stress          amplitude ∝ A      (flat in T)
T ≫ t_a :  seismicity rate tracks stressing rate  amplitude ∝ A / T  (falls with T)
```

Long-period constituents are therefore **not** automatically favoured, and the
monotonic ranking (nodal > Sa > Ssa > Mm > Mf > semidiurnal) is dead.

**Replacement prediction — the response spectrum has a peak near T ≈ t_a.**
A band-pass shape requires two timescales: t_a governing the long-period roll-off,
and a separate nucleation duration governing the short-period cutoff. If so the
spectrum identifies **two** fault properties, and a **peak at a predicted location**
is far more diagnostic than a slope — artifacts readily produce slopes, but rarely
peaks where theory says one should be.

⚠ Heimisson & Avouac state their model is **not valid near T ≈ t_a** — exactly the
region of interest. The band-pass reconstruction is inferred from abstracts, not
from the papers. **Obtain the primary sources before building on this.**

**Prediction 2 (β vs. GPS strain rate / pore pressure) survives and is
strengthened.** EPSL (2025) finds tidal modulation amplitude for San Andreas LFEs
is controlled mainly by **background effective stress** — direct published support
for β ∝ 1/(aσ).

**Priority consequence.** Hypothesis 1 (spectrum as inverse method) is **already an
active research programme** — Beeler et al. (2018), EPSL (2025). It is established
literature we can cite rather than a claim we must defend, but it is not ours.
Novelty must move to: **long-period constituents** (published work is
diurnal/semidiurnal only), **global scale** (published work is single-fault), the
**generalised Doodson extension**, and **other bodies**.

---

## 10. Import statistical machinery from pulsar timing arrays

Detecting a tiny correlated periodic signal buried in red noise, with rigorous
false-alarm control, is the **central problem of pulsar timing array analysis** —
and the PTA community has spent two decades building machinery for exactly it:

- Red-noise-aware Bayesian model comparison
- Marginalisation over unknown noise spectral indices
- Correlated-noise-aware significance estimation
- Rigorous, well-tested false-alarm control on periodic detections

The seismological tidal-triggering literature generally uses the Schuster test,
which assumes independent samples and has no red-noise treatment. Earthquake
catalogues are strongly red (aftershock clustering) — precisely the condition PTA
methods were built for.

**This looks like an under-exploited transfer.** The statistical problem shape is
the same; the field is different, so the tooling has not crossed over. Worth a
dedicated session in the next research pass.

Related transferable domains: asteroseismic mode identification (power spectra with
correlated noise), and gravitational-wave matched filtering (known-waveform search
with rigorous detection statistics — relevant if the transfer function gives us a
predicted response waveform to match against).

---

## 11. Use hydrological loading as a second probe, not just a confounder ★★

*Added second pass. Supersedes the framing of hydrological loading as a pure
threat.*

**The reframe.**

Doc 07 flagged seasonal hydrological loading as the project's most dangerous
confounder, because it is annual and so are Sa and Ssa. That framing was
incomplete.

Hydrological loading is a **second periodic forcing with independently measurable
amplitude and phase** — from GRACE/GRACE-FO gravity, GLDAS land surface models,
well levels, and snow load. We do not have to guess it. We can measure it.

Which means it is not only a nuisance to subtract. It is an **instrument**.

**Why this is powerful.**

The transfer function R(ω), φ(ω) is the object we actually want (hypothesis 1).
A single forcing probes it at one set of frequencies. **Two forcings with
independently known amplitudes probe it at two well-separated frequency bands** —
tidal (hours to weeks) and hydrological (annual) — which constrains the shape far
better than either alone.

Critically, the annual band is where tidal amplitude is weak and hydrological
amplitude is strong. They are complementary rather than competing.

**Published precedent.** *Tidal and hydrological seismicity modulations reveal pore
fluid diffusion during earthquake nucleation* (Science Advances, sciadv.ady6350)
does exactly this joint treatment to infer pore fluid diffusion. So the approach is
validated; our extension is scale and constituent coverage.

**It also resolves the confounder cleanly.** Modelling both forcings jointly with
measured hydrological amplitude means an annual signal gets **attributed**, not
merely flagged as ambiguous. The confounder stops being a threat to validity and
becomes a fitted term.

**Failure condition.** Hydrological loading models have real uncertainty, and
GRACE resolution (~300 km) is coarse relative to fault-scale processes. If the
hydrological amplitude is badly wrong, it will contaminate rather than calibrate.

**Priority: raise to 1.** Low cost, high value, and it converts the project's
largest identified threat into an asset.
