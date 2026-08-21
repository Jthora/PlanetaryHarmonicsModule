# Bibliography

Tiered per [00-framing.md](00-framing.md). Tier A = established physics,
B = established method / novel application, C = exploratory or contested.

---

## Tidal triggering — positive results

**Tanaka, S., Ohtake, M. & Sato, H. (2002).** Evidence for tidal triggering of
earthquakes as revealed from statistical analysis of global data. *JGR: Solid
Earth* 107(B10), 2211. — Tier A. Establishes the Schuster-test methodology.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001JB001577

**Cochran, E.S., Vidale, J.E. & Tanaka, S. (2004).** Earth tides can trigger
shallow thrust fault earthquakes. *Science* 306(5699), 1164–1166. — Tier A.
Rate ×3; p < 10⁻⁴; best fit μ = 0.4. Strongest positive result in the field.
https://www.science.org/doi/10.1126/science.1103961

**Métivier, L. et al. (2009).** Evidence of earthquake triggering by the solid
earth tides. *EPSL* 278(3–4), 370–375. — Tier A. 442,412 NEIC events, ~99%
confidence, uplift phase.
https://www.ipgp.fr/~lalmetiv/metivier_etal_epsl2009.pdf

**Scholz, C.H., Tan, Y.J. & Albino, F. (2019).** The mechanism of tidal triggering
of earthquakes at mid-ocean ridges. *Nature Communications* 10, 2526. — Tier A.
Inverted phase; magma chamber mechanism.
https://www.nature.com/articles/s41467-019-10605-2

**Ide, S., Yabe, S. & Tanaka, Y. (2016).** Earthquake potential revealed by tidal
influence on earthquake size–frequency statistics. *Nature Geoscience* 9, 834–837.
— Tier A/B. b-value modulation rather than rate modulation.
https://www.nature.com/articles/ngeo2796

**Beaucé, E. et al. (2023).** Enhanced tidal sensitivity of seismicity before the
2019 magnitude 7.1 Ridgecrest, California earthquake. *GRL* 50, e2023GL104375. —
Tier A/B. Sensitivity rise ~1.5 yr prior; basis for the β(x,t) target.
https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL104375

---

## Mechanism and theory

**Beeler, N.M. & Lockner, D.A. (2003).** Why earthquakes correlate weakly with the
solid Earth tides: effects of periodic stress on the rate and probability of
earthquake occurrence. *JGR: Solid Earth* 108(B8), 2391. — Tier A. **Most important
methodological paper for this project.** Nucleation-time argument; 10⁵–10⁶ event
requirement.
https://www.researchgate.net/publication/228542429

**Dieterich, J. (1994).** A constitutive law for rate of earthquake production and
its application to earthquake clustering. *JGR* 99(B2), 2601–2618. — Tier A.
Source of the rate-and-state seismicity equation. *To obtain — doc 09 §1.*

**King, G.C.P., Stein, R.S. & Lin, J. (1994).** Static stress changes and the
triggering of earthquakes. *BSSA* 84(3). — Tier A. Coulomb failure stress
conventions. *To obtain — doc 09 §5.*

---

## Tremor and slow slip — high sensitivity

**Rubinstein, J.L. et al. (2008).** Tidal modulation of nonvolcanic tremor.
*Science* 319(5860), 186–189. — Tier A. Tremor pulsing at 12.4 h and 24–25 h.
https://www.science.org/doi/abs/10.1126/science.1150558

**Thomas, A.M. et al. (2013).** Evidence for tidal triggering of high-amplitude
rapid tremor reversals and tremor streaks in northern Cascadia. *GRL* 40(16). —
Tier A.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50832

**Thomas, A.M. et al. (2012).** Tidal triggering of low frequency earthquakes near
Parkfield, California. *JGR* 117(B5). — Tier A. Also flags spurious-correlation
risk from other stress components.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB009036

**Hawthorne, J.C. & Rubin, A.M. (2010).** Tidal modulation of slow slip in
Cascadia. *JGR* 115(B9). — Tier A.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JB007502

---

## Negative, null, and contested

**Southern California study.** *GJI* 205(2), 681. Tidal stress triggering of
earthquakes in Southern California — no significant correlation found. Tier A.
https://academic.oup.com/gji/article/205/2/681/685563

**Bendick, R. & Bilham, R. (2017).** Do weak global stresses synchronize
earthquakes? *GRL* 44, 8320–8327. — Tier C. 32-yr LOD/M7+ cycle; questioned on
effective cycle count.
https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017GL074934

**Holt, E.W. & Newman, E. (2025).** Tidal triggering of magnitude 7+ earthquakes by
Jupiter. arXiv:2508.07064. — **Not a supporting citation.** Retained as the
project's canonical multiple-comparisons cautionary example (doc 00).
https://arxiv.org/abs/2508.07064

**Tidally synchronized solar dynamo: a rebuttal.** arXiv:2206.14809. — Tier C.
Cite alongside Scafetta / Abreu et al. if solar-planetary coupling is discussed at
all.
https://arxiv.org/pdf/2206.14809

---

## Forecasting methodology and tooling

**RECAST** — neural temporal point process, GRU encoder-decoder; matches or beats
ETAS on Southern California given sufficient catalogue size.
https://github.com/keliankaz/recast

**EarthquakeNPP** (2024). A benchmark for earthquake forecasting with neural point
processes. arXiv:2410.08226.
https://arxiv.org/pdf/2410.08226

**pyCSEP / CSEP** — N-test, S-test, M-test, L/PL-test; T-test of equal predictive
ability via information gain per earthquake.
https://cseptesting.org/

**Ogata, Y. (1988, 1998).** ETAS formulation and space–time extension. *To obtain —
doc 09 §7.*

---

## Reference works to obtain

- **IERS Conventions (2010), Ch. 7** — solid Earth tide, Love numbers,
  frequency-dependent corrections, pole tide. Authoritative. *doc 09 §3*
- **Farrell, W.E. (1972).** Deformation of the Earth by surface loads.
  *Rev. Geophys.* 10(3). Load Love numbers and Green's functions. *doc 09 §4*
- **Cartwright & Tayler (1971); Cartwright & Edden (1973).** Tidal potential
  harmonic development. *doc 09 §2*
- **Hartmann, T. & Wenzel, H.-G. (1995).** HW95 tidal potential catalogue — more
  complete than CTE, includes planetary terms. *doc 09 §2*
- **Agnew, D.C.** SPOTL: Some Programs for Ocean-Tide Loading. *doc 09 §4*
- **Schuster, A. (1897).** On lunar and solar periodicities of earthquakes.
  Origin of the Schuster test. *doc 09 §6*
- **Nakamura, Y. et al.** Apollo Passive Seismic Experiment deep moonquake
  catalogue. *doc 09 §10*

---

## Data sources

| Source | Contents |
|---|---|
| USGS ComCat | Global earthquake catalogue |
| GCMT | Focal mechanisms (strike/dip/rake) |
| IERS EOP | Polar motion, length of day |
| FES2014 / TPXO9 / GOT | Global ocean tide models |
| GRACE / GRACE-FO, GLDAS | Hydrological loading (annual confounder) |
| Apollo PSE | Deep moonquake catalogue |
| InSight | Marsquake catalogue |
