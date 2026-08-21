# 11 — Literature Access Without Institutional Credentials

Project constraint: **no institutional access.** Everything must be obtained free
and legitimately. This is a working record of what actually yields results.

## What worked

**Author self-archived pages — highest yield by far.**
Both critical-path papers (Ader et al. 2014 *GJI*, Heimisson & Avouac 2020 *GRL*)
were obtained from `web.gps.caltech.edu/~avouac/publications/`. Most active
researchers maintain a publications page with PDFs; journal self-archiving policies
generally permit it.

Method: search `"<author surname>" publications <institution>`, or
`"<paper title>" filetype:pdf`.

Note: the Caltech host has an incomplete TLS chain — browser fetch tools may fail
on certificate verification where `curl -sL` succeeds.

## Routes worth trying, roughly in order

| Route | Best for | Notes |
|---|---|---|
| **Author / lab pages** | Anything | Highest yield. Try first. |
| **USGS Publications Warehouse** (`pubs.usgs.gov`) | Beeler, Cochran, Hardebeck, Michael | **USGS-authored work is public domain.** Covers a large share of this field. |
| **arXiv** | Theory, statistics, PTA methods | Full text, no barrier |
| **ESS Open Archive** (`essopenarchive.org`) | **AGU journals — GRL, JGR** | AGU's own preprint server. Directly relevant, often overlooked. |
| **EarthArXiv** (`eartharxiv.org`) | Geoscience preprints | Growing coverage |
| **Unpaywall API** (`api.unpaywall.org/v2/{doi}?email=`) | Any DOI | Returns legal OA locations. Free, scriptable. |
| **OpenAlex** (`api.openalex.org`) | Any DOI | Metadata plus OA links. Free, no key. |
| **Semantic Scholar API** | Any DOI | `openAccessPdf` field. Free. |
| **CORE** (`core.ac.uk`) | Repository aggregation | Large index |
| **NASA ADS** (`ui.adsabs.harvard.edu`) | Astronomy, geophysics | Often links free full text |
| **HAL** (`hal.science`) | French institutions | IPGP, IGN — relevant for Métivier et al. |
| **ResearchGate** | Anything | Author-uploaded; hit or miss, often abstract only |
| **PubMed Central** | Anything NIH-adjacent | Rare here |
| **Institutional repositories** | Theses | Ader's PhD thesis likely contains the 2014 GJI material in fuller form |

## Particularly high-value for this project

**USGS Publications Warehouse.** Beeler & Lockner (2003), Beeler et al. (2018), and
much of the Coulomb-stress and rate-and-state literature is USGS-authored and
therefore **public domain**. This may be the single most valuable route for our
critical path.

**ESS Open Archive.** AGU preprints — GRL and JGR are where much of the tidal
triggering literature lives.

**PhD theses.** Often contain a fuller derivation than the published paper, with
appendices the journal cut. Ader (Caltech), Heimisson (Stanford), and Beaucé are
all worth checking. Caltech THESIS and Stanford's repository are both open.

## Data — all free

None of the datasets in doc 09 §10 require institutional access:

- **USGS ComCat** — public API
- **GCMT** — public
- **IRIS/EarthScope** — public, open waveform and catalogue services
- **IERS EOP** — public
- **GRACE / GRACE-FO, GLDAS** — NASA, public
- **Apollo PSE** — NASA PDS Geosciences Node, public
- **FES2014** — free for research with registration; **TPXO** free for academic use;
  **GOT** NASA, free
- **DE440 / DE441 SPICE kernels** — NAIF, public

Data access is not a constraint on this project. Only paywalled *papers* are, and
mostly they are obtainable by the routes above.

## Standing practice

1. Try author page first, then USGS/ESS/arXiv, then Unpaywall by DOI.
2. When a paper is obtained, **record the working URL in
   [bibliography.md](bibliography.md)** so it need not be rediscovered.
3. When a paper cannot be obtained, note that explicitly in the research log rather
   than reconstructing its content from abstracts. The second pass shows why —
   an abstract-based reconstruction was directionally right and structurally wrong.
