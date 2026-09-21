# Gompertz-like hazard from load×gain coupling of damaged subsystems: a computational biogerontology object

**Thesis #13** (series label NP-07). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Under what network-coupling assumptions does near-linear local damage produce Gompertz-like hazard in a simulated subsystem graph, and which of those assumptions are identifiable from public demographic schedules versus remaining free gain parameters?

The calculations use six abstract nodes. Local damage is exactly linear when the damage-rate gain is zero. An exponential map of mean load is an exact Gompertz hazard (slope 0.0840 per model year, doubling time 8.25 years). The same damage under an additive map that matches the endpoints is not Gompertz-like (weighted residual 5776 against a χ² budget of 83.68). On a synthetic age schedule the hazard-map gain is free: its profile is flat, and only the product of gain and load velocity is closed. A mean-load channel identifies that gain. Damage-rate coupling is a curvature parameter with a closed profile. Ring and star graphs reproduce the uncoupled hazard once the gain (and, for the star, the baseline hazard) may move.

The draws are a synthetic surrogate, not a download of the Human Mortality Database. Related work on systemic amplification and on network Gompertz laws is cited as related work. This deposit does not claim organism rejuvenation, plasma-exchange efficacy, or clinical age reversal.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/identifiability.py` | Seeded damage, hazard, Fisher and profile sketches (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Damage, coupling scan, spectra, profiles, and graphs |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/identifiability.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the sketches. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Gompertz-like hazard from load×gain coupling of damaged subsystems: a computational biogerontology object [Internet]. Thesis #13 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-13-gompertz-load-gain-coupling

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
