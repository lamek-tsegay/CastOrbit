# SOURCES.md — CastOrbit

Every external source this project draws numbers, models, or data from. Where
`data/satellite_specs.json` already tags a value's provenance and confidence,
this file points to that tag rather than repeating it — the JSON is the record
of record for hardware figures; this file is the record of record for
everything else, and the index that ties both together.

---

## The validation target

**Baruah, Y., Roy, S., Sinha, S., Palmerio, E., Pal, S., Oliveira, D. M., &
Nandy, D. (2024).** *The Loss of Starlink Satellites in February 2022: How
Moderate Geomagnetic Storms Can Adversely Affect Assets in Low‑Earth Orbit.*
Space Weather, 22, e2023SW003716. DOI:
[10.1029/2023SW003716](https://doi.org/10.1029/2023SW003716)

The paper this whole project validates against. Event parameters, the two
bounding-case decay targets, the Swarm C comparison, and the Cd = 1.0
validation convention all come from here — see `data/event_feb2022.json` for
the exact values extracted and `PHYSICS.md §8` / `README.md` for how they're
used. Author list verified against the paper's NSF Public Access Repository
record, independent of the publisher page.

**Corroborating loss-count sources** (`data/event_feb2022.json` →
`corroborating_sources`), used only to show that published counts vary and 38
is the figure this project adopts, not to re-derive anything:

- Kataoka, R., et al. (2022). *J. Space Weather Space Clim.*, 12, 41. — 38 lost.
- Zhang, S.-R., et al. (2022). *Space Weather*, 20, e2022SW003168. — 40 lost.
- Guarnieri, F., et al. (2023). arXiv:2307.02923. — 32 lost (NORAD tracking analysis).

---

## The atmosphere model

**Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M. Jr.,
Mlynczak, M. G., et al. (2021).** *NRLMSIS 2.0: A whole‑atmosphere empirical
model of temperature and neutral species densities.* Earth and Space Science,
8(3), e2020EA001321. DOI:
[10.1029/2020EA001321](https://doi.org/10.1029/2020EA001321)

**Emmert, J. T., Jones, M. Jr., Siskind, D. E., Drob, D. P., Picone, J. M.,
Stevens, M. H., et al. (2022).** *NRLMSIS 2.1: An empirical model of nitric
oxide incorporated into MSIS.* Journal of Geophysical Research: Space Physics,
127, e2022JA030896. DOI:
[10.1029/2022JA030896](https://doi.org/10.1029/2022JA030896)

The density model underlying every `rho` in this simulator
(`sim/atmosphere.py`). This project calls it through **pymsis 0.12.0**
(pinned in `requirements.txt` — the version the `PHYSICS.md §4.1` reference
table was generated against, not floated):

**Lucas, G. (2022).** *pymsis* [Computer software]. DOI:
[10.5281/zenodo.5348502](https://doi.org/10.5281/zenodo.5348502).
Source: [github.com/SWxTREC/pymsis](https://github.com/SWxTREC/pymsis).
Docs: [swxtrec.github.io/pymsis](https://swxtrec.github.io/pymsis/).

---

## Space weather indices

**CelesTrak `SW-All.csv`** — daily and 3‑hourly geomagnetic (Ap/ap, Kp) and
solar flux (F10.7) indices, 1957‑10‑01 onward. Committed at
`data/SW-All.csv` so results are reproducible against a fixed input rather
than whatever CelesTrak serves on a given day.
Source: [celestrak.org/SpaceData/SW-All.csv](https://celestrak.org/SpaceData/SW-All.csv).

Methodology reference for the underlying index compilation: Vallado, D. A.,
and T. S. Kelso, *"Using EOP and Solar Weather Data for Satellite
Operations,"* 15th AIAA/AAS Astrodynamics Specialist Conference, Lake Tahoe,
CA, 2005 August 7–11.
([celestrak.org/publications/AAS/05-406/](https://celestrak.org/publications/AAS/05-406/))

CelesTrak is maintained by Dr. T. S. Kelso.

**CelesTrak Starlink GP/TLE data** — `data/starlink_snapshot.txt`, the
current active-constellation element set for Phase 5's globe view. Not tied to
any date the validation depends on; fetched per the `curl` recipe in
`setup.sh`.
Source: [celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle](https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle).

---

## Satellite hardware

Every field in **`data/satellite_specs.json`** carries its own `source` and
`confidence` tag (`published` / `derived` / `estimated` / `disputed`) —
that file is authoritative for hardware figures, not this list. Summarized
here for a single index of everything the project cites:

| Field | Confidence | Source, as tagged |
|---|---|---|
| v1.5 mass, knife-edge area, max area, Cd (validation) | published/estimated/disputed | Baruah et al. 2024, §5 (above) |
| v1.5 thruster type | published | SpaceX; *SpaceNews*, 2023‑02‑28 |
| v1.5 thrust (0.071 N), Isp (1666 s) | **derived** | Computed from SpaceX's published *relative* scaling claims against the v2 mini (2.4× thrust, 1.5× Isp) — not directly published; labelled DERIVED everywhere it appears in code and in `README.md` |
| v1.5/v2 mini bus dimensions, v1.5/v2 mini nominal ram area, v2 mini mass | estimated/disputed | Secondary reporting, not independently verified by this project — `satellite_specs.json`'s own `action` fields note where sources disagree by 4–10× and instruct sweeping rather than picking a value |
| v2 mini thruster type, thrust, Isp, power, solar array, bus dimensions | published | SpaceX |
| Drag coefficient baseline (Cd = 2.2) | estimated | Standard convention for satellites in free-molecular flow; used in NASA DAS, STK, and orbit-determination work generally, not a Starlink-specific figure |
| Swarm C reference case | published | Baruah et al. 2024, §5 (above) |

The "secondary reporting" and "Source A" / "Source B" entries are exactly as
vague as `satellite_specs.json` itself records them — the file's own
`_comment` and `action` fields explain why (typically a several-fold spread
between two uncredited secondary sources). No more specific citation exists to
give without inventing one; sweeping the disputed range, which
`sim/montecarlo.py` does, is the response to that gap, not resolving it here.

---

## Numerical method

RK4 (`sim/integrator.py`) and the atmosphere/critical-altitude bisection
(`sim/critical.py`) are standard textbook methods, hand-implemented per
`PHYSICS.md §7`, `§4` — no external library or paper is cited for these
because none is used; `scipy.integrate.solve_ivp` is explicitly excluded
(`docs/ARCHITECTURE.md §4`).
