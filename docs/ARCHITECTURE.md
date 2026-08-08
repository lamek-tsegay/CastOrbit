# ARCHITECTURE.md — CastOrbit

Standing brief for this repository. Read this first, every session.
Read `docs/V2_BRIEF.md` **before** this file — it defines the current scope and
wins wherever the two disagree about *what the project is*.
Read `PHYSICS.md` before writing any simulation code; it wins wherever
anything disagrees about *how the physics works*. `docs/PHYSICS_V2.md`
documents the equations V2 adds (adaptive stepping, the climatology
atmosphere, disposal delta-v) and does not modify anything in `PHYSICS.md`.

---

## 1. What this is

A satellite design tool built on a validated orbital-decay engine.

**The one-line pitch:**
> Describe a satellite and a mission. CastOrbit sizes it, builds it, flies it,
> and tells you whether it complies.

**The physics does not change; the scope of what it is pointed at does.** The
engine in `sim/` — hand-written RK4, NRLMSIS atmosphere, the `PHYSICS.md` §3.2
equation — is the same engine that reproduced the February 2022 Starlink loss.
V2 generalises what it is aimed at: arbitrary altitude, inclination, epoch and
target, with ram area *derived from geometry* rather than assumed.

That last point is the substance of V2, not a feature. The README's central
finding is that `rho`, `Cd` and `A` are inseparable from a decay curve — only
the product is observable. Pinning `A` from geometry the user specified
constrains one of the three by construction rather than by assumption
(`V2_BRIEF.md` §2).

### V1 is the credibility, not legacy code

V1 was a reconstruction of one event: did the February 2022 batch survive?
That is complete, validated, and tagged `v1.0-validated`. **The Baruah
reproduction is now a regression test.** It stays in the test suite, it stays
in the README, and any change that moves it outside its acceptance band stops
work until understood — that is a bug, not a scope change.

It is still **not** a satellite tracker or a 3D visualisation of public TLE
data. It is now deliberately a design tool, which V1 explicitly was not — but
a design tool whose numbers come from an engine that was checked against a
peer-reviewed result, which is the distinction that matters. Above ~550 km
"will it decay?" stops being the interesting question and the compliance
question replaces it (`V2_BRIEF.md` §3).

### Success criteria

This project succeeds if a skeptical aerospace engineer can:

1. Ask "did you write the integrator?" — and the answer is yes, with tests
2. Ask "how accurate is it?" — and get a number, not a shrug
3. Ask "what happens if the drag area is wrong?" — and see it swept
4. Ask "where did this generated design's mass come from?" — and get a real
   spacecraft it was interpolated from, not a plausible number
5. Find the limitations section before they have to point one out

It fails if it is visually impressive and numerically unverifiable. **When
trading polish against verifiability, choose verifiability.** The specific
failure mode V2 must avoid: a studio UI over plausible numbers with no engine
underneath (`V2_BRIEF.md` §8).

---

## 2. Architecture

```
Python (physics)  →  JSON  →  React (playback only)
```

**The frontend never simulates anything.** No physics in JavaScript. Not for
convenience, not for interactivity, not for "just this one slider." Python emits
precomputed runs; the browser plays them back and draws charts.

This is the single most important structural rule in the project. It is what
makes a three-day build feasible, and it keeps all physics in one testable place.

### When the user moves a slider

The parameter space is precomputed. Sliders select among existing runs and
interpolate between them for display. If a combination truly has not been
computed, show that clearly rather than faking a result.

---

## 3. Repository layout

```
/
├── PHYSICS.md              equations — authoritative, read before coding
├── README.md               written last
│
├── sim/                    Python. All physics lives here.
│   ├── constants.py        MU, R_E, G0
│   ├── atmosphere.py       pymsis wrapper, SW-All.csv loader, density cache
│   ├── integrator.py       RK4, hand-written
│   ├── dynamics.py         the da/dt and dm/dt equations
│   ├── critical.py         critical density and critical altitude solver
│   ├── satellite.py        config dataclasses, ThrusterMode enum
│   ├── mission.py          V2 entry point: arbitrary altitude/inclination/epoch
│   ├── disposal.py         V2 delta-v, propellant, compliance verdict
│   ├── montecarlo.py       batch runs, outcome taxonomy
│   ├── sweeps.py           insertion altitude, ram area, safe-mode timing
│   ├── export.py           JSON emission for the frontend
│   └── validate.py         the four tests from PHYSICS.md §8
│
├── tests/                  pytest. Tests 1–3 are unit tests.
│
├── data/
│   ├── SW-All.csv                  CelesTrak space weather indices
│   ├── satellite_specs.json        hardware, each field source + confidence tagged
│   ├── disposal_rules.json         post-mission disposal rules, each cited
│   ├── event_feb2022.json          validation target parameters
│   └── starlink_snapshot.txt       TLE snapshot, offline use
│
├── out/                    generated JSON — gitignored, regenerable
│
├── docs/
│   ├── V2_BRIEF.md         current scope -- read before this file
│   ├── ARCHITECTURE.md     this file
│   ├── PHYSICS_V2.md       equations added in V2
│   ├── SOURCES.md          every external source cited
│   └── baruah_2024.pdf     the validation paper
│
└── web/                    React + Vite frontend
    ├── src/
    └── package.json
```

---

## 4. Stack

**Python:** numpy, scipy (interpolation only — *not* for the main integrator),
pymsis, pytest. Keep it minimal.

**Frontend:** React + Vite, `react-globe.gl` for the globe, `recharts` for plots.

Do not add dependencies beyond these without a clear reason. Every added package
is a build failure waiting to happen on day three.

### Explicitly forbidden

- `scipy.integrate.solve_ivp` for the main propagation loop. Write RK4 by hand.
  This is the first thing anyone will ask about.
- Physics in JavaScript. See §2.
- CesiumJS. Heavier than needed.
- Any orbital mechanics library (poliastro, hapsira, orekit). The model is one
  ODE; a library adds a day of dependency debugging and removes the thing that
  makes this yours.
- 6DOF anything.
- Tuning parameters to make validation match. Report the discrepancy instead.

---

## 5. Build order

Each phase has a gate. **Do not start the next phase until the gate passes.**

**Phases 1–6 below are complete** — V1, tagged `v1.0-validated`. They are kept
here as the record of what was built and what each gate actually required, not
as work remaining. **Phases 7–12 are the current plan and live in
[`V2_BRIEF.md`](V2_BRIEF.md) §7**, which is authoritative for them; they are
not duplicated here.

### Phase 1 — Physics core ✅

Implement `constants`, `dynamics`, `integrator`, `atmosphere`, `critical`.

**Gate:** Validation tests 1, 2, and 3 from `PHYSICS.md` §8 pass. The critical
density table in `PHYSICS.md` §4.1 is reproduced exactly.

### Phase 2 — Validation ✅

Implement the Baruah et al. reproduction (test 4). Run at `Cd = 1.0`.

**Gate:** A decay curve exists for both the 1.00 m² and 4.48 m² cases, with the
altitude at 08:58 UT on 5 Feb 2022 printed and compared against the published
targets. Whatever the numbers are, they are recorded.

If they disagree badly, **investigate before proceeding** — but do not adjust
parameters to force agreement. A documented 30% discrepancy with a hypothesis
about its cause is a stronger result than a suspicious exact match.

### Phase 3 — Batches and sweeps ✅

Monte Carlo over 49 satellites. The three sweeps in `PHYSICS.md` §9.

**Gate:** Survival-fraction-versus-insertion-altitude curves exist for both quiet
and storm conditions, as matplotlib plots. No frontend yet.

### Phase 4 — JSON export ✅

Define and emit the run format. See §6.

**Gate:** A JSON file exists in `out/` that fully describes a batch, and a Python
script can read it back and reproduce the plots from Phase 3.

### Phase 5 — Frontend ✅

React app. Globe, altitude chart, sweep chart, validation view.

**Gate:** It loads the JSON and displays it. Nothing is hardcoded.

### Phase 6 — Writeup ✅

`README.md` with the equations, the validation table, and the limitations list
from `PHYSICS.md` §10.

**This phase is not optional.** The writeup is the artifact that travels; the
demo is what people look at for ninety seconds. If time runs short, cut frontend
polish, not the writeup.

---

## 6. JSON contract

Two files, `out/batch.json` and `out/sweeps.json`, emitted by
`sim/export.py`. Deliberately separate rather than one stretched schema — a
Monte Carlo batch and a parameter sweep are different shapes of data and
forcing them into one file just makes both harder to consume.

### `out/batch.json`

The fleet reproduction: both Cd conventions compared throughout the README's
central finding, run side by side rather than picking one. Every satellite
carries its own `cd_times_area_m2` — the one parameter that's actually
observable from a decay curve (see the README) — rather than leaving a
frontend to recompute `cd * ram_area_m2` itself.

```json
{
  "meta": {
    "scenario": "fleet_reproduction_feb2022",
    "generated": "ISO-8601",
    "sim_version": "git sha",
    "atmosphere_model": "NRLMSIS 2.1 via pymsis 0.12.0",
    "epoch": "2022-02-03T17:43:00+00:00",
    "window_end": "2022-02-08T00:00:00+00:00",
    "n_satellites": 49,
    "insertion_altitude_km": 210.0,
    "reentry_altitude_km": 100.0
  },
  "observed": { "lost": 38, "survived": 11, "source": "data/event_feb2022.json" },
  "runs": [
    {
      "label": "cd2.2",
      "description": "this project's own baseline Cd",
      "config": {
        "cd": 2.2, "density_scale": 1.0,
        "ram_area_range_m2": [1.00, 4.48],
        "effective_drag_range_m2": [2.20, 9.86]
      },
      "outcome_counts": { "REENTERED": 49, "INDETERMINATE": 0, "...": 0 },
      "critical_altitude": {
        "_note": "counterfactual -- see sim/export.py",
        "times": ["ISO-8601", "..."],
        "h_crit_km": [181.3, 180.9, "..."]
      },
      "satellites": [
        {
          "id": 0,
          "params": {
            "mass_kg": 227.0, "ram_area_m2": 2.1, "cd": 2.2,
            "cd_times_area_m2": 4.62, "thrust_n": 0.0
          },
          "outcome": "REENTERED",
          "outcome_time": "ISO-8601",
          "trajectory": {
            "t_s": [0, 600, 1200],
            "h_km": [210.0, 209.4, 208.7],
            "rho": [1.47e-10, 1.51e-10, 1.55e-10]
          }
        }
      ]
    },
    { "label": "cd1.0", "...": "same shape, Baruah et al.'s convention" }
  ]
}
```

### `out/sweeps.json`

The three §9 sweep curves plus the analytic critical-altitude band each is
checked against. See `sim/sweeps.py`'s `plot_from_payload` for the exact keys
— it plots from this schema alone, which is also how the Phase 4 read-back
check works (below).

Downsample trajectories for export — 10 minute spacing is plenty for display.
Keep full resolution in Python.

**Include `critical_altitude` as a time series.** Drawing the critical altitude
as a moving line on the altitude chart, with trajectories crossing it and then
diverging, is the single clearest visual in the project. For a fleet in safe
mode it's a counterfactual (F = 0 the whole window has no balance point) —
label it as one, as `sim/export.py` does.

**The read-back check is only real if the plotting code cannot see the
in-memory objects.** `sim/sweeps.py` plots from a plain payload dict, and
`sim/export.py`'s `replot_sweeps_from_json` / `replot_batch_from_json` feed it
one parsed straight off disk — `tests/test_export.py` exercises both against
freshly-written JSON, not a stale file sitting in `out/`.

---

## 7. Frontend spec

Dark, dense, information-first. The reference aesthetic is satellitemap.space:
technical, unornamented, every pixel carrying data.

**Four views:**

1. **Globe** — 49 satellites, coloured by outcome. Click one for its parameters
   and cause of loss. This is the arresting image.
2. **Altitude chart** — 49 curves over time, with the critical altitude line
   overlaid. Survivors bend up, losses bend down. This is the *explanatory*
   image, and it matters more than the globe.
3. **Sweeps** — survival fraction versus insertion altitude, quiet versus storm.
   This is the engineering argument.
4. **Validation** — the Baruah comparison table, discrepancies shown plainly.
   This is the credibility.

Playback scrubber with a UT clock, shared across views.

Colour by outcome, consistently everywhere:

| Outcome | Colour |
|---|---|
| `REACHED_SHELL` | green |
| `REENTERED` | red |
| `PROPELLANT_EXHAUSTED` | amber |
| `INDETERMINATE` | grey |

**Do not spend time on:** shaders, custom satellite meshes, camera animations,
loading screens, transitions. Every hour there is an hour not spent on physics.

---

## 8. Two traps

Both are documented in `PHYSICS.md`. They are repeated here because they will
each cost hours if missed.

**Cd = 1.0 for validation.** Baruah et al. use unity, stated explicitly as a
simplification. Running validation at 2.2 doubles the drag and the model
reenters far too early. Validate at 1.0; use 2.2 for the project's own runs and
report both.

**Safe mode is why they died.** A model with thrusters fighting drag predicts
most satellites survive — at 210 km the thrust term generally beats the drag
term. SpaceX commanded the fleet into a low-drag attitude with thrusters off.
Thruster state is an explicit time-dependent input, never a constant.

---

## 9. Working style

- Small commits, each leaving tests passing.
- Every physics function gets a docstring naming its `PHYSICS.md` section.
- When a published number is used, cite it inline.
- When something is a guess, say so in the code and in the writeup.
- If a validation result is disappointing, report it. The honest version is more
  persuasive than the flattering one, and it is the version that survives being
  questioned.

**When time runs short, cut in this order:** globe polish → globe → extra
sweeps → Monte Carlo size. Never cut: the four validation tests, the Baruah
comparison, the limitations section.

---

## 10. Environment notes (added in Phase 1)

The system Python on this machine is 3.9, which is too old for current `pymsis`
wheels. The virtualenv is built against Python 3.12 via `uv`:

```bash
uv venv --python python3.12 .venv
uv pip install --python .venv/bin/python pymsis==0.12.0 numpy scipy pytest matplotlib
.venv/bin/python -m pytest tests/ -q
```

`pymsis` is pinned to 0.12.0 — the version the `PHYSICS.md` §4.1 reference
values were generated with.
