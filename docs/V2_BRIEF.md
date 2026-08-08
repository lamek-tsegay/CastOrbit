# V2_BRIEF.md — CastOrbit, second scope

**Read this before `docs/ARCHITECTURE.md`.** Where the two disagree about
*what the project is*, this document wins. Where they disagree about *how the
physics works*, `PHYSICS.md` wins and always has.

---

## 1. What changed

V1 was a reconstruction: one event, one satellite type, one question — did the
February 2022 batch survive? That is complete, validated, and tagged
`v1.0-validated`.

V2 is a design tool built on the same engine:

> **Describe a satellite and a mission. CastOrbit sizes it, builds it,
> flies it, and tells you whether it complies.**

The physics does not change. The scope of what it is pointed at does.

### What V1 becomes

| V1 asset | Role in V2 |
|---|---|
| `sim/` — integrator, atmosphere, dynamics | The engine. Unchanged core. |
| Baruah reproduction | **The credibility.** Never breaks. |
| Monte Carlo, sweeps | Generalised to arbitrary designs |
| Frontend | Becomes the Validation view inside a larger studio |

**The Feb 2022 reproduction is not legacy code.** It is the reason anyone
should believe a generated lifecycle. It stays in the test suite, it stays in
the README, and any change that breaks it is a bug, not a scope change.

---

## 2. The V2 loop

```
description  ->  design spec  ->  geometry  ->  ram area
                      |                             |
                      v                             v
                 mass budget  --------------->  lifecycle sim
                                                    |
                                                    v
                                              compliance verdict
```

Each arrow is computed. None is generated, guessed, or hardcoded.

**The critical link is geometry -> ram area.** In V1, `A` was an assumed
parameter swept across a disputed range. In V2 it is *derived* from geometry
the user specified, projected onto the velocity vector. That directly attacks
the degeneracy documented in the README: rho, Cd and A are inseparable from a
decay curve, but if geometry pins A, one of the three is constrained by
construction rather than assumed.

This is the difference between a design tool and a mesh generator with a
spec table beside it.

---

## 3. The altitude problem, and the job it creates

The V1 engine was built for 200-300 km, where drag dominates and decay takes
days. Generated designs will mostly live higher, where drag barely acts.

Natural decay to 100 km, Bc = m/(Cd*A) = 100 kg/m^2, quiet conditions,
computed with this project's own atmosphere model:

| Altitude | da/dt | Time to reentry |
|---|---|---|
| 300 km | -211 km/yr | ~1 year |
| 400 km | -23.6 km/yr | ~13 years |
| 500 km | -3.5 km/yr | ~100+ years |
| 705 km | -0.2 km/yr | millennia |

(Linear estimates. Real decay accelerates, so these overstate — but the order
of magnitude holds.)

**So "will it decay?" is the wrong question above ~550 km.** The right one is
regulatory:

> Post-mission disposal rules require LEO satellites to reenter within a
> defined window (the FCC adopted a five-year rule; the long-standing
> international guideline was 25 years). **Verify the current rule and its
> applicability before implementing — do not take this paragraph as
> authoritative.**

That reframes the lifecycle sim as a compliance calculator:

- Below ~550 km: does natural decay comply? Compute the time.
- Above: it cannot. Compute the **delta-v to lower perigee enough that it
  does**, convert to propellant mass via the rocket equation, and check that
  against the design's own propellant budget.

A design that cannot afford its own disposal is a real and common failure,
and this tool can catch it. That is the compliance verdict.

---

## 4. Propagator constraint — read before Phase 7

RK4 at 10 s steps over 25 years is roughly 8e7 steps per satellite. Infeasible.

But `da/dt` in `PHYSICS.md` §3.2 is already an orbit-averaged secular rate, not
an instantaneous one. It is smooth and slow when drag is weak, so the step can
scale with the physics:

- **Adaptive stepping.** Choose `dt` such that the fractional change in `a`
  per step stays under a fixed tolerance (1e-4 is a reasonable start). At
  700 km this permits steps of days; at 200 km it will collapse to seconds
  on its own.
- **Do not replace RK4.** Same integrator, variable step. "Did you write the
  integrator?" must keep its answer.
- **Regression requirement.** The adaptive propagator must reproduce the
  Baruah results to within 0.1% of the fixed-step values. That test is the
  gate for this work.

---

## 5. Mass closure — scope honestly

Real mass estimating relationships are weeks of work. Do not start there.

**Stage 1 (build now):** interpolate from a table of published satellites —
mass, power, array area, bus dimensions across smallsat to large-GEO. Tag
every output `estimated` with its source, exactly as `data/satellite_specs.json`
already does. Interpolation between real spacecraft is defensible and honest.

**Stage 2 (later):** component-level MERs with a convergence loop.

**Never:** an LLM producing a mass number. The LLM's only job is turning prose
into a spec the sizing code consumes.

A generated design must carry its provenance. If dry mass came from
interpolating two real satellites, the UI says so.

---

## 6. Rules carried forward from V1

These are unchanged and non-negotiable.

- **No physics in JavaScript.** Python computes, the browser displays.
- **Nothing hardcoded, nothing faked.** Every displayed number traces to an
  engine output. If a view needs a value, add it to the export.
- **No `solve_ivp`** in the propagation loop.
- **Do not tune to match.** Report discrepancies with a hypothesis.
- **Uncertain values get swept, not chosen.**
- **Commit granularly.** Every commit compiles; every test-bearing commit passes.
- **Verify UI in a real browser**, not by reading JSX.

### New rule

- **The Baruah reproduction is a regression test.** Any change that moves it
  outside its acceptance band stops work until understood.

---

## 7. Phase plan

Same discipline: each phase has a gate, no phase starts before the previous
one passes.

### Phase 7 — Generalise the engine
Arbitrary altitude, inclination, epoch, target. Adaptive stepping per §4.
Delta-v and propellant-mass calculation for disposal.

**Gate:** Baruah reproduction still passes within 0.1% of V1 fixed-step
values. A 705 km design propagates 25 years in reasonable wall time. Disposal
delta-v computed and checked against a hand calculation.

### Phase 8 — Geometry and projected area
Parametric bus + panels from a spec. Projected-area solver: geometry plus
attitude gives ram area.

**Gate:** For a Starlink v1.5-like spec, the solver's knife-edge and broadside
areas bracket the published 1.00-4.48 m^2 range. Analytic check: a flat plate
at angle theta to the velocity vector gives `A*cos(theta)` to within 1%.

### Phase 9 — Mass closure, stage 1
Interpolation table, provenance tags.

**Gate:** Reproduces the dry mass of three real satellites held out of the
table, within 25%. Every output tagged with its source.

### Phase 10 — LLM front door
Prose to spec, constrained JSON. Smallest phase in the project.

**Gate:** Ten varied descriptions produce schema-valid specs. Invalid or
underspecified input asks a question rather than inventing values.

### Phase 11 — Studio UI
Three-panel layout per `docs/mockup/`. Validation view survives intact.

**Gate:** Every number traces to an engine field. Provenance visible in UI.

### Phase 12 — Writeup
README updated. Validation stays the credibility section.

---

## 8. Realistic scoping

Phases 7 and 8 are the ones that matter and the ones that make the tool honest.
9 through 11 can ship in a rougher form.

**If time is short, build 7 and 8 properly and stop.** A generalised, validated
engine with a real projected-area solver is a coherent thing to show. A
half-finished mass-closure loop underneath a pretty studio is not.

Phase 11 without Phases 7-9 is the Lovable mockup: plausible numbers, no
engine. That is the failure mode this document exists to prevent.
