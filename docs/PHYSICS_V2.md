# PHYSICS_V2.md — equations added in V2

Companion to [`PHYSICS.md`](../PHYSICS.md), which is **unchanged and remains
authoritative**. Every equation V1 used is still specified there and still
governs: the state is `[a, m]`, the core equation is §3.2, the integrator is
§7, and the limitations in §10 all still apply.

This file documents only what V2 adds, and the reason it is a separate file is
that none of it changes how the existing physics works. Where the two files
could be read as disagreeing, `PHYSICS.md` wins.

---

## V2.1 Adaptive step control

The equation and the integrator are unchanged. `rk4_step` is the same
classical four-stage RK4; only the choice of `dt` is new.

`da/dt` in `PHYSICS.md` §3.2 is an orbit-averaged secular rate. It is smooth
and slow when drag is weak, so the step can scale with the physics:

```
dt = tol * a / |da/dt|
```

clamped to `[dt_min, dt_max]`, with `tol = 1e-4` by default.

**This is a rate bound, not a truncation-error estimate**, and the distinction
matters. It bounds how far the state moves per step, which for this ODE is the
quantity that governs accuracy: the physics is driven by `rho(h)`, and the way
this integration goes wrong is a step that carries the satellite a long way
down the density profile before the derivative is re-evaluated. RK4's own
truncation error is far below the rate bound — verified by halving `tol` and
confirming the answer does not move
(`tests/test_adaptive.py::test_tightening_tolerance_does_not_move_the_answer`).

Two properties are worth recording because they are not obvious:

- **The controller is self-limiting.** Bounding `|Δa/a|` also bounds how much
  `rho` can change within the step, so the start-of-step rate estimate stays
  accurate. Measured `max_frac_change` is ~1.03e-4 against a 1e-4 tolerance on
  both the exponential stand-in and real NRLMSIS, with zero rejections. The
  reject-and-halve path is a genuine safety net but is almost never
  load-bearing.
- **`dt_max` matters as much as `tol`.** At 700 km the tolerance alone permits
  steps of years. The cap is what keeps the integration resolving the
  atmosphere it is flying through, and it must be set below the timescale on
  which the density model varies — 1800 s on a real `DensityGrid`, the grid's
  own node spacing.

**Validation.** Against the `PHYSICS.md` §8 Test 4 configuration and real space
weather, the adaptive driver reproduces the V1 fixed-step results to 0.0012%
(reentry time, 4.48 m²) and 0.0001% (altitude at reference, 1.00 m²), using
89–180× fewer steps.

---

## V2.2 Long-horizon atmosphere

`PHYSICS.md` §6 specifies the atmosphere for a five-day reconstruction with
real indices. A 25-year compliance run cannot use it, and the binding reason is
not cost:

> `data/SW-All.csv` ends a few weeks after the present. A disposal run reaching
> 2051 spends ~99% of its time past the end of the file, and solar cycle 27 is
> not predictable.

So V2 does not forecast. `ClimatologyDensity` holds the indices constant at a
stated activity level and reports a profile; the caller runs the levels
separately and quotes a band. Activity levels are the p05 / p50 / p95
percentiles of the **observed** F10.7 and Ap record in `SW-All.csv` (1957
onward, ~6 solar cycles), read at runtime.

The profile is the mean over 12 months and 8 longitudes at fixed indices,
which removes the seasonal and local-time cycles.

**Averaging density is exact here, and that is not a coincidence.** `da/dt` is
linear in `rho` (`PHYSICS.md` §3.2), so the mean density gives the mean secular
rate exactly, provided `a` does not move much within an averaging period. Over
a year at 700 km it does not.

**That argument does not extend to the solar cycle**, which is why the cycle is
swept rather than averaged. Decay *time* is strongly non-linear in the duration
spent at high density. A single run at mean F10.7 would be a precise answer to
the wrong question.

Altitude interpolation is log-linear on a 2 km grid to 2000 km — the upper edge
of the LEO region as 47 CFR § 25.283(e) defines it. Maximum interpolation error
against fresh pymsis calls is 0.125%, inside the `PHYSICS.md` §6.3 requirement
of 1%.

---

## V2.3 Disposal manoeuvre

Standard two-body results, written out rather than imported, for the same
reason the integrator is: "did you write it?" should keep its answer.

**Vis-viva.** Speed at radius `r` on an orbit of semi-major axis `a`:

```
v = sqrt(MU * (2/r - 1/a))
```

This is the one place in the codebase where `a` is a genuine ellipse
semi-major axis rather than also being the orbital radius; everywhere else the
circular assumption makes them equal (`PHYSICS.md` §2).

**Hohmann transfer between circular orbits** `r1 → r2`, via a transfer ellipse
with `a_t = (r1 + r2)/2`:

```
dv1 = |sqrt(MU*(2/r1 - 1/a_t)) - sqrt(MU/r1)|
dv2 = |sqrt(MU/r2) - sqrt(MU*(2/r2 - 1/a_t))|
```

**Single-burn perigee lowering** from a circular orbit at `r_c` to perigee
`r_p`, with `a_t = (r_c + r_p)/2`:

```
dv = sqrt(MU/r_c) - sqrt(MU*(2/r_c - 1/a_t))
```

**Tsiolkovsky**, solved for propellant mass, with `m0` the mass at the start of
the burn (end of mission, not launch):

```
m_prop = m0 * (1 - exp(-dv / (Isp * g0)))
```

### Verification

Hand-calculated for 705 km → 400 km circular, `MU = 3.986004418e14`,
`R_E = 6378.137 km`:

| Quantity | Value |
|---|---|
| `r1`, `r2` | 7083.137 km, 6778.137 km |
| `a_t` | 6930.637 km |
| `v1 = sqrt(MU/r1)` | 7501.6374 m/s |
| `vt1 = sqrt(MU*(2/r1 - 1/a_t))` | 7418.6462 m/s |
| **`dv1`** | **82.9911 m/s** |
| `v2 = sqrt(MU/r2)` | 7668.5582 m/s |
| `vt2 = sqrt(MU*(2/r2 - 1/a_t))` | 7752.4676 m/s |
| **`dv2`** | **83.9094 m/s** |
| **total** | **166.9006 m/s** |

Independently cross-checked against the first-order expansion for small
`|Δr|/r`, `dv ≈ (v/2)·|Δr|/r = 161.5 m/s` — 3.2% low, which is the expected
sign and magnitude at `|Δr|/r ≈ 0.043`.

---

## V2.4 Known limitation — disposal is priced on circular orbits

**This is the largest modelling compromise in V2 and it is not hidden in a
docstring.**

`PHYSICS.md` §1 excludes orbital eccentricity. Real post-mission disposal is
almost never circular: an operator fires once near apogee, drops perigee into
the thermosphere, and lets drag at perigee do the rest over many elliptical
revolutions. That is both cheaper and faster than the circular equivalent.

The compliance verdict therefore uses a **two-burn transfer to a lower circular
orbit** — the only disposal orbit this engine can actually propagate, so the
decay time behind the verdict is computed rather than asserted. The
single-burn perigee-lowering delta-v is also reported, because it is what an
operator would really fly and the delta-v itself is exact, but **no decay time
is attached to it** and it never enters a verdict.

Consequence, stated plainly:

> The verdict is a **conservative bound**. A design this module calls compliant
> is compliant. A design it calls non-compliant may still be recoverable with
> an elliptical profile the engine cannot yet represent.

`ComplianceResult.margin_note` carries this caveat with every result that
involves a burn, so it travels into the JSON export and the UI rather than
living only here.

Lifting this limitation means adding eccentricity to the state vector, which
is a change to `PHYSICS.md` §1 and §3, not a V2 addition.

### Eccentricity is the highest-value remaining physics extension

Of everything still missing from this model, adding `e` to the state buys the
most, and it is worth stating why rather than leaving it as a line in a
limitations list.

**It converts the disposal calculation from a bound into an answer.** Every
other limitation in `PHYSICS.md` §10 makes the existing numbers less precise.
This one makes a whole manoeuvre — the single-burn perigee drop that is what
operators actually fly — impossible to evaluate at all. The delta-v for it is
already computed and already correct (`perigee_lowering_dv`); what is missing
is only the ability to propagate the resulting ellipse and attach a decay time
to it. That is the gap between "this design needs at most 167 m/s" and "this
design needs 140 m/s and reenters in 3.2 years".

**The secular formulation should extend to it rather than being replaced.**
`PHYSICS.md` §3.2 is already an orbit-averaged rate, not an instantaneous one,
and the standard King-Hele treatment of drag on an eccentric orbit is of the
same kind: average the drag over one revolution and get coupled secular rates
in `a` and `e`, with the atmosphere sampled around the orbit rather than at a
single altitude. The state becomes `[a, e, m]`, `rk4_step` is untouched, the
adaptive controller extends by bounding fractional change in both `a` and `e`,
and the circular case must fall out exactly at `e = 0` — which is a sharp
regression test, since the entire V1 validation has to survive it unchanged.

The real work is in the atmosphere, not the integrator: a per-revolution
density average means sampling `rho` at several true anomalies per step
instead of once, so `DensityGrid` and `ClimatologyDensity` both need a path
that returns density along an orbit rather than at a point. That also happens
to be the natural place to lift the single-point latitude sampling of
`PHYSICS.md` §10.2, which is the other limitation that would become tractable
at the same time.

Not attempted in V2. Recorded here so the next scope has the argument already
made and knows what the regression test is.

---

## V2.5 Regulatory thresholds are data, not physics

`data/disposal_rules.json` holds each rule with its citation, quoted text,
jurisdiction and verification date. Two corrections to the summary in
`V2_BRIEF.md` §3 came out of checking it against 47 CFR § 25.283(e):

1. **The five-year clock starts at end of mission, not at launch.** A satellite
   with a 7-year design life at 705 km has 12 years from launch, not 5.
2. **The FCC rule is scoped to disposal by uncontrolled atmospheric re-entry
   from below 2000 km.** Controlled de-orbit and re-orbit to a disposal orbit
   are different paths with different requirements.

A third correction concerns the altitude table in `V2_BRIEF.md` §3 itself. Its
`da/dt` column reproduces well against this engine, but its **time-to-reentry**
column is a linear extrapolation and overstates by roughly 6–9×, because decay
accelerates sharply as the satellite falls. Computed properly at
`Bc = m/(Cd·A) = 100 kg/m²` and mean solar activity:

| Altitude | `V2_BRIEF` §3 estimate | Computed |
|---|---|---|
| 300 km | ~1 year | 0.17 yr |
| 400 km | ~13 years | 1.68 yr |
| 500 km | ~100+ years | 11.5 yr |
| 705 km | millennia | >200 yr |

This changes the answer to the question the brief poses. The brief's "below
~550 km, does natural decay comply?" turns out to describe the **25-year**
threshold (543 km at mean solar activity), not the FCC five-year one, which
lands near **454 km**:

| Rule | low solar | mean | high |
|---|---|---|---|
| FCC 5-year | 403 km | 454 km | 571 km |
| IADC 25-year | 475 km | 543 km | 701 km |

Pinned in `tests/test_disposal.py::test_the_briefs_550km_figure_is_the_25_year_threshold`.

---

## V2.6 Projected area from geometry (Phase 8)

At 200–700 km the flow is free-molecular: the mean free path is kilometres, so
molecules strike the satellite and are re-emitted without colliding with each
other on the way in. There is no boundary layer and no wake pressure recovery.
Drag is therefore proportional to the area geometrically projected onto the
plane perpendicular to the velocity vector, with the accommodation physics
bundled into `Cd` (`PHYSICS.md` §2).

**Convex body, projected onto the plane perpendicular to `v`:**

```
A_proj = (1/2) * sum_i A_i |n_i . v|
```

For a box with edges `L`, `W`, `T` along the body axes the six faces pair up:

```
A_proj = |vx|*(W*T) + |vy|*(L*T) + |vz|*(L*W)
```

giving the closed-form extremes used to validate the numerical direction
search:

```
min = min(W*T, L*T, L*W)                       (face-on to the smallest face)
max = sqrt((W*T)^2 + (L*T)^2 + (L*W)^2)        (corner-on, along the diagonal)
```

**Flat plate** of area `A` at angle `theta` between its normal and `v`:
`A_proj = A*cos(theta)`, the Phase 8 gate's analytic check. Plates carry a
1 cm thickness rather than zero, deliberately — a zero-thickness plate presents
exactly zero area edge-on, which would make a feathered array look drag-free.

### Two solvers, and why the difference is worth keeping

`projected_area` sums each part's analytic projection. `projected_area_union`
projects every vertex, takes each part's convex hull in the projection plane,
and rasterises the union. The first is fast and is an **upper bound**: two
parts overlapping in projection get counted twice. The second is correct under
occlusion. Their difference is the self-shadowing fraction — a quantity to
report, not an error to remove.

Same pattern as the `montecarlo.py` / `critical.py` agreement the project
already relies on: two paths through entirely different code reaching the same
number is the cheapest real correctness check available.

### Findings

**1 — The published 4.48 m² maximum is a face-on chassis, not a broadside.**
Chassis face-on is 2.8×1.3 = 3.64 m² or 3.0×1.5 = 4.50 m² for the two sourced
dimension pairs; the latter agrees with Baruah's published maximum to 0.45%. A
genuine broadside including the solar array is ~15 m², more than three times
larger. It is not even the chassis's own geometric maximum — corner-on gives
4.65 m², 3.3% higher. So 4.48 m² is best read as *chassis face-on, array
feathered*, which is consistent with the knife-edge attitude the February 2022
fleet was actually commanded into.

**2 — A fourth independent line on the effective drag parameter.** The README
documents three diagnostics converging on `rho*Cd*A`. This is a fourth, and
the first that never touches a decay curve:

| Route | `Cd` | `A` | `Cd*A` |
|---|---|---|---|
| Baruah et al. | 1.0 (stated simplification) | 1.00 m² | 1.00 m² |
| Geometry + free-molecular convention | 2.2 | 0.405 m² | 0.89 m² |

Two different splits of a product that is not separately observable, landing
11% apart. The geometry-derived knife-edge range (0.27–0.61 m²) also
reproduces the *secondary-source* range recorded in `satellite_specs.json`
(0.3–0.7 m²) rather than Baruah's 1.00 m² — which is a stated lower bound on a
swept range, not a measurement.

This is what V2_BRIEF.md §2 meant by pinning `A` from geometry to attack the
degeneracy rather than document it.

### Limitations specific to this solver

- **Attitude is an input, never an output.** The solver says what area an
  attitude presents; it never predicts what attitude a satellite settles into.
  Aerodynamic torque and attitude dynamics remain excluded (`PHYSICS.md` §1).
- **`NOMINAL` attitude is a conservative convention, not a measurement.** It is
  taken as the maximum-area orientation. Real operating attitude is set by
  pointing requirements. Sweep it if an answer depends on it.
- **The array chord is the weakest input.** Span is derived from a published
  ~11 m end-to-end figure; chord is assumed equal to the chassis width and is
  not independently sourced. It is tagged `estimated` and swept.
- **Rigid plates.** No panel flexure.
