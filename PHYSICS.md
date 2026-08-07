# PHYSICS.md — CastOrbit

Physics specification for the CastOrbit satellite decay simulator.

This document defines every equation to be implemented. Implement these
equations as written. Do not substitute alternative formulations without
noting the change here.

---

## 1. Scope

We model the altitude of a satellite in a near-circular low Earth orbit as a
function of time, under two competing effects:

- **Atmospheric drag**, which removes orbital energy and lowers the orbit
- **Continuous low thrust** from an electric propulsion system, which adds
  orbital energy and raises the orbit

The model answers one question: **given an insertion altitude, a satellite
configuration, and a space weather history, does the satellite climb to its
operational shell or decay and reenter?**

### Explicitly out of scope

The following are deliberately excluded. They are documented here so the
limitations section of the writeup is accurate.

| Excluded | Why |
|---|---|
| Orbital eccentricity | Circular approximation only; real insertion was elliptical |
| J2 and higher gravity harmonics | Negligible effect on energy budget over days |
| Attitude dynamics | Attitude enters only as a choice of ram area `A` |
| Thruster duty cycling | Thrust is either full-on or off |
| Solar radiation pressure | Orders of magnitude below drag at 210 km |
| Lift, or any non-drag aerodynamic force | Free-molecular flow, drag dominates |
| Latitude/longitude variation of density | Density sampled at a fixed representative latitude |
| Satellite-to-satellite interaction | Each satellite is independent |

---

## 2. Notation and constants

| Symbol | Meaning | Units |
|---|---|---|
| `a` | Orbital semi-major axis (radius, since circular) | m |
| `h` | Altitude above Earth's surface, `h = a - R_e` | m |
| `m` | Satellite mass | kg |
| `A` | Ram cross-sectional area (area facing velocity vector) | m² |
| `Cd` | Drag coefficient | dimensionless |
| `F` | Thruster force | N |
| `Isp` | Specific impulse | s |
| `rho` | Atmospheric mass density at altitude `h` | kg/m³ |
| `v` | Orbital speed | m/s |

Constants (use these exact values):

```python
MU    = 3.986004418e14   # Earth gravitational parameter, m^3/s^2
R_E   = 6378.137e3       # Earth equatorial radius, m
G0    = 9.80665          # standard gravity, m/s^2
```

---

## 3. Core equation

### 3.1 Derivation

For a circular orbit, specific orbital energy is

```
eps = -MU / (2*a)
```

and orbital speed is

```
v = sqrt(MU / a)
```

Differentiating the energy expression:

```
d(eps)/dt = (MU / (2 * a^2)) * da/dt
```

so

```
da/dt = (2 * a^2 / MU) * d(eps)/dt          ... (1)
```

Any tangential acceleration `a_t` changes specific energy at a rate
`d(eps)/dt = a_t * v`.

**Drag** produces a tangential acceleration opposing motion:

```
a_drag = -(1/2) * rho * v^2 * (Cd * A / m)
```

Substituting into (1) with `v = sqrt(MU/a)` and simplifying:

```
da/dt |_drag = -rho * (Cd * A / m) * sqrt(MU * a)
```

**Thrust** produces a tangential acceleration along motion:

```
a_thrust = F / m
```

Substituting into (1):

```
da/dt |_thrust = 2 * (F / m) * a^(3/2) / sqrt(MU)
```

### 3.2 The equation to implement

```
da/dt = 2 * (F/m) * a**1.5 / sqrt(MU)  -  rho(h, t) * (Cd*A/m) * sqrt(MU*a)
```

with mass depletion from propellant consumption:

```
dm/dt = -F / (Isp * G0)
```

When the thruster is off (see §5), set `F = 0` in both equations.

This is the entire simulator. Two coupled ODEs, state vector `[a, m]`.

### 3.3 Termination

Stop integration and record the outcome when any of these occur:

| Condition | Outcome label |
|---|---|
| `h < 100e3` (100 km) | `REENTERED` — unrecoverable |
| `h >= target_shell_altitude` | `REACHED_SHELL` — success |
| Propellant mass exhausted | continue with `F = 0`, outcome resolves later |
| Simulation time limit reached | `INDETERMINATE` — report separately |

The 100 km threshold is the "unrecoverable" altitude used by Baruah et al.
(2024). Use the same value so results are comparable.

---

## 4. Critical altitude

Setting `da/dt = 0` in §3.2 and solving for density:

```
2 * F * a**1.5 / sqrt(MU)  =  rho * (Cd*A) * sqrt(MU*a)
```

```
rho_crit = 2 * F * a / (Cd * A * MU)
```

**Interpretation.** If the local density is *below* `rho_crit`, thrust wins and
the satellite climbs. If *above*, drag wins and it sinks. Because density rises
steeply as altitude falls, sinking increases drag, which increases sinking —
the process runs away in both directions. This is the knife edge the whole
project is about.

To get the **critical altitude** `h_crit`, invert the atmosphere model: find the
altitude at which `rho(h) == rho_crit`. Use bisection over `h` in `[100e3, 600e3]`.
Density is monotonically decreasing in altitude, so bisection is safe.

### 4.1 Verified sanity values

These were computed with `pymsis` 0.12.0 at lat 53.22°, 2022-02-04,
F10.7 = 127, F10.7A = 110, `F = 0.071 N`, `a = R_E + 210 km`.
**Your implementation should reproduce these.**

Critical densities:

| Cd | A (m²) | rho_crit (kg/m³) |
|---|---|---|
| 1.0 | 1.00 | 2.347e-09 |
| 1.0 | 4.48 | 5.239e-10 |
| 2.2 | 1.00 | 1.067e-09 |
| 2.2 | 4.48 | 2.381e-10 |

Actual densities from NRLMSIS 2.1:

| Altitude | quiet (ap=5) | storm (ap=56, ap3=80) |
|---|---|---|
| 180 km | 4.118e-10 | 4.685e-10 |
| 200 km | 2.030e-10 | 2.360e-10 |
| 210 km | 1.468e-10 | 1.726e-10 |
| 230 km | 8.010e-11 | 9.632e-11 |
| 260 km | 3.495e-11 | 4.347e-11 |

**Read the implication.** For the worst-case configuration
(`Cd = 2.2`, `A = 4.48 m²`), `rho_crit = 2.381e-10`. During the storm, density
reaches that value at approximately **200 km**. The satellites were inserted at
**210 km**. The margin was about 10 km.

This falls directly out of the formula and is a headline result. Compute it,
don't hardcode it.

> **Implementation note added in Phase 1.** The storm column reproduces exactly
> with daily `Ap = 56` under pymsis's *default* options. The `ap3 = 80` in the
> heading has no effect on those numbers: NRLMSIS 2.1 reads only the daily Ap
> unless the `geomagnetic_activity` option is set to `-1`. See §6.2.

---

## 5. Safe mode — the physical crux

**This is the most important modeling decision in the project.**

During the February 2022 storm, SpaceX commanded the satellites into a
low-drag "knife-edge" survival attitude. In that configuration they were not
raising their orbits. Drag acted alone.

The naive model — thruster on, fighting drag continuously — will predict that
most satellites survive, because at 210 km the thrust term generally beats the
drag term. That prediction is wrong, and the reason it is wrong is
operational rather than physical.

Implement thruster state as an explicit, time-dependent input:

```python
class ThrusterMode(Enum):
    NOMINAL   = "nominal"     # F = rated thrust, A = nominal ram area
    SAFE_MODE = "safe_mode"   # F = 0,            A = knife-edge ram area
```

A run is defined partly by **when** the satellite enters safe mode and **whether
it ever exits**. Make this a first-class simulation parameter, not a constant.

The interesting question the simulator can then answer:
*given the storm, how late could safe mode have been exited and still recovered?*

---

## 6. Atmosphere model

Use `pymsis` (NRLMSIS 2.1). Confirmed working: `pip install pymsis`, version
0.12.0, prebuilt wheel, no compilation required.

### 6.1 Call signature

```python
from pymsis import msis
import numpy as np

output = msis.calculate(dates, lons, lats, alts, f107s, f107as, aps)
# output[..., 0] is total mass density in kg/m^3
```

**`alts` is in kilometres, not metres.** Convert.

### 6.2 The `aps` array

`aps` must be a sequence of 7-element sequences, one per date. The elements are:

| Index | Meaning | Source column in `SW-All.csv` |
|---|---|---|
| 0 | Daily Ap | `AP_AVG` |
| 1 | 3-hour ap, current | `AP1`–`AP8` by UT slot |
| 2 | 3-hour ap, 3 h before | previous slot |
| 3 | 3-hour ap, 6 h before | two slots back |
| 4 | 3-hour ap, 9 h before | three slots back |
| 5 | Mean of eight 3-h ap, 12–33 h before | slots 4–11 back |
| 6 | Mean of eight 3-h ap, 36–57 h before | slots 12–19 back |

**Verified `SW-All.csv` format** (CelesTrak, 25372 rows, 1957-10-01 onward):

```
DATE,BSRN,ND,KP1..KP8,KP_SUM,AP1..AP8,AP_AVG,CP,C9,ISN,
F10.7_OBS,F10.7_ADJ,F10.7_DATA_TYPE,
F10.7_OBS_CENTER81,F10.7_OBS_LAST81,F10.7_ADJ_CENTER81,F10.7_ADJ_LAST81
```

Column mapping for `msis.calculate`:

| pymsis argument | Column |
|---|---|
| `f107s` | `F10.7_OBS` — **of the previous day** (NRLMSIS convention) |
| `f107as` | `F10.7_OBS_CENTER81` — **centred**, not `LAST81` |
| `aps` | `AP_AVG` and `AP1`–`AP8`, assembled per the table above |

`AP1`–`AP8` map to UT slots 00–03, 03–06, 06–09, 09–12, 12–15, 15–18, 18–21, 21–00.

Three traps:

- **Use `CENTER81`, not `LAST81`.** `LAST81` is the trailing average used for
  real-time forecasting. For a historical reconstruction the centred average is
  correct, and it is available because the future data exists.
- **`f107s` takes the previous day's value.** This is the NRLMSIS convention and
  is easy to miss.
- **The file has CRLF line endings.** Strip them.

**A fourth trap, found in Phase 1:**

- **The 3-hourly `ap` history is ignored unless you ask for it.** NRLMSIS 2.1
  uses only `aps[0]` (the daily Ap) under pymsis's default options. Elements
  `aps[1:]` — everything the table above so carefully assembles — are inert
  until the model option `geomagnetic_activity` is set to `-1`
  (`msis.create_options(geomagnetic_activity=-1)`). At 210 km on 2022-02-04,
  switching it on changes density by ~1.2%, and giving it a realistic spiky
  history rather than a flat one changes it by a further ~3.4%. The §4.1
  reference table was computed with the default, so reproducing that table and
  running a storm with 3-hourly structure are two different call configurations.
  `sim/atmosphere.py` exposes this as `storm_time=`.

**Verified values for the validation window** (read directly from the file):

| Date | AP_AVG | AP1–AP8 | F10.7_OBS | F10.7_OBS_CENTER81 |
|---|---|---|---|---|
| 2022-02-03 | 26 | 7, 32, 48, 56, 27, 15, 15, 12 | 126.5 | 109.1 |
| 2022-02-04 | 32 | 27, 27, 22, 22, 27, 56, 48, 27 | 129.6 | 108.8 |
| 2022-02-05 | 11 | 12, 15, 7, 9, 15, 12, 12, 6 | 125.9 | 108.5 |

Note the storm structure: ap peaks at 56 during 09–12 UT on 3 Feb, **before** the
18:13 UT launch, then peaks at 56 again during 15–18 UT on 4 Feb. This matches
the paper's account — launch occurred ~6 h after CME1's passage, with a second
shock arriving ~6 h after launch. The satellites were inserted into an already
disturbed thermosphere and then hit again.

**Write a unit test for the loader** that reproduces the table above. An
off-by-one in the history window is silent and will corrupt every result.

### 6.3 Caching

Density is queried at every integration step. Querying `pymsis` per-step is
slow. Precompute a density grid over `(altitude, time)` at the start of each
run and interpolate. Suggested grid: 1 km altitude spacing over
100–600 km, 30 minute time spacing. Verify the interpolation error is under 1%
against direct calls before relying on it.

> **Phase 1 note.** A single uniform grid does *not* meet the 1% target: the
> NRLMSIS inputs are piecewise constant in time and step at every 3 h UT
> boundary, and linear interpolation across a step has a peak error equal to
> the step regardless of spacing (measured 5.6% near 550 km). The implemented
> cache builds one interpolator per 3 h block, so no cell spans a
> discontinuity. Measured worst case is then 0.20%. Interpolation is linear in
> `log(rho)` rather than `rho`, since density is near-exponential in altitude.

---

## 7. Numerical integration

Use **RK4** with fixed step. Implement it directly — do not call
`scipy.integrate.solve_ivp` for the main loop. The first question anyone asks
about this project is whether the integrator is yours.

Suggested step: 10 s. **Verify by halving**: run at 10 s and 5 s; if the final
altitude differs by more than 0.1%, reduce the step further.

State vector: `[a, m]`. Both must remain positive; assert this each step.

---

## 8. Validation

Four tests. All four must pass and their results must appear in the writeup,
including any that fail.

### Test 1 — Energy conservation (analytic)

Set `rho = 0` and `F = 0`. The orbit must be exactly stable: `da/dt = 0`.
After a simulated week, `a` must be unchanged to within floating point error
(relative change < 1e-12).

*Catches: sign errors, spurious terms, integrator bugs.*

### Test 2 — Pure thrust spiral (analytic)

Set `rho = 0`, `F` nonzero, and hold `m` constant. Then

```
da/dt = 2 * (F/m) * a**1.5 / sqrt(MU)
```

has a closed-form solution. Separating variables:

```
a(t) = a0 / (1 - (F/m) * t * sqrt(a0/MU))**2
```

Numerical and analytic results must agree to within 0.01% over a
one-day propagation.

*Catches: incorrect exponents, wrong constant factors.*

### Test 3 — Critical density (analytic)

Place a satellite at exactly `h_crit` computed from §4. Run for one hour with
constant density. `da/dt` must be zero to within 1e-9 relative.

*Catches: inconsistency between the decay equation and the critical-altitude solver.*

### Test 4 — Baruah et al. (2024) reproduction

The published reference case. **Parameters must match the paper exactly:**

| Parameter | Value | Note |
|---|---|---|
| Epoch | 2022-02-03 18:13 UT | Launch time |
| Initial altitude | 210 km | |
| Inclination | 53.22° | Used for density latitude |
| Mass | 227 kg | Paper's value |
| Cd | **1.0** | **Not 2.2 — see below** |
| Ram area | 1.00 m² and 4.48 m² | Two bounding cases |
| Thrust | 0 N | Safe mode |

Published targets:

| Case | Target |
|---|---|
| `A = 4.48 m²` | Reaches ~100 km by 08:58 UT, 5 Feb 2022 |
| `A = 1.00 m²` | Reaches ~203.24 km by the same time |

> **Cd trap.** The paper uses `Cd = 1.0`, stated explicitly as a simplification
> for comparative analysis. Running validation at `Cd = 2.2` produces roughly
> double the drag, the model reenters far too early, and hours get lost hunting
> a bug that does not exist. **Validate at Cd = 1.0.** Then rerun at 2.2 for the
> project's own results and report the difference — that comparison is itself a
> finding.

> **Model difference.** Baruah et al. use JB2008; this project uses NRLMSIS 2.1.
> Exact agreement is not expected. Report the discrepancy rather than tuning
> parameters to hide it. Agreement within roughly 20% on decay timing is a good
> result and should be stated as such.

---

## 9. Monte Carlo

A batch is 49 satellites. Each draws independently:

| Parameter | Distribution | Justification |
|---|---|---|
| Ram area `A` | Uniform, 1.00–4.48 m² | Published bounding range; attitude not individually known |
| Mass `m` | Normal, mean 227 kg, sigma 3% | Manufacturing and propellant load variation |
| Insertion altitude | Normal, mean 210 km, sigma 2 km | Injection dispersion |
| Insertion time offset | Uniform, ±30 min | Deployment sequence spread |
| Thrust `F` | Normal, rated, sigma 2% | Unit-to-unit performance |
| `Cd` | Uniform, 2.0–2.4 | Free-molecular flow uncertainty |
| Safe mode entry | Scenario parameter | Not random — swept deliberately |

Outcome taxonomy, reported as counts and percentages:

- `REACHED_SHELL` — climbed to operational altitude
- `REENTERED` — fell below 100 km
- `PROPELLANT_EXHAUSTED` — ran dry before reaching shell
- `INDETERMINATE` — still in flight at simulation end

### Required sweeps

1. **Insertion altitude sweep.** 190–320 km in 10 km steps. Plot survival
   fraction against insertion altitude, one curve for quiet conditions and one
   for the storm. The gap between those curves is the central engineering result.

2. **Ram area sensitivity.** Survival fraction against ram area. Shows how much
   the published area uncertainty actually matters — turning a data gap into a
   quantified result.

3. **Safe mode timing.** Survival fraction against safe-mode exit time.

---

## 10. Known limitations

Reproduce this list in the writeup. Stating these first is the difference
between a project that survives questioning and one that does not.

1. **Circular orbit assumption.** The real insertion was elliptical with ~210 km
   perigee. Drag is concentrated near perigee, so a circular model at perigee
   altitude overestimates the time-averaged drag, while a model at mean altitude
   underestimates it. Quantify the direction of this bias if time allows.

2. **Single-point density sampling.** Density is evaluated at one
   latitude/longitude rather than integrated around the orbit. Real storm-time
   density enhancement is strongly non-uniform, which the paper cites as a
   possible reason 11 satellites survived.

3. **Tangential thrust only.** Real orbit-raising uses steering laws that are not
   purely tangential.

4. **Empirical atmosphere model error.** NRLMSIS carries roughly 15–30% density
   error during geomagnetic storms. This is the dominant uncertainty in the
   entire model and exceeds every parameter spread in §9.

5. **Ram area is inferred, not known.** SpaceX does not publish effective drag
   area. The 1.00–4.48 m² range comes from Baruah et al. and is itself an
   estimate from assumed geometry.

6. **No attitude dynamics.** Attitude is a choice of `A`, not a simulated state.
   Real satellites tumbled.

---

## 11. Sources

- Baruah, Y., et al. (2024). *The Loss of Starlink Satellites in February 2022:
  How Moderate Geomagnetic Storms Can Adversely Affect Assets in Low-Earth
  Orbit.* Space Weather, 22, e2023SW003716.
  Primary source for event parameters and validation targets.
- Emmert, J. T., et al. NRLMSIS 2.0/2.1. Atmosphere model underlying `pymsis`.
- CelesTrak `SW-All.csv`. Space weather indices.
- Satellite hardware specifications: see `data/satellite_specs.json`, where each
  field carries a `source` and a `confidence` tag. Note that the v1.5 thruster
  figures (71 mN, 1666 s) are **derived** from SpaceX's published relative
  scaling claims against the v2 mini, not directly published. Label them as
  derived wherever they appear.
