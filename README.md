# CastOrbit

A physics-grounded simulator of Starlink satellite survival in very low Earth
orbit during geomagnetic storms, validated against the published February 2022
Starlink loss before being pointed at anything else.

**The one-line result:** the model's drag term runs 15–20% low against reality
at 210 km, confirmed three independent ways — and once that single correction
is applied, the model's picture of the February 2022 loss is consistent with
what happened.

This is not a satellite tracker or a 3D visualisation of public TLE data. Every
number below was computed from the equations in [PHYSICS.md](PHYSICS.md) and
checked against a peer-reviewed result: Baruah, Y., et al. (2024). *The Loss of
Starlink Satellites in February 2022: How Moderate Geomagnetic Storms Can
Adversely Affect Assets in Low-Earth Orbit.* Space Weather, 22, e2023SW003716
([DOI: 10.1029/2023SW003716](https://doi.org/10.1029/2023SW003716)).

---

## The central finding

The equation this simulator integrates ([PHYSICS.md §3.2](PHYSICS.md#3-core-equation))
is

```
da/dt = 2*(F/m)*a**1.5/sqrt(MU) - rho*(Cd*A/m)*sqrt(MU*a)
```

Look at the drag term: `rho`, `Cd` and `A` never appear separately, only as the
product `rho * Cd * A`. Nothing in a decay curve alone can tell you whether
that product is low because the density model is low, because the drag
coefficient is off, or because the ram area is smaller than assumed — only the
product is observable from an altitude history. That degeneracy turns out to
be the whole story here.

**Three independent diagnostics find the same ~15–20% deficit in that
product**, each holding a different piece fixed:

**1 — Baruah's 4.48 m² bounding case.** Reproducing the paper's own convention
(Cd = 1.0, mass 227 kg, thrusters off) at the maximum ram area, the model
reenters at 45.81 h after launch against the published 38.75 h — 18.2% late.
Bisecting on a uniform multiplier applied to NRLMSIS density (holding Cd and A
at the paper's values) finds the model needs **×1.181** to hit the paper's
reentry time exactly.

**2 — Baruah's 1.00 m² bounding case.** Same convention, minimum ram area: the
model reaches 204.49 km at the reference time against the paper's 203.24 km,
decaying 5.51 km against a published 6.76 km — 18.5% too little. The same
bisection, run independently on this completely different trajectory shape
(this satellite never gets close to reentering), needs **×1.204**.

Two configurations 4.48× apart in drag area, one reentering and one barely
decaying, agree on the correction to within **2.0%**. That is the signature of
a single systematic offset in `rho·Cd·A`, not a bug specific to one
configuration — see `density_scale_diagnostic` in
[`sim/validate.py`](sim/validate.py).

**3 — The fleet loss count, from a completely different angle.** The first two
diagnostics adjust density while holding Cd and A at Baruah's own values. The
fleet reproduction does the opposite: it holds density fixed (scale 1.0) and
asks what ram-area range, combined with this project's own baseline drag
coefficient (Cd = 2.2, not Baruah's 1.0), reproduces the observed 38 losses out
of 49 satellites kept in safe mode (F = 0) from deployment to 2022-02-08.

Naively combining Cd = 2.2 with Baruah's published area range (1.00–4.48 m²)
loses **all 49** satellites — unsurprising, since that range was calibrated
against Cd = 1.0, and 2.2× the drag coefficient with the same area is roughly
2.2× the drag. Solving for the area range that *does* reproduce 38 losses at
Cd = 2.2 gives a uniform range of **0.52–2.34 m²**, i.e. the published range
scaled by **×0.523**.

That scaled range, converted back to the effective drag parameter that
actually enters the equation:

```
Cd * A  =  2.2 * [0.52, 2.34]  =  [1.15, 5.15] m²
```

Compare that to Baruah's own effective drag range, `Cd * A = 1.0 * [1.00,
4.48] = [1.00, 4.48] m²`. The two ranges nearly coincide — and the small gap
between them (roughly 15% at both ends) is the same ~15–20% correction found
by diagnostics 1 and 2, arrived at from the fleet's *loss count* rather than
any single satellite's decay curve, using a different Cd convention entirely.
Three methods, three different pieces of evidence, one number.

**What this means:** presenting "Cd = 2.2 loses everyone, Cd = 1.0 roughly
matches" as a contradiction between two arbitrary choices would be the wrong
reading. It is the same 15–20% effective-drag correction showing up a third
time, because `Cd` and `A` were never separable from each other or from `rho`
in the first place. The project's own runs use Cd = 2.2
([`satellite_specs.json`](data/satellite_specs.json), standard free-molecular-flow
convention); Baruah's use Cd = 1.0, stated as a simplification. Both are
consistent with the same physical satellite once the ~18% deficit is accounted
for.

![Baruah et al. (2024) reproduction: two bounding cases](out/baruah_validation.png)

*Both Baruah bounding cases, Cd = 1.0, 227 kg, thrusters off from the
2022-02-03 18:13 UT epoch. Stars mark the published targets at the 2022-02-05
08:58 UT reference time. Solid vs dashed lines compare the real 3-hourly space
weather history against a daily-average simplification — the two are visually
indistinguishable, discussed in [Uncertainty exceeds signal](#uncertainty-exceeds-signal).*

**What NRLMSIS's own known bias predicts.** This isn't a free parameter fit —
NRLMSIS 2.1 carries a documented 15–30% density error during geomagnetic
storms ([PHYSICS.md §10.4](PHYSICS.md#10-known-limitations)). A measured
15–20% deficit sits at the low edge of that band, not outside it. The
[Swarm C secondary validation](#validation) below, at a different altitude
entirely, points the same direction but by a larger margin — evidence the
offset is real but not a single altitude-independent constant.

---

## Validation

Four tests are specified in [PHYSICS.md §8](PHYSICS.md#8-validation); a fifth,
Swarm C, is an optional secondary check the paper itself uses.

| Test | What it checks | Result |
|---|---|---|
| **1 — Energy conservation** | `rho=0, F=0` ⇒ `da/dt=0` exactly, over a simulated week | Relative change in `a`: **0.000e+00** (limit 1e-12) |
| **2 — Pure thrust spiral** | Closed-form `a(t) = a0/(1-(F/m)t√(a0/MU))²` vs numerical, 1 day | **4.07e-15** relative (limit 1e-4). Convergence order independently confirmed at coarser steps: **4.00–4.03** (see `test_2b`/`test_2c` in [`tests/test_validation.py`](tests/test_validation.py) — at 10 s the error is round-off floor, not truncation) |
| **3 — Critical density fixed point** | `da/dt=0` at the computed `h_crit`, held 1 hour | `\|da/dt\|` / thrust term: **2.1e-16** (limit 1e-9) |
| **4 — Baruah reproduction** | Two bounding cases vs published targets | **+18.2%** and **−18.5%** decay-timing error — both inside the paper's own 20% acceptance band, and both explained by the single finding above |
| **Swarm C (secondary, flagged weakest)** | 434 km, 468 kg, 0.7 m², Cd 1.0, 53.78 h window | CastOrbit: 18.31 m decay. Paper's model: 25.02 m (**×1.366** implied). Observed: 23.08 m (**×1.260** implied) |

**Why Swarm C is flagged as the weakest evidence, not dropped.** Its implied
correction (×1.366) is ~14% larger than the Starlink pair's tight 2.0%
agreement, and it depends on an inclination (87.4°) that is **not in this
repo's data files** — I used the published Swarm C orbital inclination by the
same density-sampling convention applied elsewhere, but the sensitivity is
real: sampling latitude alone moves the implied correction from ×1.218 (0°) to
×1.449 (45°), a wider spread than the gap being explained. Swarm C confirms the
*direction* of the offset at a second altitude; it should not be read as
independently pinning its *size*. See `swarm_c_validation` in
[`sim/validate.py`](sim/validate.py).

Four validation tests, all passing, live in
[`tests/test_validation.py`](tests/test_validation.py),
[`tests/test_baruah.py`](tests/test_baruah.py),
[`tests/test_atmosphere.py`](tests/test_atmosphere.py) and
[`tests/test_critical.py`](tests/test_critical.py) — 39 tests total, run with
`pytest`.

---

