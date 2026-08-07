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

