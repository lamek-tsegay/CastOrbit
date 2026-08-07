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

