"""Display-only ground-track geometry for the Phase 5 globe view.

Not part of the validated physics model. PHYSICS.md §1 explicitly excludes
latitude/longitude variation and J2 effects ("negligible... over days"), and
this project is deliberately not a "3D visualisation of public TLE data"
(ARCHITECTURE.md §1). This module exists only because the globe view needs
SOME position to plot, and the Phase 5 hard rule is that nothing is computed
in the browser -- every number the frontend shows has to come from a Python
export, so a position has to be computed somewhere, and it has to be here.

What's real: the along-track angular rate at every sample comes from the
satellite's own already-integrated semi-major axis a(t) via mean motion, so a
satellite that decays visibly sinks toward Earth on the globe exactly as it
does on the altitude chart -- this is Kepler's third law applied to a number
the simulator already produced, not a new physical assumption.

What's illustrative, not measured: each satellite's initial orbital-plane
phase is spread evenly by id (plus a nudge from its own deploy time) purely so
49 satellites launched together are visually distinguishable. Real per-
satellite RAAN is not modelled or tracked anywhere in this project. Inclination
is held at the fleet's single published value and does not precess (no J2),
consistent with PHYSICS.md's own exclusion list.
"""

from __future__ import annotations

import numpy as np

from .constants import MU, R_E

# data/event_feb2022.json: orbit.inclination_deg. One value for the whole
# fleet -- this project has never modelled per-satellite inclination.
INCLINATION_DEG = 53.22

# Mean sidereal rotation rate, degrees per second.
EARTH_ROTATION_DEG_PER_S = 360.9856 / 86400.0

# Rough LEO period at ~210 km, used only to scale each satellite's real
# deploy-time offset into a modest additional phase spread (see module note).
_NOMINAL_PERIOD_S = 5400.0


def ground_track(
    t_s: np.ndarray,
    h_km: np.ndarray,
    satellite_id: int,
    n_satellites: int,
    deploy_time_s: float,
) -> tuple[list[float], list[float]]:
    """Sub-satellite (lat_deg, lon_deg) at each sample in `t_s`.

    `t_s` and `h_km` are this satellite's own downsampled trajectory -- the
    same arrays already exported for the altitude chart. No new propagation
    happens here beyond converting an already-known orbital radius into an
    angular rate; see the module docstring for what's real versus display-only.
    """
    t = np.asarray(t_s, dtype=float)
    a_m = R_E + np.asarray(h_km, dtype=float) * 1e3
    n_rad_s = np.sqrt(MU / a_m**3)  # instantaneous mean motion, rad/s

    dt = np.diff(t, prepend=t[0])
    phase0 = 2 * np.pi * satellite_id / max(n_satellites, 1)
    phase0 += 2 * np.pi * (deploy_time_s / _NOMINAL_PERIOD_S)
    u = phase0 + np.cumsum(n_rad_s * dt)  # argument of latitude

    i = np.radians(INCLINATION_DEG)
    lat = np.arcsin(np.sin(i) * np.sin(u))
    lon_inertial = np.arctan2(np.cos(i) * np.sin(u), np.cos(u))
    lon = lon_inertial - np.radians(EARTH_ROTATION_DEG_PER_S) * t

    lat_deg = np.degrees(lat)
    lon_deg = (np.degrees(lon) + 180.0) % 360.0 - 180.0  # wrap to [-180, 180]
    return (
        [round(float(v), 3) for v in lat_deg],
        [round(float(v), 3) for v in lon_deg],
    )


def cause_of_loss(outcome: str, outcome_time_iso: str | None, cd_times_area_m2: float) -> str:
    """Deterministic, human-readable summary of an already-known outcome.

    Every fact used here (outcome, outcome_time, Cd*A) is already in the
    export; this only phrases them for the globe's click panel, so it is not
    "generating" data, just formatting it.
    """
    if outcome == "REENTERED":
        return (
            f"Reentered (fell below 100 km) at {outcome_time_iso}. "
            f"Cd·A = {cd_times_area_m2:.2f} m² in safe mode (F = 0) "
            f"put it below the recovery boundary for the rest of the window."
        )
    if outcome == "INDETERMINATE":
        return (
            f"Still above 100 km when the observation window ended. "
            f"Cd·A = {cd_times_area_m2:.2f} m² was low enough, or the "
            f"window short enough, that drag had not yet won."
        )
    if outcome == "REACHED_SHELL":
        return "Reached the operational shell altitude."
    return "Propellant exhausted before reaching the shell."
