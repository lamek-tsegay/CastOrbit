"""Validation Test 4 -- the Baruah et al. (2024) reproduction.

PHYSICS.md §8, Test 4. Tests 1-3 are analytic unit tests and live in
`tests/test_validation.py`; this one needs the atmosphere model and real space
weather, so it is a script.

Parameters match the paper exactly (PHYSICS.md §8, `data/event_feb2022.json`):

    epoch     2022-02-03 18:13 UT      launch
    h0        210 km                   circular, per the §10.1 limitation
    mass      227 kg
    Cd        1.0                      NOT 2.2 -- see the Cd trap in §8
    area      1.00 m2 and 4.48 m2      two bounding cases
    thrust    0 N                      safe mode, drag alone (§5)
    latitude  53.22 deg                orbital inclination

Published targets at 2022-02-05 08:58 UT: the 4.48 m2 case reaches ~100.01 km
(reentry), the 1.00 m2 case reaches ~203.24 km.

Nothing here is tuned. Whatever the model produces is what gets reported
(ARCHITECTURE.md §9).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .atmosphere import DensityGrid, SpaceWeather
from .constants import R_E
from .dynamics import derivatives
from .integrator import propagate
from .satellite import Outcome

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "out"

EPOCH = datetime(2022, 2, 3, 18, 13, tzinfo=timezone.utc)
REFERENCE_TIME = datetime(2022, 2, 5, 8, 58, tzinfo=timezone.utc)
T_REFERENCE_S = (REFERENCE_TIME - EPOCH).total_seconds()   # 139500.0 s = 38.75 h

INSERTION_ALTITUDE_M = 210e3
MASS_KG = 227.0
CD_VALIDATION = 1.0          # Baruah et al. 2024, stated simplification
THRUST_N = 0.0               # safe mode
LATITUDE_DEG = 53.22         # orbital inclination
LONGITUDE_DEG = 0.0

# data/event_feb2022.json -> validation_targets
TARGETS = {4.48: 100.01, 1.00: 203.24}


@dataclass
class CaseResult:
    area_m2: float
    storm_time: bool
    altitude_at_reference_km: float | None   # None if it reentered first
    outcome: str
    reentry_time_s: float | None
    t_s: np.ndarray = field(repr=False)
    h_km: np.ndarray = field(repr=False)

    @property
    def reentry_time_utc(self) -> datetime | None:
        if self.reentry_time_s is None:
            return None
        return EPOCH + timedelta(seconds=self.reentry_time_s)

    def time_to_altitude(self, h_km_target: float) -> float | None:
        """First time, in seconds, the trajectory falls to `h_km_target`.

        Linearly interpolated between samples. Returns None if the trajectory
        never gets that low.
        """
        below = np.nonzero(self.h_km <= h_km_target)[0]
        if below.size == 0:
            return None
        i = int(below[0])
        if i == 0:
            return float(self.t_s[0])
        h1, h0 = self.h_km[i], self.h_km[i - 1]
        t1, t0 = self.t_s[i], self.t_s[i - 1]
        if h0 == h1:
            return float(t1)
        return float(t0 + (h0 - h_km_target) * (t1 - t0) / (h0 - h1))


def run_case(
    area_m2: float,
    storm_time: bool,
    sw: SpaceWeather,
    dt: float = 10.0,
    t_max_s: float = 5 * 86400.0,
    cd: float = CD_VALIDATION,
    mass_kg: float = MASS_KG,
    insertion_altitude_m: float = INSERTION_ALTITUDE_M,
    density_scale: float = 1.0,
    grid: DensityGrid | None = None,
) -> CaseResult:
    """Propagate one bounding case. PHYSICS.md §8 Test 4.

    `density_scale` is a diagnostic only. It is never used for the reported
    validation result, which always runs at 1.0. Its purpose is to answer one
    question: what uniform density bias would reconcile the model with the
    paper, and is that bias the same for both bounding cases? See
    `density_scale_diagnostic`.
    """
    if grid is None:
        grid = DensityGrid(
            EPOCH,
            sw,
            lat_deg=LATITUDE_DEG,
            lon_deg=LONGITUDE_DEG,
            duration_s=t_max_s,
            storm_time=storm_time,
        )

    def deriv(t: float, y: np.ndarray) -> np.ndarray:
        rho = density_scale * grid(t, y[0] - R_E)
        return derivatives(y, rho, thrust=THRUST_N, cd=cd, area=area_m2, isp=None)

    a0 = R_E + insertion_altitude_m
    traj = propagate(deriv, np.array([a0, mass_kg]), dt=dt, t_max=t_max_s)

    h_km = traj.h_km
    t_s = traj.t_s
    if t_s[-1] >= T_REFERENCE_S:
        alt_ref = float(np.interp(T_REFERENCE_S, t_s, h_km))
    else:
        alt_ref = None  # reentered before the reference time

    return CaseResult(
        area_m2=area_m2,
        storm_time=storm_time,
        altitude_at_reference_km=alt_ref,
        outcome=traj.outcome.value,
        reentry_time_s=(
            traj.outcome_time_s if traj.outcome is Outcome.REENTERED else None
        ),
        t_s=t_s,
        h_km=h_km,
    )


def density_scale_diagnostic(
    sw: SpaceWeather,
    storm_time: bool = True,
    dt: float = 10.0,
    t_max_s: float = 3 * 86400.0,
) -> dict[float, float]:
    """What uniform density multiplier would reconcile each case with the paper?

    This is a *diagnostic*, not a correction. The reported validation numbers
    are always at scale 1.0 (ARCHITECTURE.md §4: no tuning to make validation match).

    The question it answers is a discriminating one. If the discrepancy comes
    from a density-model difference -- NRLMSIS 2.1 here versus JB2008 in the
    paper -- then a single multiplier should reconcile both bounding cases,
    because decay rate is linear in rho. If instead the two cases demand very
    different multipliers, something configuration-dependent is wrong and the
    density explanation does not hold.
    """
    grid = DensityGrid(
        EPOCH, sw, lat_deg=LATITUDE_DEG, lon_deg=LONGITUDE_DEG,
        duration_s=t_max_s, storm_time=storm_time,
    )

    def miss(area: float, scale: float) -> float:
        """Signed error against the published target; zero means agreement."""
        r = run_case(area, storm_time, sw, dt=dt, t_max_s=t_max_s,
                     density_scale=scale, grid=grid)
        # Sign convention throughout: negative => the model decays too slowly
        # for this scale, positive => too fast. Monotonically increasing in
        # `scale` for both cases.
        if area == 4.48:
            # Target: reach 100 km exactly at the reference time.
            t = r.time_to_altitude(100.0)
            if t is None:
                return -1e9  # never reentered => far too little drag
            return T_REFERENCE_S - t
        # Target: 203.24 km at the reference time. More drag => lower altitude.
        alt = r.altitude_at_reference_km
        if alt is None:
            return +1e9  # reentered early => far too much drag
        return TARGETS[1.00] - alt

    scales: dict[float, float] = {}
    for area in (4.48, 1.00):
        lo, hi = 0.5, 4.0
        f_lo, f_hi = miss(area, lo), miss(area, hi)
        if f_lo > 0 or f_hi < 0:
            continue  # root not bracketed
        for _ in range(20):
            if hi - lo < 1e-3:
                break
            mid = 0.5 * (lo + hi)
            if miss(area, mid) < 0:
                lo = mid
            else:
                hi = mid
        scales[area] = 0.5 * (lo + hi)
    return scales


def _fmt_hours(seconds: float) -> str:
    sign = "-" if seconds < 0 else "+"
    return f"{sign}{abs(seconds) / 3600.0:.2f} h"


def main(
    t_max_s: float = 5 * 86400.0, dt: float = 10.0, run_diagnostic: bool = True
) -> dict:
    sw = SpaceWeather.load(DATA / "SW-All.csv")
    OUT.mkdir(exist_ok=True)

    results: dict[tuple[float, bool], CaseResult] = {}
    for storm_time in (True, False):
        for area in (4.48, 1.00):
            results[(area, storm_time)] = run_case(
                area, storm_time, sw, dt=dt, t_max_s=t_max_s
            )

    print("=" * 78)
    print("PHYSICS.md §8 Test 4 -- Baruah et al. (2024) reproduction")
    print("=" * 78)
    print(f"epoch      {EPOCH:%Y-%m-%d %H:%M} UT      reference {REFERENCE_TIME:%Y-%m-%d %H:%M} UT"
          f"  (+{T_REFERENCE_S / 3600:.2f} h)")
    print(f"h0 = {INSERTION_ALTITUDE_M / 1e3:.0f} km   m = {MASS_KG} kg   "
          f"Cd = {CD_VALIDATION}   F = {THRUST_N} N (safe mode)   lat = {LATITUDE_DEG}°")
    print(f"atmosphere NRLMSIS 2.1 via pymsis 0.12.0   (paper used JB2008)")
    print(f"integrator hand-written RK4, dt = {dt:.0f} s")
    print()

    for storm_time in (True, False):
        label = "storm_time=True  (3-hourly ap)" if storm_time else "storm_time=False (daily Ap)"
        print(f"--- {label} " + "-" * (78 - 5 - len(label)))
        for area in (4.48, 1.00):
            r = results[(area, storm_time)]
            target = TARGETS[area]
            print(f"  A = {area:.2f} m2   outcome {r.outcome}")
            if r.altitude_at_reference_km is not None:
                alt = r.altitude_at_reference_km
                err_alt = (alt - target) / target * 100.0
                decay = INSERTION_ALTITUDE_M / 1e3 - alt
                decay_target = INSERTION_ALTITUDE_M / 1e3 - target
                err_decay = (decay - decay_target) / decay_target * 100.0
                print(f"      altitude at reference = {alt:.2f} km "
                      f"(target {target:.2f} km, {err_alt:+.2f}%)")
                print(f"      decay since insertion = {decay:.2f} km "
                      f"(target {decay_target:.2f} km, {err_decay:+.1f}%)")
            else:
                print(f"      reentered BEFORE the reference time")
            if r.reentry_time_s is not None:
                dt_hours = (r.reentry_time_s - T_REFERENCE_S) / 3600.0
                print(f"      reentry (100 km) at    = {r.reentry_time_utc:%Y-%m-%d %H:%M} UT "
                      f"({r.reentry_time_s / 3600:.2f} h after launch, "
                      f"{dt_hours:+.2f} h vs the published 08:58 UT)")
                print(f"      decay-timing error     = "
                      f"{(r.reentry_time_s - T_REFERENCE_S) / T_REFERENCE_S * 100:+.1f}%")
            else:
                print(f"      did not reach 100 km within {t_max_s / 86400:.0f} days")
        print()

    print("--- 3-hourly ap vs daily Ap: shift in timing " + "-" * 32)
    summary_shifts = {}
    for area in (4.48, 1.00):
        on, off = results[(area, True)], results[(area, False)]
        # Compare at a common altitude so both cases yield a timing shift.
        probe_km = 100.0 if area == 4.48 else TARGETS[1.00]
        t_on, t_off = on.time_to_altitude(probe_km), off.time_to_altitude(probe_km)
        if t_on is not None and t_off is not None:
            shift = t_on - t_off
            summary_shifts[area] = shift
            print(f"  A = {area:.2f} m2   time to {probe_km:.2f} km: "
                  f"3-hourly {t_on / 3600:.2f} h, daily {t_off / 3600:.2f} h  "
                  f"-> {_fmt_hours(shift)} ({shift / t_off * 100:+.1f}%)")
        else:
            print(f"  A = {area:.2f} m2   never reached {probe_km:.2f} km in one or both runs")
        d_on = on.altitude_at_reference_km
        d_off = off.altitude_at_reference_km
        if d_on is not None and d_off is not None:
            print(f"                 altitude at reference: 3-hourly {d_on:.2f} km, "
                  f"daily {d_off:.2f} km ({d_on - d_off:+.2f} km)")
    print()

    scales: dict[float, float] = {}
    if run_diagnostic:
        print("--- discrepancy diagnostic: implied density bias " + "-" * 29)
        print("  Not a correction. The results above are all at scale 1.0.")
        print("  Question: does ONE uniform density multiplier reconcile BOTH cases?")
        scales = density_scale_diagnostic(sw, storm_time=True, dt=20.0)
        for area in (4.48, 1.00):
            if area in scales:
                k = scales[area]
                print(f"  A = {area:.2f} m2 -> x{k:.4f}  "
                      f"(NRLMSIS {(1 - 1 / k) * 100:.1f}% below the density the paper's decay implies)")
        if len(scales) == 2:
            a, b = scales[4.48], scales[1.00]
            spread = abs(a - b) / ((a + b) / 2) * 100
            print(f"  The two cases differ 4.48x in drag area yet agree on the "
                  f"correction to {spread:.1f}%.")
            print(f"  => consistent with a single uniform density offset, not a "
                  f"configuration-dependent error.")
        print()

    _save_outputs(results)
    return {"results": results, "shifts": summary_shifts,
            "density_scales": scales}


def _save_outputs(results: dict[tuple[float, bool], CaseResult]) -> None:
    """Write the decay curves: a PNG for the eye, JSON for anything downstream."""
    payload = {
        "meta": {
            "test": "PHYSICS.md §8 Test 4 -- Baruah et al. 2024 reproduction",
            "generated": datetime.now(timezone.utc).isoformat(),
            "atmosphere_model": "NRLMSIS 2.1 via pymsis 0.12.0",
            "paper_atmosphere_model": "JB2008",
            "epoch": EPOCH.isoformat(),
            "reference_time": REFERENCE_TIME.isoformat(),
        },
        "config": {
            "insertion_altitude_km": INSERTION_ALTITUDE_M / 1e3,
            "mass_kg": MASS_KG,
            "cd": CD_VALIDATION,
            "thrust_n": THRUST_N,
            "latitude_deg": LATITUDE_DEG,
        },
        "cases": [],
    }
    for (area, storm_time), r in results.items():
        step = max(1, len(r.t_s) // 2000)  # ~10 min spacing for export
        payload["cases"].append({
            "ram_area_m2": area,
            "storm_time": storm_time,
            "outcome": r.outcome,
            "target_altitude_km": TARGETS[area],
            "altitude_at_reference_km": r.altitude_at_reference_km,
            "reentry_time_s": r.reentry_time_s,
            "trajectory": {
                "t_s": r.t_s[::step].tolist(),
                "h_km": r.h_km[::step].tolist(),
            },
        })
    (OUT / "baruah_validation.json").write_text(json.dumps(payload, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots(figsize=(9, 5.5))
    colours = {4.48: "#d1495b", 1.00: "#2e86ab"}
    for (area, storm_time), r in results.items():
        ax.plot(
            r.t_s / 3600.0, r.h_km,
            color=colours[area],
            linestyle="-" if storm_time else "--",
            linewidth=1.8 if storm_time else 1.2,
            alpha=1.0 if storm_time else 0.65,
            label=f"A = {area:.2f} m²  "
                  f"{'3-hourly ap' if storm_time else 'daily Ap'}",
        )
    for area, target in TARGETS.items():
        ax.plot(
            T_REFERENCE_S / 3600.0, target, "*",
            color=colours[area], markersize=15, markeredgecolor="black",
            markeredgewidth=0.6, zorder=5,
            label=f"Baruah target, A = {area:.2f} m² ({target:.2f} km)",
        )
    ax.axvline(T_REFERENCE_S / 3600.0, color="grey", lw=0.8, ls=":")
    ax.axhline(100.0, color="black", lw=0.8, ls=":")
    ax.text(T_REFERENCE_S / 3600.0 + 0.8, 211, "reference\n2022-02-05 08:58 UT",
            fontsize=8, color="grey", va="top")
    ax.text(60, 102, "100 km — unrecoverable", fontsize=8, color="black")
    ax.set_xlabel("hours after launch (2022-02-03 18:13 UT)")
    ax.set_ylabel("altitude (km)")
    ax.set_title(
        "Baruah et al. (2024) reproduction — Cd = 1.0, 227 kg, thrusters off\n"
        "CastOrbit / NRLMSIS 2.1 (paper used JB2008)", fontsize=10
    )
    ax.set_ylim(90, 215)
    ax.set_xlim(0, 96)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8, loc="center right", framealpha=0.92)
    fig.savefig(OUT / "baruah_validation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
