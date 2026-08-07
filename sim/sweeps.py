"""The three required sweeps. PHYSICS.md §9.

  1. Insertion altitude, 190-320 km, quiet vs storm  -- the central engineering
     result, and the one the whole project is arguing about.
  2. Ram area sensitivity -- turns the disputed nominal-area data gap in
     `data/satellite_specs.json` into a quantified survival penalty.
  3. Safe-mode exit timing -- "how late could safe mode have been exited and
     still recovered?" (PHYSICS.md §5).

Every sweep is run at density scale 1.0 and 1.19 and plotted as a band. The
1.19 comes from the Phase 2 validation: the Baruah bounding cases each imply a
uniform density multiplier of 1.181 and 1.204 to reconcile NRLMSIS 2.1 with the
paper's JB2008 result (`sim/validate.py`). Carrying it as a band rather than
folding it into the model keeps the uncorrected result visible.

Matplotlib only; nothing here touches the frontend.
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .atmosphere import SpaceWeather
from .constants import R_E
from .critical import critical_altitude
from .montecarlo import (
    CD_RANGE,
    N_SATELLITES,
    RAM_AREA_RANGE,
    THRUST_RATED_N,
    build_grid,
    run_batch,
    sample_batch,
)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "out"

LAUNCH = datetime(2022, 2, 3, 18, 13, tzinfo=timezone.utc)
# The clock starts half an hour before nominal deployment so the +/-30 min
# deployment spread of PHYSICS.md §9 lies entirely inside the density grid.
EPOCH = LAUNCH - timedelta(seconds=1800)

T_MAX_S = 15 * 86400.0
DT_S = 60.0
SEED = 20220203

DENSITY_SCALES = (1.00, 1.19)
CONDITIONS = ("storm", "quiet")

# §9 asks for 190-320 km in 10 km steps. The survival curve turns out to be a
# step at ~185-200 km and flat at 1.0 above it, so the specified grid resolves
# almost none of the interesting structure. The coarse grid is run as specified
# and a finer one is added across the knee.
ALTITUDES_KM = sorted(set(
    list(np.arange(190.0, 320.0 + 1e-9, 10.0))
    + list(np.arange(176.0, 214.0 + 1e-9, 2.0))
))
# satellite_specs.json v1_5.ram_area_nominal_m2 is DISPUTED at 0.5-6.0 m2.
RAM_AREAS_M2 = [0.5, 1.0, 1.5, 2.0, 2.74, 3.5, 4.48, 5.0, 5.5, 6.0]
SAFE_MODE_EXIT_H = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0, 54.0,
                    60.0, 66.0, 72.0, 84.0, 96.0, 108.0, 120.0]


def _grids(sw: SpaceWeather) -> dict[str, object]:
    return {c: build_grid(sw, EPOCH, storm=(c == "storm"), duration_s=T_MAX_S)
            for c in CONDITIONS}


def critical_altitude_band(grid, cd: float, area: float,
                           thrust: float = THRUST_RATED_N,
                           window_s: float = 2 * 86400.0) -> tuple[float, float]:
    """Min and max critical altitude over the first `window_s` of the run.

    PHYSICS.md §4. Computed from the same density grid the Monte Carlo uses, so
    the vertical line drawn on the sweep plot and the simulated transition are
    derived independently from one another -- if they disagree, one of them is
    wrong.
    """
    h_crits = []
    for t in np.arange(0.0, window_s, 3600.0):
        try:
            h_crits.append(
                critical_altitude(lambda h: float(grid.lookup(t, np.array([h]))[0]),
                                  thrust, cd, area)
            )
        except ValueError:
            continue
    if not h_crits:
        return (float("nan"), float("nan"))
    return (min(h_crits) / 1e3, max(h_crits) / 1e3)


def sweep_insertion_altitude(grids, dt: float = DT_S) -> dict:
    """PHYSICS.md §9, sweep 1. Thrusters nominal from deployment."""
    results = {}
    for cond in CONDITIONS:
        for scale in DENSITY_SCALES:
            fracs = []
            for h0 in ALTITUDES_KM:
                rng = np.random.default_rng(SEED)
                r = run_batch(
                    sample_batch(rng, h0), grids[cond], dt=dt, t_max_s=T_MAX_S,
                    density_scale=scale,
                )
                fracs.append(r.survival_fraction)
            results[(cond, scale)] = np.array(fracs)
            print(f"  sweep 1 {cond:5s} x{scale:.2f} done")
    return results


def sweep_ram_area(grids, insertion_km: float = 210.0, dt: float = DT_S) -> dict:
    """PHYSICS.md §9, sweep 2. Fixed insertion altitude, ram area swept."""
    results = {}
    for cond in CONDITIONS:
        for scale in DENSITY_SCALES:
            fracs = []
            for area in RAM_AREAS_M2:
                rng = np.random.default_rng(SEED)
                r = run_batch(
                    sample_batch(rng, insertion_km, area_m2=area), grids[cond],
                    dt=dt, t_max_s=T_MAX_S, density_scale=scale,
                )
                fracs.append(r.survival_fraction)
            results[(cond, scale)] = np.array(fracs)
            print(f"  sweep 2 {cond:5s} x{scale:.2f} done")
    return results


def sweep_safe_mode_timing(grids, insertion_km: float = 210.0,
                           dt: float = DT_S) -> dict:
    """PHYSICS.md §9, sweep 3 and the §5 question.

    Satellites are in safe mode (F = 0, knife-edge area) from deployment and
    exit at the swept time. Exit at 0 h is the same as never entering.
    """
    results = {}
    for cond in CONDITIONS:
        for scale in DENSITY_SCALES:
            fracs = []
            for exit_h in SAFE_MODE_EXIT_H:
                rng = np.random.default_rng(SEED)
                r = run_batch(
                    sample_batch(rng, insertion_km), grids[cond], dt=dt,
                    t_max_s=T_MAX_S, density_scale=scale,
                    safe_mode_exit_s=exit_h * 3600.0,
                )
                fracs.append(r.survival_fraction)
            results[(cond, scale)] = np.array(fracs)
            print(f"  sweep 3 {cond:5s} x{scale:.2f} done")
    return results


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------

COLOURS = {"storm": "#d1495b", "quiet": "#2e86ab"}


def _band_plot(ax, x, results, xlabel, title):
    for cond in CONDITIONS:
        lo = np.minimum(results[(cond, DENSITY_SCALES[0])],
                        results[(cond, DENSITY_SCALES[1])])
        hi = np.maximum(results[(cond, DENSITY_SCALES[0])],
                        results[(cond, DENSITY_SCALES[1])])
        ax.fill_between(x, lo, hi, color=COLOURS[cond], alpha=0.22, linewidth=0)
        ax.plot(x, results[(cond, 1.00)], color=COLOURS[cond], lw=1.9,
                label=f"{cond}, NRLMSIS as-is")
        ax.plot(x, results[(cond, 1.19)], color=COLOURS[cond], lw=1.2, ls="--",
                label=f"{cond}, density x1.19")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("survival fraction (reached 550 km shell)")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(-0.03, 1.05)
    ax.grid(alpha=0.25, lw=0.5)


def _unpack(block: dict) -> dict:
    """JSON's "storm|1.00" keys back into (condition, scale) tuples."""
    out = {}
    for key, values in block.items():
        cond, scale = key.split("|")
        out[(cond, float(scale))] = np.asarray(values, dtype=float)
    return out


def plot_from_payload(payload: dict, suffix: str = "") -> list[Path]:
    """Draw all three sweeps from a plain data dict.

    Takes only JSON-shaped data, never live grids or result objects, so the
    Phase 4 read-back check (`sim/export.replot_sweeps_from_json`) exercises the
    same code path as the original run. Nothing is hardcoded: axis ranges,
    annotations and the critical-altitude bands all come from `payload`.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(exist_ok=True)
    written = []
    n_sats = payload["meta"]["n_satellites"]
    crit = payload.get("critical_altitude_km", {})

    # --- sweep 1 ---------------------------------------------------------
    blk = payload["sweep_insertion_altitude"]
    x, s1 = blk["x_km"], _unpack(blk["survival"])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _band_plot(ax, x, s1, "insertion altitude (km)",
               f"Survival vs insertion altitude — {n_sats} satellites, "
               f"thrusters nominal\nCd 2.0–2.4, ram area 1.00–4.48 m², "
               f"Feb 2022 epoch (storm) vs ap=5 (quiet)")
    for cond in CONDITIONS:
        if cond not in crit:
            continue
        lo, hi = crit[cond]
        ax.axvspan(lo, hi, color=COLOURS[cond], alpha=0.13, linewidth=0)
        ax.axvline(hi, color=COLOURS[cond], lw=1.0, ls=":")
        ax.annotate(
            f"analytic $h_{{crit}}$ ({cond})\n{lo:.1f}–{hi:.1f} km",
            xy=(hi, 0.5), xytext=(hi + 12, 0.30 if cond == "storm" else 0.12),
            fontsize=7.5, color=COLOURS[cond],
            arrowprops=dict(arrowstyle="->", color=COLOURS[cond], lw=0.8),
        )
    insertion = payload["meta"].get("actual_insertion_km", 210.0)
    ax.axvline(insertion, color="black", lw=0.9, ls="-.", alpha=0.7)
    ax.text(insertion + 1, 0.02, f"actual insertion\n{insertion:.0f} km", fontsize=7.5)
    ax.set_xlim(min(x), 260)
    ax.legend(fontsize=8, loc="lower right")
    p = OUT / f"sweep_insertion_altitude{suffix}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); written.append(p)

    # --- sweep 2 ---------------------------------------------------------
    blk = payload["sweep_ram_area"]
    x, s2 = blk["x_m2"], _unpack(blk["survival"])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _band_plot(ax, x, s2, "nominal ram area (m²)",
               f"Survival vs ram area — {n_sats} satellites inserted at "
               f"{insertion:.0f} km, thrusters nominal\n"
               "the disputed parameter, quantified")
    lo, hi = blk.get("published_range_m2", list(RAM_AREA_RANGE))
    ax.axvspan(lo, hi, color="grey", alpha=0.15, linewidth=0)
    ax.text((lo + hi) / 2, 0.05, f"Baruah range\n{lo:.2f}–{hi:.2f} m²",
            fontsize=7.5, ha="center")
    for (a, b), colour, label in blk.get("source_ranges_m2", []):
        ax.axvspan(a, b, color=colour, alpha=0.12, linewidth=0)
        ax.text((a + b) / 2, 0.93, label, fontsize=7, ha="center")
    ax.legend(fontsize=8, loc="center left")
    p = OUT / f"sweep_ram_area{suffix}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); written.append(p)

    # --- sweep 3 ---------------------------------------------------------
    blk = payload["sweep_safe_mode_timing"]
    x, s3 = blk["x_hours"], _unpack(blk["survival"])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _band_plot(ax, x, s3, "safe-mode exit time (hours after deployment)",
               f"Survival vs safe-mode exit time — {n_sats} satellites at "
               f"{insertion:.0f} km\nhow late could the fleet have been recovered?")
    ax.legend(fontsize=8, loc="lower left")
    p = OUT / f"sweep_safe_mode_timing{suffix}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig); written.append(p)

    return written


def _serialise(results: dict) -> dict:
    return {f"{cond}|{scale:.2f}": v.tolist() for (cond, scale), v in results.items()}


def main(dt: float = DT_S) -> dict:
    sw = SpaceWeather.load(DATA / "SW-All.csv")
    OUT.mkdir(exist_ok=True)
    grids = _grids(sw)

    print("PHYSICS.md §9 sweeps — 49 satellites, "
          f"dt = {dt:.0f} s, t_max = {T_MAX_S / 86400:.0f} days")
    s1 = sweep_insertion_altitude(grids, dt)
    s2 = sweep_ram_area(grids, dt=dt)
    s3 = sweep_safe_mode_timing(grids, dt=dt)

    from .export import ATMOSPHERE_MODEL, sim_version, write_json

    mid_area, mid_cd = float(np.mean(RAM_AREA_RANGE)), float(np.mean(CD_RANGE))
    payload = {
        "meta": {
            "generated": datetime.now(timezone.utc).isoformat(),
            "sim_version": sim_version(),
            "n_satellites": N_SATELLITES,
            "epoch": EPOCH.isoformat(),
            "dt_s": dt,
            "t_max_s": T_MAX_S,
            "density_scales": list(DENSITY_SCALES),
            "target_shell_km": 550.0,
            "actual_insertion_km": 210.0,
            "atmosphere_model": ATMOSPHERE_MODEL,
            "quiet_definition": "same epoch and F10.7, ap forced to 5 (PHYSICS.md §4.1)",
            "thruster_mode": "nominal from deployment (sweeps 1-2); swept (sweep 3)",
        },
        "critical_altitude_km": {
            cond: list(critical_altitude_band(grids[cond], mid_cd, mid_area))
            for cond in CONDITIONS
        },
        "sweep_insertion_altitude": {
            "x_km": ALTITUDES_KM, "survival": _serialise(s1)},
        "sweep_ram_area": {
            "x_m2": RAM_AREAS_M2,
            "survival": _serialise(s2),
            "published_range_m2": list(RAM_AREA_RANGE),
            "source_ranges_m2": [
                [[0.5, 1.5], "#8ac926", "Source A\n0.5–1.5 m²"],
                [[5.0, 6.0], "#ff924c", "Source B\n5–6 m²"],
            ],
        },
        "sweep_safe_mode_timing": {
            "x_hours": SAFE_MODE_EXIT_H, "survival": _serialise(s3)},
    }
    write_json(OUT / "sweeps.json", payload)
    plot_from_payload(payload)

    _report(grids, s1, s2, s3)
    return payload


def _report(grids, s1, s2, s3) -> None:
    def crossing(x, y, level=0.5):
        """First x at which survival rises through `level`."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        for i in range(1, len(x)):
            if y[i - 1] < level <= y[i]:
                if y[i] == y[i - 1]:
                    return x[i]
                return x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
        return float("nan")

    print("\n" + "=" * 78)
    print("SWEEP 1 — insertion altitude (50% survival crossing)")
    print("=" * 78)
    mid_area, mid_cd = float(np.mean(RAM_AREA_RANGE)), float(np.mean(CD_RANGE))
    for cond in CONDITIONS:
        lo, hi = critical_altitude_band(grids[cond], mid_cd, mid_area)
        for scale in DENSITY_SCALES:
            print(f"  {cond:5s} x{scale:.2f}: 50% at {crossing(ALTITUDES_KM, s1[(cond, scale)]):.1f} km")
        print(f"        analytic h_crit (Cd={mid_cd}, A={mid_area:.2f} m²): "
              f"{lo:.1f}–{hi:.1f} km")
    st = crossing(ALTITUDES_KM, s1[("storm", 1.00)])
    qt = crossing(ALTITUDES_KM, s1[("quiet", 1.00)])
    sc = crossing(ALTITUDES_KM, s1[("storm", 1.19)])
    print(f"\n  storm vs quiet shift        : {st - qt:+.1f} km")
    print(f"  density x1.19 shift (storm) : {sc - st:+.1f} km")
    print(f"  -> the density-model uncertainty moves the threshold "
          f"{abs(sc - st) / max(abs(st - qt), 1e-9):.1f}x further than the storm does")

    def falling_crossing(x, y, level=0.5):
        """First x at which a falling survival curve drops through `level`."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        for i in range(1, len(x)):
            if y[i - 1] >= level > y[i]:
                return x[i - 1] + (y[i - 1] - level) * (x[i] - x[i - 1]) / (y[i - 1] - y[i])
        return None

    print("\n" + "=" * 78)
    print("SWEEP 2 — ram area at 210 km")
    print("=" * 78)
    for cond in CONDITIONS:
        for scale in DENSITY_SCALES:
            y = s2[(cond, scale)]
            drop = falling_crossing(RAM_AREAS_M2, y)
            where = f"50% at {drop:.2f} m²" if drop else "never falls below 50%"
            print(f"  {cond:5s} x{scale:.2f}: survival {y[0]:.2f} -> {y[-1]:.2f} "
                  f"across 0.5-6.0 m²; {where}")

    print("\n" + "=" * 78)
    print("SWEEP 3 — safe-mode exit time at 210 km")
    print("=" * 78)
    for cond in CONDITIONS:
        for scale in DENSITY_SCALES:
            y = s3[(cond, scale)]
            xs = np.asarray(SAFE_MODE_EXIT_H, float)
            last = falling_crossing(xs, y)
            where = f"50% at {last:.1f} h" if last else "never falls below 50%"
            print(f"  {cond:5s} x{scale:.2f}: survival {y[0]:.2f} at 0 h -> "
                  f"{y[-1]:.2f} at {xs[-1]:.0f} h; {where}")
    print()


if __name__ == "__main__":
    main()
