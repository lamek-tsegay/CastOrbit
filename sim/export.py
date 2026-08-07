"""JSON emission for the frontend, and the read-back check.

ARCHITECTURE.md §6. Two files, deliberately separate rather than one stretched
schema:

  out/batch.json   the fleet reproduction: both Cd conventions run side by
                   side, each satellite tagged with its own effective drag
                   parameter Cd*A, plus the critical-altitude time series
                   the altitude chart overlays.
  out/sweeps.json  the three §9 sweep curves.

The Phase 4 gate is that a Python script can read these back and reproduce the
Phase 3 plots. That is only a real check if the plotting code cannot see the
in-memory objects, so `sim/sweeps.py` plots from a plain payload dict and
`replot_sweeps_from_json` feeds it one parsed straight off disk; the same
applies to `replot_batch_from_json` here.
"""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .constants import R_E
from .critical import critical_altitude
from .groundtrack import cause_of_loss, ground_track
from .satellite import Outcome

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out"

ATMOSPHERE_MODEL = "NRLMSIS 2.1 via pymsis 0.12.0"

# The two fleet-reproduction conventions compared throughout README.md's
# central finding: this project's own baseline Cd (satellite_specs.json,
# aerodynamics.drag_coefficient_baseline) against Baruah et al.'s own
# validation convention. Both keep the drawn ram area for the whole window
# (data/event_feb2022.json: 1.00-4.48 m2 IS the safe-mode bounding range, not
# a nominal-flight range) -- see sim/validate.py:run_fleet.
FLEET_RUNS = (
    ("cd2.2", 2.2, "this project's own baseline Cd (aerodynamics convention)"),
    ("cd1.0", 1.0, "Baruah et al.'s validation convention"),
)


def sim_version() -> str:
    """Current git sha, for `meta.sim_version` (ARCHITECTURE.md §6)."""
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT, capture_output=True, text=True, timeout=5,
        )
        if sha.returncode == 0:
            dirty = subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=ROOT, capture_output=True, text=True, timeout=5,
            )
            suffix = "-dirty" if dirty.stdout.strip() else ""
            return sha.stdout.strip() + suffix
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def write_json(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def critical_altitude_series(
    grid, thrust_n: float, cd: float, area_m2: float,
    t_s: np.ndarray, epoch: datetime,
) -> dict:
    """Critical altitude over time. PHYSICS.md §4, ARCHITECTURE.md §6.

    For a fleet in safe mode this is a counterfactual and must be labelled as
    one: with F = 0 there is no balance point at all. What is drawn is the
    altitude at which thrust *would* balance drag if orbit-raising resumed --
    i.e. the altitude below which the fleet can no longer save itself. That is
    what makes the trajectories crossing it meaningful.
    """
    times, h_crit = [], []
    for t in t_s:
        try:
            h = critical_altitude(
                lambda hh: float(grid.lookup(float(t), np.array([hh]))[0]),
                thrust_n, cd, area_m2,
            )
            h_crit.append(round(h / 1e3, 4))
        except ValueError:
            h_crit.append(None)
        times.append((epoch + timedelta(seconds=float(t))).isoformat())
    return {
        "_note": (
            "Counterfactual: the fleet is in safe mode with F = 0, where no "
            "balance point exists. This is the altitude at which the rated "
            "thrust would balance drag if orbit-raising were resumed, i.e. the "
            "recovery boundary, evaluated at this run's own Cd and the "
            "median drawn ram area."
        ),
        "thrust_n": thrust_n,
        "cd": cd,
        "ram_area_m2": area_m2,
        "times": times,
        "h_crit_km": h_crit,
    }


def _build_run(
    result, grid, epoch: datetime, label: str, cd: float, description: str,
    density_scale: float, ram_area_range_m2: tuple[float, float],
    thrust_rated_n: float, sample_stride_s: float = 600.0,
) -> dict:
    """One fleet-reproduction run: config, outcomes, critical altitude, and
    every satellite's downsampled trajectory tagged with its own Cd*A.
    """
    batch = result.batch
    t_hist, h_hist = result.t_hist, result.h_hist_km

    if t_hist.size > 1:
        step = max(1, int(round(sample_stride_s / (t_hist[1] - t_hist[0]))))
    else:
        step = 1
    idx = np.arange(0, t_hist.size, step)
    if idx[-1] != t_hist.size - 1:
        idx = np.append(idx, t_hist.size - 1)

    t_ds = t_hist[idx]
    rho_ds = np.array([grid.lookup(float(t), h_hist[i] * 1e3)
                       for i, t in zip(idx, t_ds)])

    n_sats = len(batch)
    satellites = []
    for k in range(n_sats):
        outcome = result.outcomes[k]
        outcome_value = outcome.value if isinstance(outcome, Outcome) else str(outcome)
        t_out = result.outcome_time_s[k]
        if np.isfinite(t_out):
            keep = t_ds <= t_out
            if not keep.any():
                keep[0] = True
        else:
            keep = np.ones(t_ds.size, dtype=bool)
        area = float(batch.area_m2[k])
        outcome_time_iso = (
            (epoch + timedelta(seconds=float(t_out))).isoformat()
            if np.isfinite(t_out) else None
        )
        t_kept = t_ds[keep]
        h_kept = h_hist[idx][keep, k]
        lat_deg, lon_deg = ground_track(
            t_kept, h_kept, satellite_id=k, n_satellites=n_sats,
            deploy_time_s=float(batch.deploy_time_s[k]),
        )
        satellites.append({
            "id": k,
            "params": {
                "mass_kg": round(float(batch.mass_kg[k]), 4),
                "ram_area_m2": round(area, 4),
                "cd": round(float(cd), 4),
                "cd_times_area_m2": round(cd * area, 4),
                "thrust_n": 0.0,  # safe mode throughout
                "rated_thrust_n": round(float(batch.thrust_n[k]), 6),
                "insertion_altitude_km": round(
                    float(batch.insertion_altitude_m[k]) / 1e3, 4),
                "deploy_time_s": round(float(batch.deploy_time_s[k]), 1),
            },
            "outcome": outcome_value,
            "outcome_time": outcome_time_iso,
            "cause": cause_of_loss(outcome_value, outcome_time_iso, cd * area),
            "trajectory": {
                "t_s": [round(float(v), 1) for v in t_kept],
                "h_km": [round(float(v), 4) for v in h_kept],
                "rho": [float(f"{v:.6g}") for v in rho_ds[keep, k]],
                "lat_deg": lat_deg,
                "lon_deg": lon_deg,
            },
        })

    counts = result.counts()
    lo, hi = ram_area_range_m2
    return {
        "label": label,
        "description": description,
        "config": {
            "cd": cd,
            "density_scale": density_scale,
            "ram_area_range_m2": [lo, hi],
            "effective_drag_range_m2": [round(cd * lo, 4), round(cd * hi, 4)],
        },
        "outcome_counts": counts,
        "critical_altitude": critical_altitude_series(
            grid, thrust_rated_n, cd, float(np.median(batch.area_m2)), t_ds, epoch
        ),
        "satellites": satellites,
    }


def export_fleet_batch(sw, path: Path | None = None, dt: float = 30.0) -> Path:
    """Emit the fleet reproduction, both Cd conventions, to out/batch.json.

    README.md's central finding is that Cd*A is what's actually observable,
    not Cd or A separately -- so every satellite here carries its own
    `cd_times_area_m2` rather than making a consumer recompute it, and each
    run's config states the Cd*A range that run spans. Both runs keep the
    published, unscaled 1.00-4.48 m2 ram-area range: this file reports what
    the simulator produces, not the reconciled range solved for in
    `sim/validate.py:fleet_reproduction` (which is a diagnostic, reported in
    README.md, not a run of the simulator).
    """
    from .montecarlo import THRUST_RATED_N
    from .validate import FLEET_END, FLEET_LOST, FLEET_N, FLEET_SURVIVED, run_fleet

    grid = None
    epoch = None
    runs = []
    for label, cd, description in FLEET_RUNS:
        result, grid, t_max, epoch = run_fleet(
            sw, density_scale=1.0, area_range=(1.00, 4.48), cd=cd, dt=dt, grid=grid,
        )
        runs.append(_build_run(
            result, grid, epoch, label, cd, description,
            density_scale=1.0, ram_area_range_m2=(1.00, 4.48),
            thrust_rated_n=THRUST_RATED_N,
        ))

    payload = {
        "meta": {
            "scenario": "fleet_reproduction_feb2022",
            "generated": datetime.now(timezone.utc).isoformat(),
            "sim_version": sim_version(),
            "atmosphere_model": ATMOSPHERE_MODEL,
            "integrator": "hand-written RK4, fixed step",
            "epoch": epoch.isoformat(),
            "window_end": FLEET_END.isoformat(),
            "n_satellites": FLEET_N,
            "insertion_altitude_km": 210.0,
            "reentry_altitude_km": 100.0,
            "thruster_mode": "safe_mode for the whole window (F = 0, never exits)",
        },
        "observed": {
            "lost": FLEET_LOST,
            "survived": FLEET_SURVIVED,
            "source": "data/event_feb2022.json",
        },
        "runs": runs,
    }
    return write_json(path or (OUT / "batch.json"), payload)


def replot_sweeps_from_json(path: Path | None = None) -> list[Path]:
    """Phase 4 gate: rebuild the Phase 3 sweep plots from the JSON alone."""
    from .sweeps import plot_from_payload

    payload = load_json(path or (OUT / "sweeps.json"))
    return plot_from_payload(payload, suffix="_from_json")


def replot_batch_from_json(path: Path | None = None) -> Path:
    """Fleet altitude chart straight from batch.json -- nothing but the file.

    One panel per run, sharing a y-axis, so the two Cd conventions are visually
    comparable: this is the picture behind "Cd=2.2 loses everyone, Cd=1.0
    brackets reality" in README.md's central finding.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = load_json(path or (OUT / "batch.json"))
    colours = {
        "REACHED_SHELL": "#2a9d3f",
        "REENTERED": "#d1495b",
        "PROPELLANT_EXHAUSTED": "#e8a33d",
        "INDETERMINATE": "#8d99ae",
    }
    epoch = datetime.fromisoformat(payload["meta"]["epoch"])
    observed = payload.get("observed")

    runs = payload["runs"]
    fig, axes = plt.subplots(1, len(runs), figsize=(6.5 * len(runs), 5.5),
                             sharey=True)
    if len(runs) == 1:
        axes = [axes]

    for ax, run in zip(axes, runs):
        for sat in run["satellites"]:
            tr = sat["trajectory"]
            ax.plot(np.array(tr["t_s"]) / 3600.0, tr["h_km"],
                    color=colours.get(sat["outcome"], "grey"),
                    lw=0.9, alpha=0.75)

        crit = run["critical_altitude"]
        h_crit = np.array([np.nan if v is None else v for v in crit["h_crit_km"]],
                          dtype=float)
        t_crit = np.array([
            (datetime.fromisoformat(s) - epoch).total_seconds()
            for s in crit["times"]
        ])
        ax.plot(t_crit / 3600.0, h_crit, color="black", lw=1.5, ls="--",
                label="recovery boundary (Cd*A this run)")

        ax.axhline(100.0, color="black", lw=0.8, ls=":")
        counts = run["outcome_counts"]
        lost = counts.get("REENTERED", 0)
        cfg = run["config"]
        title = (
            f"Cd = {cfg['cd']}  (Cd*A ∈ [{cfg['effective_drag_range_m2'][0]:.2f}, "
            f"{cfg['effective_drag_range_m2'][1]:.2f}] m²)\n"
            f"{lost} lost / {payload['meta']['n_satellites'] - lost} survived"
        )
        if observed:
            title += f"  (observed {observed['lost']}/{observed['survived']})"
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel(f"hours after {payload['meta']['epoch'][:16].replace('T', ' ')} UT")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=7.5, loc="lower left")

    axes[0].set_ylabel("altitude (km)")
    axes[0].set_ylim(90, 220)
    fig.suptitle(
        f"{payload['meta']['n_satellites']} satellites, safe mode throughout, "
        "both Cd conventions", fontsize=11,
    )
    fig.tight_layout()
    out = OUT / "batch_altitude_from_json.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out
