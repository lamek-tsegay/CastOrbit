"""JSON emission for the frontend, and the read-back check.

ARCHITECTURE.md §6. Two files, deliberately separate rather than one stretched
schema:

  out/batch.json   one Monte Carlo batch: per-satellite parameters, outcome and
                   downsampled trajectory, plus the critical-altitude time
                   series that the altitude chart overlays.
  out/sweeps.json  the three §9 sweep curves.

The Phase 4 gate is that a Python script can read these back and reproduce the
Phase 3 plots. That is only a real check if the plotting code cannot see the
in-memory objects, so `sim/sweeps.py` plots from a plain payload dict and
`replot_sweeps_from_json` feeds it one parsed straight off disk.
"""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .constants import R_E
from .critical import critical_altitude
from .satellite import Outcome

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out"

ATMOSPHERE_MODEL = "NRLMSIS 2.1 via pymsis 0.12.0"


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
            "recovery boundary."
        ),
        "thrust_n": thrust_n,
        "cd": cd,
        "ram_area_m2": area_m2,
        "times": times,
        "h_crit_km": h_crit,
    }


def export_batch(
    result,
    grid,
    epoch: datetime,
    scenario: str,
    cd: float,
    density_scale: float,
    insertion_altitude_km: float,
    target_shell_km: float,
    thrust_rated_n: float,
    observed: dict | None = None,
    path: Path | None = None,
    sample_stride_s: float = 600.0,
) -> Path:
    """Emit one batch to JSON in the ARCHITECTURE.md §6 shape."""
    batch = result.batch
    t_hist, h_hist = result.t_hist, result.h_hist_km

    # Downsample to ~10 min spacing for display; full resolution stays in Python.
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

    satellites = []
    for k in range(len(batch)):
        outcome = result.outcomes[k]
        t_out = result.outcome_time_s[k]
        # Trajectory is truncated at the satellite's own termination.
        if np.isfinite(t_out):
            keep = t_ds <= t_out
            if not keep.any():
                keep[0] = True
        else:
            keep = np.ones(t_ds.size, dtype=bool)
        satellites.append({
            "id": k,
            "params": {
                "mass_kg": round(float(batch.mass_kg[k]), 4),
                "ram_area_m2": round(float(batch.area_m2[k]), 4),
                "thrust_n": 0.0,  # safe mode throughout
                "rated_thrust_n": round(float(batch.thrust_n[k]), 6),
                "cd": round(float(cd), 4),
                "insertion_altitude_km": round(
                    float(batch.insertion_altitude_m[k]) / 1e3, 4),
                "deploy_time_s": round(float(batch.deploy_time_s[k]), 1),
            },
            "outcome": outcome.value if isinstance(outcome, Outcome) else str(outcome),
            "outcome_time": (
                (epoch + timedelta(seconds=float(t_out))).isoformat()
                if np.isfinite(t_out) else None
            ),
            "trajectory": {
                "t_s": [round(float(v), 1) for v in t_ds[keep]],
                "h_km": [round(float(v), 4) for v in h_hist[idx][keep, k]],
                "rho": [float(f"{v:.6g}") for v in rho_ds[keep, k]],
            },
        })

    counts = result.counts()
    payload = {
        "meta": {
            "scenario": scenario,
            "generated": datetime.now(timezone.utc).isoformat(),
            "sim_version": sim_version(),
            "atmosphere_model": ATMOSPHERE_MODEL,
            "integrator": "hand-written RK4, fixed step",
        },
        "config": {
            "epoch": epoch.isoformat(),
            "n_satellites": len(batch),
            "cd": cd,
            "density_scale": density_scale,
            "insertion_altitude_km": insertion_altitude_km,
            "target_shell_km": target_shell_km,
            "thruster_mode": "safe_mode for the whole window (F = 0, never exits)",
            "reentry_altitude_km": 100.0,
        },
        "outcome_counts": counts,
        "critical_altitude": critical_altitude_series(
            grid, thrust_rated_n, cd, float(np.median(batch.area_m2)), t_ds, epoch
        ),
        "satellites": satellites,
    }
    if observed is not None:
        payload["validation"] = observed

    return write_json(path or (OUT / "batch.json"), payload)


def replot_sweeps_from_json(path: Path | None = None) -> list[Path]:
    """Phase 4 gate: rebuild the Phase 3 sweep plots from the JSON alone."""
    from .sweeps import plot_from_payload

    payload = load_json(path or (OUT / "sweeps.json"))
    return plot_from_payload(payload, suffix="_from_json")


def replot_batch_from_json(path: Path | None = None) -> Path:
    """Altitude chart straight from batch.json -- nothing but the file."""
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

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for sat in payload["satellites"]:
        tr = sat["trajectory"]
        ax.plot(np.array(tr["t_s"]) / 3600.0, tr["h_km"],
                color=colours.get(sat["outcome"], "grey"), lw=0.9, alpha=0.75)

    crit = payload["critical_altitude"]
    t_h = [i * 0 for i in range(0)]  # placeholder, replaced below
    t_axis = np.array(payload["satellites"][0]["trajectory"]["t_s"])
    h_crit = np.array([np.nan if v is None else v for v in crit["h_crit_km"]],
                      dtype=float)
    n = min(h_crit.size, len(crit["times"]))
    t_crit = np.linspace(0, t_axis.max() if t_axis.size else 0, n)
    # The critical-altitude series shares the batch's downsampled time base.
    t_full = np.array([
        (datetime.fromisoformat(s) - datetime.fromisoformat(payload["config"]["epoch"])
         ).total_seconds() for s in crit["times"][:n]
    ])
    ax.plot(t_full / 3600.0, h_crit[:n], color="black", lw=1.6, ls="--",
            label="recovery boundary (critical altitude if thrust resumed)")

    ax.axhline(100.0, color="black", lw=0.8, ls=":")
    ax.text(1, 103, "100 km — unrecoverable", fontsize=8)
    counts = payload["outcome_counts"]
    ax.set_title(
        f"{payload['config']['n_satellites']} satellites, safe mode throughout "
        f"— Cd = {payload['config']['cd']}, density x{payload['config']['density_scale']}\n"
        f"{counts.get('REENTERED', 0)} reentered / "
        f"{payload['config']['n_satellites'] - counts.get('REENTERED', 0)} survived"
        + (f"  (observed {payload['validation']['lost']}/"
           f"{payload['validation']['survived']})" if "validation" in payload else ""),
        fontsize=10,
    )
    ax.set_xlabel(f"hours after {payload['config']['epoch'][:16].replace('T', ' ')} UT")
    ax.set_ylabel("altitude (km)")
    ax.set_ylim(90, 220)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8, loc="lower left")
    out = OUT / "batch_altitude_from_json.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out
