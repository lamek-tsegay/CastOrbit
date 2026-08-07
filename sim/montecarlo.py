"""Monte Carlo batches and the outcome taxonomy.

PHYSICS.md §9. A batch is 49 satellites, each drawing its own parameters.

The ensemble is propagated *vectorised*: the state is a (2, n) array
`[a_vector, m_vector]` and every satellite advances together through the same
`rk4_step` and the same `derivatives` used by the scalar validation runs. There
is no second integrator and no second copy of the equations -- only a mask that
freezes satellites which have terminated or have not yet deployed. That keeps
the "did you write the integrator?" answer to exactly one function
(ARCHITECTURE.md §4) while making 10^7 derivative evaluations tractable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

import numpy as np

from .atmosphere import DensityGrid, SpaceWeather
from .constants import REENTRY_ALTITUDE, R_E
from .dynamics import derivatives
from .integrator import rk4_step
from .satellite import Outcome

# --- §9 sampling distributions --------------------------------------------
N_SATELLITES = 49
RAM_AREA_RANGE = (1.00, 4.48)        # uniform, published bounding range
MASS_MEAN_KG, MASS_SIGMA_FRAC = 227.0, 0.03
INSERTION_SIGMA_KM = 2.0
DEPLOY_SPREAD_S = 3600.0             # uniform +/- 30 min, see `sample_batch`
THRUST_RATED_N, THRUST_SIGMA_FRAC = 0.071, 0.02   # 0.071 N is DERIVED
CD_RANGE = (2.0, 2.4)                # uniform, free-molecular flow uncertainty

ISP_S = 1666.0                       # DERIVED (satellite_specs.json v1_5.isp_s)
KNIFE_EDGE_AREA_M2 = 1.00            # safe-mode ram area (§5, "estimated")
TARGET_SHELL_M = 550e3               # ARCHITECTURE.md §6 JSON contract

# PHYSICS.md §4.1 convention for the counterfactual quiet atmosphere.
QUIET_AP = 5.0


@dataclass
class Batch:
    """One drawn ensemble. Arrays are all length n."""

    area_m2: np.ndarray
    mass_kg: np.ndarray
    insertion_altitude_m: np.ndarray
    deploy_time_s: np.ndarray
    thrust_n: np.ndarray
    cd: np.ndarray

    def __len__(self) -> int:
        return int(self.area_m2.size)


@dataclass
class BatchResult:
    outcomes: np.ndarray                  # array of Outcome
    outcome_time_s: np.ndarray
    final_altitude_m: np.ndarray
    t_hist: np.ndarray = field(repr=False)
    h_hist_km: np.ndarray = field(repr=False)   # (n_samples, n_sats)
    batch: Batch = field(repr=False, default=None)

    def counts(self) -> dict[str, int]:
        return {
            o.value: int(np.sum(self.outcomes == o))
            for o in Outcome
        }

    @property
    def survival_fraction(self) -> float:
        """Fraction reaching the operational shell. PHYSICS.md §9."""
        return float(np.mean(self.outcomes == Outcome.REACHED_SHELL))

    @property
    def reentry_fraction(self) -> float:
        return float(np.mean(self.outcomes == Outcome.REENTERED))


def sample_batch(
    rng: np.random.Generator,
    insertion_altitude_km: float = 210.0,
    n: int = N_SATELLITES,
    area_m2: float | None = None,
    area_range: tuple[float, float] = RAM_AREA_RANGE,
) -> Batch:
    """Draw one batch per the PHYSICS.md §9 distribution table.

    `deploy_time_s` is drawn uniform on [0, 3600] s rather than +/-1800 s about
    the epoch. These are the same thing: the simulation clock is started half an
    hour before nominal deployment (see `run_batch`), so the spread is +/-30 min
    about the nominal time while every query stays inside the density grid.

    `area_m2` overrides the ram-area draw with a fixed value, which is what the
    ram-area sweep needs (PHYSICS.md §9, sweep 2).
    """
    if area_m2 is None:
        area = rng.uniform(area_range[0], area_range[1], n)
    else:
        area = np.full(n, float(area_m2))
    return Batch(
        area_m2=area,
        mass_kg=rng.normal(MASS_MEAN_KG, MASS_SIGMA_FRAC * MASS_MEAN_KG, n),
        insertion_altitude_m=rng.normal(
            insertion_altitude_km * 1e3, INSERTION_SIGMA_KM * 1e3, n
        ),
        deploy_time_s=rng.uniform(0.0, DEPLOY_SPREAD_S, n),
        thrust_n=rng.normal(
            THRUST_RATED_N, THRUST_SIGMA_FRAC * THRUST_RATED_N, n
        ),
        cd=rng.uniform(CD_RANGE[0], CD_RANGE[1], n),
    )


def build_grid(
    sw: SpaceWeather,
    epoch: datetime,
    storm: bool,
    duration_s: float,
    lat_deg: float = 53.22,
) -> DensityGrid:
    """Storm or quiet atmosphere over the batch window.

    "Storm" is the real Feb 2022 SW-All.csv history. "Quiet" is the same epoch
    and the same solar flux with geomagnetic activity forced to ap = 5, the
    PHYSICS.md §4.1 quiet value. Holding F10.7 fixed is deliberate: the gap
    between the two curves is then attributable to geomagnetic activity alone,
    which is what the §9 sweep is asking about.
    """
    return DensityGrid(
        epoch,
        sw,
        lat_deg=lat_deg,
        duration_s=duration_s,
        storm_time=True,
        ap_override=None if storm else QUIET_AP,
    )


def run_batch(
    batch: Batch,
    grid: DensityGrid,
    dt: float = 60.0,
    t_max_s: float = 15 * 86400.0,
    safe_mode_exit_s: float | None = None,
    density_scale: float = 1.0,
    target_shell_m: float = TARGET_SHELL_M,
    sample_every: int = 60,

) -> BatchResult:
    """Propagate the whole ensemble. PHYSICS.md §3.3 termination, §5 safe mode.

    Args:
        safe_mode_exit_s: None means the thrusters are on from deployment
            (NOMINAL throughout). A finite value means SAFE_MODE -- F = 0 and
            the knife-edge ram area -- until that time, then NOMINAL. `inf`
            means safe mode is never exited, which is what actually happened in
            February 2022.
        density_scale: uniform multiplier on rho. Used to carry the validated
            density uncertainty band (see `sim/validate.py`), never to tune.

    Propellant exhaustion is *not* modelled: no propellant mass is published for
    the v1.5 bus in `data/satellite_specs.json`, and inventing one would put a
    fabricated number in a headline result. Mass is depleted correctly via
    dm/dt with the DERIVED Isp, but there is no dry-mass floor, so the
    PROPELLANT_EXHAUSTED outcome is unreachable by construction and always
    reports zero. That is a stated limitation, not an observation.
    """
    n = len(batch)
    a0 = R_E + batch.insertion_altitude_m
    y = np.vstack([a0.copy(), batch.mass_kg.copy()])

    knife_area = np.full(n, KNIFE_EDGE_AREA_M2)
    zero_thrust = np.zeros(n)
    exit_s = np.inf if safe_mode_exit_s is None else float(safe_mode_exit_s)
    never_safe = safe_mode_exit_s is None

    def deriv(t: float, state: np.ndarray) -> np.ndarray:
        if never_safe or t >= exit_s:
            thrust, area = batch.thrust_n, batch.area_m2
        else:
            thrust, area = zero_thrust, knife_area
        rho = density_scale * grid.lookup(t, state[0] - R_E)
        return derivatives(state, rho, thrust, batch.cd, area, isp=ISP_S)

    outcomes = np.full(n, Outcome.INDETERMINATE, dtype=object)
    outcome_time = np.full(n, np.nan)
    active = np.ones(n, dtype=bool)

    t_hist, h_hist = [0.0], [(y[0] - R_E) / 1e3]
    n_steps = int(np.ceil(t_max_s / dt))

    for step in range(1, n_steps + 1):
        y_new = rk4_step(deriv, (step - 1) * dt, y, dt)
        t = step * dt

        # Satellites that have terminated, or have not yet been deployed, keep
        # their previous state.
        moving = active & (t > batch.deploy_time_s)
        y = np.where(moving, y_new, y)

        assert np.all(y[0] > 0.0), f"non-positive semi-major axis at t={t}"
        assert np.all(y[1] > 0.0), f"non-positive mass at t={t}"

        h = y[0] - R_E
        reentered = active & (h < REENTRY_ALTITUDE)
        reached = active & (h >= target_shell_m)

        for mask, outcome in ((reentered, Outcome.REENTERED),
                              (reached, Outcome.REACHED_SHELL)):
            if mask.any():
                outcomes[mask] = outcome
                outcome_time[mask] = t
                active &= ~mask

        if step % sample_every == 0 or step == n_steps:
            t_hist.append(t)
            h_hist.append(h / 1e3)

        if not active.any():
            break

    return BatchResult(
        outcomes=outcomes,
        outcome_time_s=outcome_time,
        final_altitude_m=y[0] - R_E,
        t_hist=np.array(t_hist),
        h_hist_km=np.array(h_hist),
        batch=batch,
    )
