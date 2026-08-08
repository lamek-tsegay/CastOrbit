"""Hand-written RK4 integrator and the propagation loops.

PHYSICS.md §7. Classical fourth-order Runge-Kutta, written out directly.
`scipy.integrate.solve_ivp` is explicitly forbidden for the main propagation
loop (ARCHITECTURE.md §4).

Two drivers share the one `rk4_step`:

  * `propagate`          -- fixed step. V1. Unchanged, and the reference the
                            adaptive driver is measured against.
  * `propagate_adaptive` -- variable step, same RK4 stages. V2_BRIEF.md §4.

The adaptive driver exists because a 25-year compliance run at the V1 fixed
10 s step is ~8e7 steps per satellite. It does not replace RK4 and does not
change the equations; it only chooses `dt`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .constants import REENTRY_ALTITUDE, R_E
from .satellite import Outcome


def rk4_step(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    y: np.ndarray,
    dt: float,
) -> np.ndarray:
    """One classical RK4 step of the system dy/dt = f(t, y).

    PHYSICS.md §7.
    """
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


class Trajectory:
    """Result of a propagation: sampled state history plus an outcome label."""

    def __init__(
        self,
        t_s: np.ndarray,
        a_m: np.ndarray,
        mass_kg: np.ndarray,
        outcome: Outcome,
        outcome_time_s: float | None,
        stats: "StepStats | None" = None,
    ):
        self.t_s = t_s
        self.a_m = a_m
        self.mass_kg = mass_kg
        self.outcome = outcome
        self.outcome_time_s = outcome_time_s
        # None for the fixed-step driver, where there is nothing to report.
        self.stats = stats

    @property
    def h_m(self) -> np.ndarray:
        """Altitude above the Earth's surface, m. PHYSICS.md §2."""
        return self.a_m - R_E

    @property
    def h_km(self) -> np.ndarray:
        return self.h_m / 1e3

    def __repr__(self) -> str:
        return (
            f"Trajectory(outcome={self.outcome.value}, "
            f"t_end={self.t_s[-1]:.0f}s, h_end={self.h_km[-1]:.3f}km, "
            f"n={len(self.t_s)})"
        )


def propagate(
    deriv: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    dt: float,
    t_max: float,
    t0: float = 0.0,
    shell_altitude_m: float | None = None,
    reentry_altitude_m: float = REENTRY_ALTITUDE,
    dry_mass_kg: float | None = None,
    sample_every: int = 1,
) -> Trajectory:
    """Integrate [a, m] forward with fixed-step RK4 until a stop condition.

    PHYSICS.md §7 (integration) and §3.3 (termination).

    Termination follows §3.3: reentry below `reentry_altitude_m`, success at
    or above `shell_altitude_m`, otherwise INDETERMINATE at `t_max`.
    Propellant exhaustion does not stop the run -- per §3.3 it is the caller's
    job to switch to F = 0 and let the outcome resolve later -- but if the
    state reaches `dry_mass_kg` and nothing else terminates the run, the
    outcome is reported as PROPELLANT_EXHAUSTED.

    `deriv(t, y)` must return d[a, m]/dt. Density lookup, thruster mode
    scheduling and mass-flow cutoff all live inside that closure, which keeps
    this loop pure.
    """
    if dt <= 0:
        raise ValueError("dt must be positive")
    n_steps = int(np.ceil((t_max - t0) / dt))

    t = float(t0)
    y = np.asarray(y0, dtype=float).copy()

    ts = [t]
    a_hist = [y[0]]
    m_hist = [y[1]]

    outcome = Outcome.INDETERMINATE
    outcome_time: float | None = None
    exhausted = False

    for step in range(1, n_steps + 1):
        y = rk4_step(deriv, t, y, dt)
        t = t0 + step * dt

        # PHYSICS.md §7: both state components must remain positive.
        assert y[0] > 0.0, f"non-positive semi-major axis a={y[0]!r} at t={t}"
        assert y[1] > 0.0, f"non-positive mass m={y[1]!r} at t={t}"

        if dry_mass_kg is not None and y[1] <= dry_mass_kg:
            exhausted = True

        h = y[0] - R_E
        terminated = False
        if h < reentry_altitude_m:
            outcome, outcome_time, terminated = Outcome.REENTERED, t, True
        elif shell_altitude_m is not None and h >= shell_altitude_m:
            outcome, outcome_time, terminated = Outcome.REACHED_SHELL, t, True

        if terminated or step % sample_every == 0 or step == n_steps:
            ts.append(t)
            a_hist.append(y[0])
            m_hist.append(y[1])

        if terminated:
            break
    else:
        if exhausted:
            outcome = Outcome.PROPELLANT_EXHAUSTED

    return Trajectory(
        np.array(ts), np.array(a_hist), np.array(m_hist), outcome, outcome_time
    )


@dataclass
class StepStats:
    """What the adaptive controller actually did. Reported, not assumed.

    A variable-step run is only trustworthy if you can see the step sizes it
    chose, so these travel with the trajectory rather than being logged and
    discarded.
    """

    n_accepted: int = 0
    n_rejected: int = 0
    dt_min_used: float = float("inf")
    dt_max_used: float = 0.0
    max_frac_change: float = 0.0     # largest |da/a| actually taken
    n_at_dt_max: int = 0             # steps capped by dt_max, not by tolerance
    n_at_dt_min: int = 0             # steps floored by dt_min -- tolerance NOT met

    @property
    def tolerance_respected(self) -> bool:
        """False if any step hit the dt floor, i.e. accuracy was traded away."""
        return self.n_at_dt_min == 0

    def as_dict(self) -> dict:
        return {
            "n_accepted": self.n_accepted,
            "n_rejected": self.n_rejected,
            "dt_min_used_s": None if self.n_accepted == 0 else self.dt_min_used,
            "dt_max_used_s": self.dt_max_used,
            "max_frac_change": self.max_frac_change,
            "n_at_dt_max": self.n_at_dt_max,
            "n_at_dt_min": self.n_at_dt_min,
            "tolerance_respected": self.tolerance_respected,
        }


def _crossing_time(t0: float, h0: float, t1: float, h1: float, h_target: float) -> float:
    """Linear interpolation for the time at which h passes `h_target`.

    The fixed-step driver reports the time of the step that *ended* past the
    boundary, which is accurate to one step because its steps are 10 s. The
    adaptive driver's steps can be days long, so the crossing has to be
    resolved inside the step or the reported reentry time is meaningless.
    """
    if h1 == h0:
        return t1
    frac = (h0 - h_target) / (h0 - h1)
    return t0 + frac * (t1 - t0)


def propagate_adaptive(
    deriv: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    t_max: float,
    t0: float = 0.0,
    tol: float = 1e-4,
    dt_min: float = 1.0,
    dt_max: float = 86400.0,
    reject_factor: float = 2.0,
    shell_altitude_m: float | None = None,
    reentry_altitude_m: float = REENTRY_ALTITUDE,
    dry_mass_kg: float | None = None,
) -> Trajectory:
    """Integrate [a, m] with variable-step RK4 until a stop condition.

    V2_BRIEF.md §4. Same `rk4_step`, same `deriv`, same termination rules as
    `propagate` (PHYSICS.md §3.3) -- only the step size is chosen rather than
    fixed.

    **Step control.** `dt` is picked so the fractional change in the
    semi-major axis over the step stays under `tol`:

        dt = tol * a / |da/dt|

    clamped to [`dt_min`, `dt_max`]. This is a *rate*-based controller, not a
    truncation-error estimator: it bounds how far the state moves per step,
    which for this ODE is the quantity that matters, because the physics is
    driven by rho(h) and a step that moves the satellite a long way down the
    density profile is the way this integration goes wrong. RK4's own
    truncation error is then far below the rate bound; that claim is tested by
    halving `tol` and checking the answer does not move
    (`tests/test_adaptive.py`).

    Because the estimate uses the derivative at the *start* of the step, it
    underestimates `dt` during runaway decay, when |da/dt| grows within the
    step. Any step whose achieved |da/a| exceeds `reject_factor * tol` is
    therefore thrown away and retried at half the size, so the bound holds on
    what was actually integrated rather than on what was predicted.

    `dt_max` matters as much as `tol`. At 700 km the tolerance alone permits
    steps of years, which would step straight over the entire problem; the cap
    is what keeps the run resolving the atmosphere it is flying through.
    Callers propagating with time-varying space weather must set `dt_max`
    below the timescale on which their density model changes.

    Returns a Trajectory whose `stats` records the steps actually taken. Check
    `stats.tolerance_respected` -- if steps hit `dt_min`, the tolerance was not
    met and the result is not trustworthy at face value.
    """
    if tol <= 0:
        raise ValueError("tol must be positive")
    if dt_min <= 0 or dt_max < dt_min:
        raise ValueError(f"need 0 < dt_min <= dt_max, got {dt_min}, {dt_max}")

    t = float(t0)
    y = np.asarray(y0, dtype=float).copy()

    ts = [t]
    a_hist = [y[0]]
    m_hist = [y[1]]

    outcome = Outcome.INDETERMINATE
    outcome_time: float | None = None
    exhausted = False
    stats = StepStats()

    while t < t_max:
        f = deriv(t, y)
        rate = abs(float(f[0]))

        # Guard the rate->step division. A vanishing da/dt (drag exactly
        # balancing thrust, PHYSICS.md §4) means the state is not moving and
        # the step is limited only by dt_max.
        dt = dt_max if rate <= 0.0 else tol * y[0] / rate
        capped_by_max = dt >= dt_max
        dt = min(max(dt, dt_min), dt_max, t_max - t)

        while True:
            # A *trial* step is allowed to fail. RK4's interior stages evaluate
            # `deriv` at states the step never actually visits, and during a
            # steep decay those can be far below the surface, where a density
            # model may legitimately overflow or return nonsense. That is a
            # signal the step is too long, not a bug -- reject and halve. Only
            # a failure at dt_min is unrecoverable.
            try:
                y_new = rk4_step(deriv, t, y, dt)
            except (OverflowError, FloatingPointError, ZeroDivisionError):
                y_new = None

            usable = (
                y_new is not None
                and np.all(np.isfinite(y_new))
                and y_new[0] > 0.0
                and y_new[1] > 0.0
            )
            frac = abs(y_new[0] - y[0]) / y[0] if usable else float("inf")

            if usable and frac <= reject_factor * tol:
                break
            if dt <= dt_min:
                if not usable:
                    raise FloatingPointError(
                        f"step at t={t} produced an unusable state even at "
                        f"dt_min={dt_min}: y={y_new!r}"
                    )
                break  # tolerance not met, but recorded via n_at_dt_min below
            stats.n_rejected += 1
            dt = max(dt * 0.5, dt_min)
            capped_by_max = False

        t_new = t + dt

        stats.n_accepted += 1
        stats.dt_min_used = min(stats.dt_min_used, dt)
        stats.dt_max_used = max(stats.dt_max_used, dt)
        stats.max_frac_change = max(stats.max_frac_change, frac)
        if capped_by_max:
            stats.n_at_dt_max += 1
        if dt <= dt_min:
            stats.n_at_dt_min += 1

        if dry_mass_kg is not None and y_new[1] <= dry_mass_kg:
            exhausted = True

        h_prev = y[0] - R_E
        h = y_new[0] - R_E
        terminated = False
        if h < reentry_altitude_m:
            outcome = Outcome.REENTERED
            outcome_time = _crossing_time(t, h_prev, t_new, h, reentry_altitude_m)
            terminated = True
        elif shell_altitude_m is not None and h >= shell_altitude_m:
            outcome = Outcome.REACHED_SHELL
            outcome_time = _crossing_time(t, h_prev, t_new, h, shell_altitude_m)
            terminated = True

        t, y = t_new, y_new
        ts.append(t)
        a_hist.append(y[0])
        m_hist.append(y[1])

        if terminated:
            break
    else:
        if exhausted:
            outcome = Outcome.PROPELLANT_EXHAUSTED

    return Trajectory(
        np.array(ts), np.array(a_hist), np.array(m_hist), outcome, outcome_time,
        stats=stats,
    )
