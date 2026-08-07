"""Hand-written RK4 integrator and the propagation loop.

PHYSICS.md §7. Fixed step, classical fourth-order Runge-Kutta, written out
directly. `scipy.integrate.solve_ivp` is explicitly forbidden for the main
propagation loop (ARCHITECTURE.md §4).
"""

from __future__ import annotations

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
    ):
        self.t_s = t_s
        self.a_m = a_m
        self.mass_kg = mass_kg
        self.outcome = outcome
        self.outcome_time_s = outcome_time_s

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
