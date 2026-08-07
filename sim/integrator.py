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
