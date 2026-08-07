"""The equations of motion.

PHYSICS.md §3. Two coupled ODEs, state vector [a, m]. This is the entire
simulator; everything else is bookkeeping around these three functions.
"""

from __future__ import annotations

import math

import numpy as np

from .constants import G0, MU

_SQRT_MU = math.sqrt(MU)


def da_dt(a: float, m: float, rho: float, thrust: float,
          cd: float, area: float) -> float:
    """Rate of change of semi-major axis, m/s.

    PHYSICS.md §3.2:

        da/dt = 2*(F/m)*a**1.5/sqrt(MU) - rho*(Cd*A/m)*sqrt(MU*a)

    Args:
        a: semi-major axis (== orbital radius, circular), m
        m: satellite mass, kg
        rho: local atmospheric mass density, kg/m^3
        thrust: thruster force F, N (0 in safe mode, PHYSICS.md §5)
        cd: drag coefficient
        area: ram cross-sectional area A, m^2
    """
    thrust_term = 2.0 * (thrust / m) * a**1.5 / _SQRT_MU
    drag_term = rho * (cd * area / m) * math.sqrt(MU * a)
    return thrust_term - drag_term


def dm_dt(thrust: float, isp: float | None) -> float:
    """Propellant mass depletion rate, kg/s.

    PHYSICS.md §3.2:  dm/dt = -F / (Isp * G0)

    `isp=None` (or non-finite) means mass is held constant. That is not a
    physical case; it exists because validation Test 2 (PHYSICS.md §8) requires
    a fixed-mass thrust spiral to compare against the closed-form solution.

    """
    if thrust == 0.0:
        return 0.0
    if isp is None or not math.isfinite(isp):
        return 0.0
    return -thrust / (isp * G0)


def derivatives(
    state: np.ndarray,
    rho: float,
    thrust: float,
    cd: float,
    area: float,
    isp: float | None = None,
) -> np.ndarray:
    """Full state derivative d[a, m]/dt for the RK4 integrator.

    PHYSICS.md §3.2.
    """
    a, m = float(state[0]), float(state[1])
    return np.array([da_dt(a, m, rho, thrust, cd, area), dm_dt(thrust, isp)])
