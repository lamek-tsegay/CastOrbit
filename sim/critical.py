"""Critical density and critical altitude -- the knife edge.

PHYSICS.md §4. Setting da/dt = 0 in §3.2 and solving for density gives

    rho_crit = 2 * F * a / (Cd * A * MU)

Below `rho_crit` thrust wins and the satellite climbs; above it drag wins and
it sinks. Because density rises steeply as altitude falls, the process runs
away in both directions.
"""

from __future__ import annotations

from typing import Callable

from .constants import MU, R_E


def critical_density(a: float, thrust: float, cd: float, area: float) -> float:
    """Density at which thrust exactly balances drag, kg/m^3.

    PHYSICS.md §4:  rho_crit = 2*F*a / (Cd*A*MU)

    Args:
        a: semi-major axis, m
        thrust: F, N (this is the *nominal* thrust; in safe mode F = 0 and
            there is no critical density -- drag always wins)
        cd: drag coefficient
        area: ram area A, m^2
    """
    if thrust <= 0.0:
        raise ValueError(
            "critical density is undefined for F <= 0: with no thrust the "
            "satellite always decays (PHYSICS.md §5)"
        )
    return 2.0 * thrust * a / (cd * area * MU)


def critical_altitude(
    rho_of_h: Callable[[float], float],
    thrust: float,
    cd: float,
    area: float,
    h_lo: float = 100e3,
    h_hi: float = 600e3,
    tol_m: float = 1e-3,
    max_iter: int = 200,
) -> float:
    """Altitude, m, at which the local density equals the critical density.

    PHYSICS.md §4. Bisection over h in [100 km, 600 km]; density is
    monotonically decreasing in altitude so bisection is safe.

    Note that `rho_crit` is evaluated self-consistently at a = R_E + h rather
    than held at some fixed reference radius: h_crit is the altitude at which a
    satellite *located there* has da/dt = 0. The a-dependence of rho_crit is
    weak (linear) next to the near-exponential altitude dependence of rho, so
    the root is still bracketed cleanly.

    Args:
        rho_of_h: callable mapping altitude in metres to density in kg/m^3.

    Raises:
        ValueError: if the root is not bracketed on [h_lo, h_hi].
    """

    def excess(h: float) -> float:
        """rho(h) - rho_crit(h). Positive => drag wins, negative => thrust wins."""
        return rho_of_h(h) - critical_density(R_E + h, thrust, cd, area)

    f_lo, f_hi = excess(h_lo), excess(h_hi)
    if f_lo <= 0.0:
        raise ValueError(
            f"no critical altitude above {h_lo / 1e3:.0f} km: thrust already wins "
            f"there (rho - rho_crit = {f_lo:.3e}). Configuration is drag-safe."
        )
    if f_hi >= 0.0:
        raise ValueError(
            f"no critical altitude below {h_hi / 1e3:.0f} km: drag still wins "
            f"there (rho - rho_crit = {f_hi:.3e})."
        )

    for _ in range(max_iter):
        h_mid = 0.5 * (h_lo + h_hi)
        if h_hi - h_lo < tol_m:
            return h_mid
        if excess(h_mid) > 0.0:
            h_lo = h_mid
        else:
            h_hi = h_mid
    return 0.5 * (h_lo + h_hi)
