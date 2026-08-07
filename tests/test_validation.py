"""Validation tests 1-3 from PHYSICS.md §8.

These are the analytic tests: they need no atmosphere model and no data files,
so they are pure unit tests (ARCHITECTURE.md §3). Test 4 -- the Baruah et al. (2024)
reproduction -- is Phase 2 and lives in `sim/validate.py`.

Each test prints its measured error so the numbers can be lifted straight into
the writeup (ARCHITECTURE.md §5, Phase 6).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sim.constants import MU, R_E
from sim.critical import critical_altitude, critical_density
from sim.dynamics import da_dt, derivatives
from sim.integrator import propagate, rk4_step

A0 = R_E + 210e3          # insertion radius for the Feb 2022 case, m
MASS = 227.0              # kg, Baruah et al. 2024
THRUST = 0.071            # N, DERIVED (satellite_specs.json, v1_5.thrust_n)
WEEK = 7 * 86400.0
DAY = 86400.0


def test_1_energy_conservation(capsys):
    """PHYSICS.md §8, Test 1 -- energy conservation.

    rho = 0 and F = 0 => da/dt = 0 exactly. After a simulated week the
    semi-major axis must be unchanged to within floating point error
    (relative change < 1e-12).

    Catches: sign errors, spurious terms, integrator bugs.
    """
    def deriv(t, y):
        return derivatives(y, rho=0.0, thrust=0.0, cd=2.2, area=1.0, isp=None)

    traj = propagate(
        deriv, np.array([A0, MASS]), dt=10.0, t_max=WEEK, sample_every=8640
    )

    rel = abs(traj.a_m[-1] - A0) / A0
    rel_mass = abs(traj.mass_kg[-1] - MASS) / MASS
    with capsys.disabled():
        print(f"\n  Test 1: da/dt=0 over 1 week, relative change in a = {rel:.3e}")

    assert rel < 1e-12, f"semi-major axis drifted by {rel:.3e} relative"
    assert rel_mass == 0.0, "mass changed with the thruster off"


def test_2_pure_thrust_spiral(capsys):
    """PHYSICS.md §8, Test 2 -- pure thrust spiral against the closed form.

        a(t) = a0 / (1 - (F/m)*t*sqrt(a0/MU))**2

    rho = 0, F nonzero, m held constant. Must agree to within 0.01% over a
    one-day propagation.

    Catches: incorrect exponents, wrong constant factors.
    """
    def deriv(t, y):
        # isp=None holds mass constant, which is what the closed form assumes.
        return derivatives(y, rho=0.0, thrust=THRUST, cd=2.2, area=1.0, isp=None)

    traj = propagate(
        deriv, np.array([A0, MASS]), dt=10.0, t_max=DAY, sample_every=8640
    )

    t_end = traj.t_s[-1]
    analytic = A0 / (1.0 - (THRUST / MASS) * t_end * math.sqrt(A0 / MU)) ** 2
    numeric = traj.a_m[-1]
    rel = abs(numeric - analytic) / analytic

    with capsys.disabled():
        print(
            f"  Test 2: 1 day thrust spiral, "
            f"numeric {numeric / 1e3:.6f} km vs analytic {analytic / 1e3:.6f} km, "
            f"relative error = {rel:.3e} ({rel * 100:.2e}%)"
        )
        print(f"          altitude gain over the day = "
              f"{(numeric - A0) / 1e3:.3f} km")

    assert rel < 1e-4, f"thrust spiral off by {rel * 100:.4f}%, limit 0.01%"
    assert traj.mass_kg[-1] == MASS, "mass was not held constant"


def test_3_critical_density_is_a_fixed_point(capsys):
    """PHYSICS.md §8, Test 3 -- critical density.

    Place a satellite at exactly h_crit computed from §4, run for one hour with
    constant density, and da/dt must be zero to within 1e-9 relative.

    Catches: inconsistency between the decay equation and the critical-altitude
    solver. The atmosphere here is a deliberately crude analytic exponential --
    the test is about self-consistency between `critical.py` and `dynamics.py`,
    not about NRLMSIS.
    """
    cd, area = 2.2, 4.48

    # Exponential atmosphere, roughly NRLMSIS-like near 200 km. Any smooth
    # monotonically decreasing rho(h) exercises the solver identically.
    def rho_of_h(h_m: float) -> float:
        return 1.468e-10 * math.exp(-(h_m - 210e3) / 37.0e3)

    h_crit = critical_altitude(rho_of_h, THRUST, cd, area)
    a_crit = R_E + h_crit
    rho_crit = critical_density(a_crit, THRUST, cd, area)

    # The solver's root really is a root of the decay equation.
    rate = da_dt(a_crit, MASS, rho_crit, THRUST, cd, area)
    thrust_term = 2.0 * (THRUST / MASS) * a_crit**1.5 / math.sqrt(MU)
    rel_rate = abs(rate) / thrust_term

    # And it is a fixed point of the integrator over an hour at constant rho.
    def deriv(t, y):
        return derivatives(y, rho=rho_crit, thrust=THRUST, cd=cd, area=area, isp=None)

    traj = propagate(deriv, np.array([a_crit, MASS]), dt=10.0, t_max=3600.0)
    rel_drift = abs(traj.a_m[-1] - a_crit) / a_crit

    with capsys.disabled():
        print(
            f"  Test 3: h_crit = {h_crit / 1e3:.4f} km "
            f"(Cd={cd}, A={area} m2, F={THRUST} N)\n"
            f"          |da/dt| / thrust term   = {rel_rate:.3e}\n"
            f"          relative drift in a over 1 h = {rel_drift:.3e}"
        )

    assert rel_rate < 1e-9, f"da/dt at h_crit is {rel_rate:.3e} relative, limit 1e-9"
    assert rel_drift < 1e-9, f"a drifted {rel_drift:.3e} relative over an hour"


def test_rk4_step_matches_known_ode():
    """RK4 itself, checked against dy/dt = y on a problem with a known answer.

    Guards the integrator independently of the orbital dynamics: a fourth-order
    method must show ~16x error reduction when the step is halved.
    """
    def f(t, y):
        return y

    def final_error(dt):
        y = np.array([1.0])
        t, n = 0.0, int(round(1.0 / dt))
        for i in range(n):
            y = rk4_step(f, t, y, dt)
            t = (i + 1) * dt
        return abs(y[0] - math.e) / math.e

    e_coarse, e_fine = final_error(0.1), final_error(0.05)
    assert e_coarse < 1e-6
    assert 10.0 < e_coarse / e_fine < 25.0, "convergence is not fourth order"


def test_step_halving_convergence(capsys):
    """PHYSICS.md §7 -- verify the 10 s step by halving it.

    Run at 10 s and 5 s on a decaying case with drag active; the final altitude
    must differ by less than 0.1%.
    """
    cd, area, rho = 1.0, 4.48, 1.468e-10

    def deriv(t, y):
        return derivatives(y, rho=rho, thrust=0.0, cd=cd, area=area, isp=None)

    results = {}
    for dt in (10.0, 5.0):
        traj = propagate(deriv, np.array([A0, MASS]), dt=dt, t_max=DAY)
        results[dt] = traj.a_m[-1] - R_E

    rel = abs(results[10.0] - results[5.0]) / results[5.0]
    with capsys.disabled():
        print(
            f"  Step halving: h(10 s) = {results[10.0] / 1e3:.6f} km, "
            f"h(5 s) = {results[5.0] / 1e3:.6f} km, relative diff = {rel:.3e}"
        )
    assert rel < 1e-3, f"step-size sensitivity {rel:.3e} exceeds 0.1%"


def test_drag_lowers_and_thrust_raises():
    """Sign convention guard for PHYSICS.md §3.2."""
    assert da_dt(A0, MASS, 1.468e-10, 0.0, 2.2, 1.0) < 0.0
    assert da_dt(A0, MASS, 0.0, THRUST, 2.2, 1.0) > 0.0


def test_critical_density_undefined_without_thrust():
    """In safe mode (F = 0) there is no balance point -- drag always wins."""
    with pytest.raises(ValueError, match="undefined for F <= 0"):
        critical_density(A0, 0.0, 1.0, 1.0)
