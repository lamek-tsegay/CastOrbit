"""Adaptive-step RK4 driver. V2_BRIEF.md §4.

The adaptive driver is only allowed to exist if it gives the same answers as
the fixed-step driver that was validated in V1. These tests hold it to that,
analytically where a closed form exists and against `propagate` where one
does not.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sim.constants import MU, R_E
from sim.dynamics import derivatives
from sim.integrator import propagate, propagate_adaptive
from sim.satellite import Outcome

DAY = 86400.0
A0 = R_E + 300e3
MASS = 260.0
THRUST = 0.071


def exp_atmosphere(scale_m: float, rho_ref: float = 2.5e-11, h_ref: float = 300e3):
    """Exponential stand-in for NRLMSIS, for tests that must stay fast.

    `rho_ref` is NRLMSIS's rough quiet-time value at 300 km, so the resulting
    decay timescales are physically plausible rather than arbitrary.

    Clamped at a 100 km floor exactly as `DensityGrid.lookup` clamps, and for
    the same reason: RK4's interior stages can evaluate below the model's
    valid range, and unbounded exponential growth there drives the state
    negative inside a single step. Production density lookups already clamp,
    so an unclamped test atmosphere would be testing a hazard the real
    simulator does not have.
    """
    def rho_of_h(h_m: float) -> float:
        h = max(float(h_m), 100e3)
        return rho_ref * math.exp(-(h - h_ref) / scale_m)
    return rho_of_h


def test_zero_rate_takes_maximum_steps():
    """rho = 0, F = 0 => da/dt = 0 exactly, so nothing bounds the step but dt_max.

    PHYSICS.md §8 Test 1 as seen by the step controller. Also guards the
    divide-by-rate: a vanishing derivative must not produce a NaN or an
    infinite step.
    """
    def deriv(t, y):
        return derivatives(y, rho=0.0, thrust=0.0, cd=2.2, area=1.0, isp=None)

    traj = propagate_adaptive(
        deriv, np.array([A0, MASS]), t_max=7 * DAY, tol=1e-4, dt_max=DAY
    )

    rel = abs(traj.a_m[-1] - A0) / A0
    assert rel < 1e-12, f"a drifted by {rel:.3e} with no forces acting"
    assert traj.stats.n_accepted == 7, "should take exactly 7 one-day steps"
    assert traj.stats.n_at_dt_max == 7
    assert traj.stats.tolerance_respected


def test_thrust_spiral_matches_closed_form():
    """PHYSICS.md §8 Test 2, run adaptively instead of at a fixed 10 s.

        a(t) = a0 / (1 - (F/m)*t*sqrt(a0/MU))**2

    Same 0.01% bar the fixed-step version is held to.
    """
    def deriv(t, y):
        return derivatives(y, rho=0.0, thrust=THRUST, cd=2.2, area=1.0, isp=None)

    traj = propagate_adaptive(
        deriv, np.array([A0, MASS]), t_max=DAY, tol=1e-4, dt_max=DAY
    )

    t_end = traj.t_s[-1]
    analytic = A0 / (1.0 - (THRUST / MASS) * t_end * math.sqrt(A0 / MU)) ** 2
    rel = abs(traj.a_m[-1] - analytic) / analytic

    assert rel < 1e-4, f"adaptive thrust spiral off by {rel * 100:.4f}%"
    assert traj.mass_kg[-1] == MASS, "mass was not held constant"


def test_step_scales_with_the_physics():
    """The whole point of §4: weak drag buys long steps, strong drag forces short ones.

    Uses an exponential-atmosphere stand-in rather than NRLMSIS so the test is
    fast and depends on nothing but the controller.
    """
    def make_deriv(h_scale_m: float):
        rho_of_h = exp_atmosphere(h_scale_m)
        def deriv(t, y):
            return derivatives(y, rho_of_h(y[0] - R_E),
                               thrust=0.0, cd=2.2, area=4.0, isp=None)
        return deriv

    low = propagate_adaptive(
        make_deriv(50e3), np.array([R_E + 250e3, MASS]),
        t_max=DAY, tol=1e-4, dt_max=DAY,
    )
    high = propagate_adaptive(
        make_deriv(50e3), np.array([R_E + 700e3, MASS]),
        t_max=DAY, tol=1e-4, dt_max=DAY,
    )

    assert high.stats.n_accepted < low.stats.n_accepted, (
        "a satellite in thinner air should need fewer steps, got "
        f"{high.stats.n_accepted} at 700 km vs {low.stats.n_accepted} at 250 km"
    )


def test_reentry_time_is_interpolated_within_the_step():
    """Steps can be hours long, so the crossing must be resolved inside one.

    Without interpolation the reported reentry time would be the end of the
    step that overshot, which for a 1 h step is a 1 h error.
    """
    rho_of_h = exp_atmosphere(40e3)

    def deriv(t, y):
        return derivatives(y, rho_of_h(y[0] - R_E),
                           thrust=0.0, cd=2.2, area=8.0, isp=None)

    traj = propagate_adaptive(
        deriv, np.array([R_E + 150e3, MASS]), t_max=30 * DAY, tol=1e-4, dt_max=DAY
    )

    assert traj.outcome is Outcome.REENTERED
    # The crossing must lie strictly inside the final step, not at its end.
    assert traj.t_s[-2] < traj.outcome_time_s <= traj.t_s[-1]
    h_at_crossing = np.interp(traj.outcome_time_s, traj.t_s, traj.a_m - R_E)
    assert abs(h_at_crossing - 100e3) < 200.0, (
        f"interpolated crossing sits at {h_at_crossing / 1e3:.3f} km, not 100 km"
    )


def test_rejection_path_catches_a_rate_that_jumps_mid_step():
    """Directly exercise reject-and-halve with a derivative that spikes.

    Worth being explicit about why this test is contrived: in *physical* decay
    the rejection path almost never fires, and that is a property of the
    controller rather than luck. Bounding |da/a| per step also bounds how far
    the satellite moves down the density profile within the step, so rho --
    and therefore the rate -- cannot change much before the step ends. Runs
    against both the exponential stand-in and real NRLMSIS come out at
    max_frac_change ~1.03e-4 against a 1e-4 tolerance, with zero rejections.

    So the safety net is real but rarely load-bearing, which means the only
    way to know it works is to trip it deliberately. This derivative jumps
    100x at a fixed time, something no atmosphere does.
    """
    # Rates chosen so the predicted step lands at ~2000 s, putting the jump
    # squarely inside a step rather than on a boundary.
    t_jump = 5000.0

    def deriv(t, y):
        rate = -0.34 if t < t_jump else -34.0
        return np.array([rate, 0.0])

    tol, reject_factor = 1e-4, 2.0
    traj = propagate_adaptive(
        deriv, np.array([R_E + 400e3, MASS]), t_max=8000.0,
        tol=tol, dt_min=1e-3, dt_max=DAY, reject_factor=reject_factor,
    )

    assert traj.stats.n_rejected > 0, "the rate jump should have forced retries"
    assert traj.stats.max_frac_change <= reject_factor * tol * (1 + 1e-9), (
        f"a step of |da/a| = {traj.stats.max_frac_change:.3e} was accepted, "
        f"above the {reject_factor * tol:.3e} bound"
    )


def test_physical_decay_holds_the_bound_without_needing_rejection():
    """The companion claim: on a real decay profile the prediction is enough.

    Documents the behaviour the test above works around -- if this ever starts
    needing rejections, the controller's self-limiting property has broken.
    """
    rho_of_h = exp_atmosphere(12e3, rho_ref=3.0e-10, h_ref=200e3)

    def deriv(t, y):
        return derivatives(y, rho_of_h(y[0] - R_E),
                           thrust=0.0, cd=2.2, area=8.0, isp=None)

    tol = 1e-4
    traj = propagate_adaptive(
        deriv, np.array([R_E + 200e3, MASS]), t_max=30 * DAY,
        tol=tol, dt_min=1e-3, dt_max=DAY,
    )

    assert traj.outcome is Outcome.REENTERED
    assert traj.stats.tolerance_respected
    assert traj.stats.n_rejected == 0
    assert traj.stats.max_frac_change < 1.1 * tol, (
        f"start-of-step prediction overshot to {traj.stats.max_frac_change:.3e}"
    )


@pytest.mark.parametrize("area_m2", [1.0, 4.0])
def test_adaptive_matches_fixed_step(area_m2):
    """The regression that matters: same answer as V1's validated driver.

    Held to 0.1%, the V2_BRIEF.md §4 gate. Exponential atmosphere again, so
    this runs in milliseconds and stays a unit test; the real Baruah
    reproduction against NRLMSIS is in `test_baruah.py`.
    """
    rho_of_h = exp_atmosphere(45e3)

    def deriv(t, y):
        return derivatives(y, rho_of_h(y[0] - R_E),
                           thrust=0.0, cd=2.2, area=area_m2, isp=None)

    y0 = np.array([R_E + 250e3, MASS])
    fixed = propagate(deriv, y0, dt=10.0, t_max=20 * DAY)
    adaptive = propagate_adaptive(deriv, y0, t_max=20 * DAY, tol=1e-4, dt_max=DAY)

    assert adaptive.outcome is fixed.outcome
    rel_a = abs(adaptive.a_m[-1] - fixed.a_m[-1]) / fixed.a_m[-1]
    assert rel_a < 1e-3, f"final a differs by {rel_a * 100:.4f}%"

    assert adaptive.stats.n_accepted < fixed.t_s.size / 50, (
        "adaptive stepping bought less than a 50x reduction in steps; "
        f"{adaptive.stats.n_accepted} vs {fixed.t_s.size}"
    )


def test_tightening_tolerance_does_not_move_the_answer():
    """The rate bound is not the accuracy limit -- RK4 truncation is far below it.

    This is the claim `propagate_adaptive`'s docstring makes, and the reason a
    rate-based controller is defensible instead of a truncation-error
    estimator. If halving tol moved the answer, the controller would be
    setting the accuracy and the docstring would be wrong.
    """
    rho_of_h = exp_atmosphere(45e3)

    def deriv(t, y):
        return derivatives(y, rho_of_h(y[0] - R_E),
                           thrust=0.0, cd=2.2, area=4.0, isp=None)

    y0 = np.array([R_E + 250e3, MASS])
    coarse = propagate_adaptive(deriv, y0, t_max=20 * DAY, tol=1e-4, dt_max=DAY)
    fine = propagate_adaptive(deriv, y0, t_max=20 * DAY, tol=1e-6, dt_max=DAY)

    rel = abs(coarse.a_m[-1] - fine.a_m[-1]) / fine.a_m[-1]
    assert rel < 1e-4, (
        f"10x tighter tolerance moved the answer by {rel:.3e}; the step "
        "controller, not RK4, is setting the accuracy"
    )


def test_dt_min_floor_is_reported_not_hidden():
    """Hitting the floor means the tolerance was not met. That must be visible.

    A silent accuracy failure is the specific way an adaptive integrator
    betrays you, so `tolerance_respected` exists to be asserted on.
    """
    # An 8 km scale height is far steeper than the real thermosphere, which is
    # the point: it drives the runaway hard enough that the controller wants
    # sub-second steps and the 1 s floor stops it getting them.
    rho_of_h = exp_atmosphere(8e3, rho_ref=3.0e-10, h_ref=200e3)

    def deriv(t, y):
        return derivatives(y, rho_of_h(y[0] - R_E),
                           thrust=0.0, cd=2.2, area=8.0, isp=None)

    traj = propagate_adaptive(
        deriv, np.array([R_E + 200e3, MASS]), t_max=30 * DAY,
        tol=1e-4, dt_min=1.0, dt_max=DAY,
    )

    assert traj.stats.n_at_dt_min > 0
    assert not traj.stats.tolerance_respected
    assert traj.stats.as_dict()["tolerance_respected"] is False
    # And the floored steps really did blow the bound -- this is what the flag
    # is warning about, not a cosmetic edge case.
    assert traj.stats.max_frac_change > 1e-3


def test_rejects_invalid_step_bounds():
    def deriv(t, y):
        return derivatives(y, rho=0.0, thrust=0.0, cd=2.2, area=1.0, isp=None)

    y0 = np.array([A0, MASS])
    with pytest.raises(ValueError, match="tol must be positive"):
        propagate_adaptive(deriv, y0, t_max=DAY, tol=0.0)
    with pytest.raises(ValueError, match="dt_min <= dt_max"):
        propagate_adaptive(deriv, y0, t_max=DAY, dt_min=100.0, dt_max=1.0)
