"""Bounding Stage 1: mass intervals and their propagation.

Gate 9 failed at 1/3. Rather than build stage 2, the mass model now refuses to
produce a point estimate it cannot support, and this module carries that
refusal through to the verdict. These tests exist to make sure the refusal
cannot be lost along the way -- a silently-dropped caveat is worse than no
caveat, because it looks like an answer.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sim.atmosphere import ClimatologyDensity, SpaceWeather
from sim.bounded import (
    Interval,
    assess_compliance_bounded,
    ballistic_coefficient_interval,
    decay_time_interval,
    natural_decay_years_interval,
)
from sim.disposal import Compliance, DECIDED_VERDICTS
from sim.mass_model import (
    MAX_BRACKET_POWER_RATIO,
    MAX_SCATTER_RATIO,
    MassEstimate,
    estimate_dry_mass,
    score_held_out,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"


@pytest.fixture(scope="module")
def atmos():
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return ClimatologyDensity.for_level(
        "mean", SpaceWeather.load(DATA), lat_deg=53.0
    )


@pytest.fixture(scope="module")
def resolvable():
    """1700 W earth observation -- a dense, homogeneous part of the table."""
    return estimate_dry_mass(1700.0, "earth_observation")


@pytest.fixture(scope="module")
def unresolvable():
    """320 W earth observation -- PROBA-V's power, where Gate 9 failed."""
    return estimate_dry_mass(320.0, "earth_observation")


# --------------------------------------------------------------------------
# The interval itself
# --------------------------------------------------------------------------

def test_interval_rejects_inversion():
    with pytest.raises(ValueError, match="inverted interval"):
        Interval(10.0, 1.0)


def test_every_estimate_carries_an_interval(resolvable, unresolvable):
    for est in (resolvable, unresolvable):
        lo, hi = est.interval_kg
        assert lo > 0 and hi >= lo


def test_resolvable_estimate_keeps_its_point_inside_its_interval(resolvable):
    lo, hi = resolvable.interval_kg
    assert resolvable.resolvable
    assert lo <= resolvable.mass_kg <= hi


def test_unresolvable_estimate_has_no_point_at_all(unresolvable):
    """Not a hedged number, not a midpoint. None."""
    assert not unresolvable.resolvable
    assert unresolvable.mass_kg is None
    assert unresolvable.refusal_reasons
    assert "cannot resolve" in unresolvable.summary


def test_the_dataclass_refuses_to_hold_a_contradiction():
    """A refusal carrying a point estimate would defeat the whole mechanism."""
    with pytest.raises(ValueError, match="must not carry a point mass"):
        MassEstimate(
            mass_kg=100.0, interval_kg=(50.0, 200.0), resolvable=False,
            tag="estimated", method="x",
        )
    with pytest.raises(ValueError, match="must carry a point mass"):
        MassEstimate(
            mass_kg=None, interval_kg=(50.0, 200.0), resolvable=True,
            tag="estimated", method="x",
        )


def test_wide_bracket_refuses(capsys):
    """PROBA-V's 42x bracket is the case that motivated the threshold."""
    est = estimate_dry_mass(320.0, "earth_observation")
    assert est.bracket_ratio > MAX_BRACKET_POWER_RATIO
    assert not est.resolvable
    assert any("apart in power" in r for r in est.refusal_reasons)
    with capsys.disabled():
        print(f"\n  320 W: bracket {est.bracket_ratio:.0f}x -> refused, "
              f"range {est.interval_kg[0]:.0f}-{est.interval_kg[1]:.0f} kg")


def test_high_scatter_refuses(unresolvable):
    assert unresolvable.scatter_ratio > MAX_SCATTER_RATIO
    assert any("vary by" in r for r in unresolvable.refusal_reasons)


def test_thin_local_sample_falls_back_to_class_scatter(capsys):
    """The safeguard against false confidence, which is the real failure mode.

    GOES-16's local window holds two spacecraft whose kg/W agree to 1.11x. Used
    naively that would produce a tight interval excluding the true mass by 50%.
    Falling back to the class-wide spread is what prevents that.
    """
    est = estimate_dry_mass(4000.0, "earth_observation")
    assert est.scatter_basis == "class"
    assert est.n_samples >= 4
    assert any("below the" in w and "local scatter" in w for w in est.warnings)
    with capsys.disabled():
        print(f"\n  4000 W: scatter basis '{est.scatter_basis}' "
              f"({est.scatter_ratio:.2f}x), range "
              f"{est.interval_kg[0]:.0f}-{est.interval_kg[1]:.0f} kg")


def test_dense_homogeneous_region_still_resolves(resolvable):
    """Bounding must not make the method useless where it actually works."""
    assert resolvable.resolvable
    assert resolvable.scatter_basis == "local"
    assert resolvable.mass_kg is not None


# --------------------------------------------------------------------------
# Gate 9, re-scored under the bound
# --------------------------------------------------------------------------

def test_gate_9_bounded_every_interval_contains_the_truth(capsys):
    """Gate 9 stays FAILED. What changed is that the misses now decline to answer.

    The bar for a refusal is different from the bar for an estimate: a withheld
    point estimate is only honest if the range offered instead actually
    contains the true mass. All three do.
    """
    rows = score_held_out()
    with capsys.disabled():
        print("\n  GATE 9 (bounded) -- 25% bar unchanged, table unchanged")
        for r in rows:
            lo, hi = r["interval_kg"]
            verdict = (
                f"{r['predicted_kg']:6.0f} kg ({r['error_frac'] * 100:+.1f}%)"
                if r["resolvable"] else "REFUSED".rjust(18)
            )
            print(f"    {r['name'][:24]:24s} actual {r['actual_kg']:6.0f} kg  "
                  f"{verdict}  range {lo:5.0f}-{hi:5.0f} kg  "
                  f"contains: {r['interval_contains_actual']}")

    for r in rows:
        assert r["interval_contains_actual"], (
            f"{r['name']}: interval {r['interval_kg']} excludes the true mass "
            f"{r['actual_kg']} kg -- a refusal that returns a wrong range is "
            "worse than a wrong point estimate"
        )

    # The two that missed the 25% bar are exactly the two now refused.
    refused = {r["id"] for r in rows if not r["resolvable"]}
    assert refused == {"proba_v", "goes_16"}
    resolved = [r for r in rows if r["resolvable"]]
    assert len(resolved) == 1 and resolved[0]["id"] == "sentinel_2"
    assert resolved[0]["within_25pct"], (
        "the one case still answered must be inside the bar"
    )


# --------------------------------------------------------------------------
# Propagation
# --------------------------------------------------------------------------

def test_ballistic_coefficient_is_monotone_in_mass(resolvable):
    bc = ballistic_coefficient_interval(Interval(*resolvable.interval_kg), 2.2, 1.0)
    lo, hi = resolvable.interval_kg
    assert bc.lo == pytest.approx(lo / 2.2)
    assert bc.hi == pytest.approx(hi / 2.2)
    assert bc.lo < bc.hi


def test_ballistic_coefficient_rejects_bad_geometry():
    with pytest.raises(ValueError, match="must be positive"):
        ballistic_coefficient_interval(Interval(100, 200), 0.0, 1.0)


def test_decay_time_interval_orders_light_before_heavy(atmos):
    """Lighter decays faster. The physical claim behind endpoint evaluation."""
    interval, detail = decay_time_interval(
        atmos, 400e3, Interval(500.0, 1500.0), cd=2.2, area_m2=3.0,
        t_max_s=200 * 365.25 * 86400.0,
    )
    assert interval is not None
    assert detail["t_at_low_mass_s"] < detail["t_at_high_mass_s"]
    assert interval.lo == detail["t_at_low_mass_s"]


def test_decay_time_returns_none_rather_than_the_horizon(atmos):
    """'Never reenters' must not be reported as 'reenters at exactly t_max'."""
    interval, detail = decay_time_interval(
        atmos, 900e3, Interval(500.0, 1500.0), cd=2.2, area_m2=1.0,
        t_max_s=5 * 365.25 * 86400.0,
    )
    assert interval is None
    assert "unbounded above" in detail["note"]


def test_decay_years_interval_inherits_the_mass_bound(atmos, resolvable):
    years, detail = natural_decay_years_interval(
        atmos, 450e3, resolvable, cd=2.2, area_m2=3.0
    )
    assert years is not None
    assert years.lo < years.hi
    assert detail["mass_resolvable"] is True


# --------------------------------------------------------------------------
# The verdict must not render on an unbounded mass
# --------------------------------------------------------------------------

def test_unresolvable_mass_produces_no_verdict(atmos, unresolvable):
    """**The requirement.** No verdict, and renderable is False."""
    r = assess_compliance_bounded(
        atmos, 705e3, unresolvable, cd=2.2, area_m2=1.0, isp_s=1666.0,
        propellant_available_kg=5.0, solar_activity="mean",
    )
    assert r.verdict is Compliance.NOT_ASSESSABLE
    assert r.renderable is False
    assert r.optimistic is None and r.pessimistic is None
    assert r.delta_v_interval_ms is None
    assert r.propellant_interval_kg is None
    # The reason the mass model gave must survive to the verdict.
    assert any("could not be resolved" in n for n in r.notes)
    for reason in unresolvable.refusal_reasons:
        assert reason in r.notes


def test_unresolvable_verdict_serialises_as_non_renderable(atmos, unresolvable):
    """A UI reads one boolean. It must be present and false in the payload."""
    d = assess_compliance_bounded(
        atmos, 705e3, unresolvable, cd=2.2, area_m2=1.0, isp_s=1666.0,
        propellant_available_kg=5.0, solar_activity="mean",
    ).as_dict()
    assert d["renderable"] is False
    assert d["verdict"] == "NOT_ASSESSABLE"
    assert d["mass"]["resolvable"] is False
    assert d["mass"]["mass_kg"] is None


def test_resolvable_mass_produces_a_verdict_with_intervals(atmos, resolvable):
    r = assess_compliance_bounded(
        atmos, 705e3, resolvable, cd=2.2, area_m2=1.0, isp_s=1666.0,
        propellant_available_kg=50.0, solar_activity="mean",
    )
    assert r.renderable
    assert r.verdict in DECIDED_VERDICTS
    assert r.optimistic is not None and r.pessimistic is not None
    assert r.delta_v_interval_ms.lo <= r.delta_v_interval_ms.hi
    assert r.propellant_interval_kg.lo <= r.propellant_interval_kg.hi
    assert r.ballistic_coefficient.lo < r.ballistic_coefficient.hi


def test_a_verdict_that_flips_across_the_mass_range_is_ambiguous(atmos, resolvable, capsys):
    """At 420 km this design complies naturally if light, needs a burn if heavy.

    The mass interval, not the design, decides. That is not a verdict.
    """
    r = assess_compliance_bounded(
        atmos, 420e3, resolvable, cd=2.2, area_m2=3.0, isp_s=1666.0,
        propellant_available_kg=50.0, solar_activity="mean",
    )
    with capsys.disabled():
        print(f"\n  420 km: low mass -> {r.optimistic.verdict.value}, "
              f"high mass -> {r.pessimistic.verdict.value}")
    assert r.optimistic.verdict is not r.pessimistic.verdict
    assert r.verdict is Compliance.AMBIGUOUS
    assert r.renderable is False
    assert any("flips across the mass range" in n for n in r.notes)


def test_refusal_verdicts_are_excluded_from_the_decided_set():
    assert Compliance.NOT_ASSESSABLE not in DECIDED_VERDICTS
    assert Compliance.AMBIGUOUS not in DECIDED_VERDICTS
    assert Compliance.COMPLIANT_NATURAL in DECIDED_VERDICTS


def test_plain_compliance_results_expose_renderable_too(atmos):
    """The flag lives on the underlying result as well, so neither path can
    hand a UI a refusal without one."""
    from sim.disposal import assess_compliance
    r = assess_compliance(
        atmos, 705e3, 260.0, 2.2, 1.0, 1666.0,
        propellant_available_kg=5.0, solar_activity="mean",
    )
    assert r.renderable is True
    assert r.as_dict()["renderable"] is True
