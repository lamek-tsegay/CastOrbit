"""Propagating the mass interval through to the compliance verdict.

Gate 9 failed at 1/3 (`docs/PHYSICS_V2.md` §V2.7). The response is not a better
mass model -- V2_BRIEF.md §5 puts that in stage 2 -- but an honest boundary
around the one we have. `sim/mass_model.py` now withholds a point estimate when
the table cannot support one. This module carries that decision downstream so
it cannot be quietly dropped between the mass and the answer.

**The chain, and how uncertainty moves along it:**

```
mass interval  ->  ballistic coefficient  ->  decay time  ->  compliance
   [lo, hi]         Bc = m/(Cd*A)            propagated       verdict at
                    monotone increasing      at both ends     both ends
```

Every step is monotone in mass, which is what makes interval propagation exact
here rather than approximate:

  * **Heavier means a higher ballistic coefficient.** `Bc = m/(Cd*A)` rises
    with mass, so a heavier satellite carries more inertia per unit of drag
    area and decelerates less.
  * **A higher ballistic coefficient means slower decay.** The drag term in
    `PHYSICS.md` §3.2 is `rho*(Cd*A/m)*sqrt(MU*a)`, inversely proportional to
    mass.
  * **Slower decay means compliance is harder.** So the *low*-mass end is the
    optimistic case and the *high*-mass end the pessimistic one.

Because the mapping is monotone, evaluating the two endpoints bounds every
interior value; there is no need to sample inside the interval and no risk of
missing a worse case in the middle.

**The refusal rule.** If the mass is unresolvable, no verdict is produced --
not a hedged one, not a most-likely one. If the mass is resolvable but the two
endpoints disagree, the verdict is AMBIGUOUS. Both are non-renderable
(`ComplianceResult.renderable`), so a UI has a single boolean to check and no
way to accidentally present a refusal as an answer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from .disposal import (
    SECONDS_PER_YEAR,
    Compliance,
    ComplianceResult,
    DisposalRule,
    assess_compliance,
    decay_time_s,
    default_rule,
)
from .mass_model import MassEstimate


@dataclass(frozen=True)
class Interval:
    """A closed range. Ordering is enforced, not assumed."""

    lo: float
    hi: float

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError(f"inverted interval [{self.lo}, {self.hi}]")

    @property
    def width(self) -> float:
        return self.hi - self.lo

    @property
    def ratio(self) -> float:
        return float("inf") if self.lo == 0 else self.hi / self.lo

    def contains(self, value: float) -> bool:
        return self.lo <= value <= self.hi

    def as_list(self) -> list[float]:
        return [self.lo, self.hi]

    def __repr__(self) -> str:
        return f"[{self.lo:.4g}, {self.hi:.4g}]"


def ballistic_coefficient_interval(
    mass: Interval, cd: float, area_m2: float
) -> Interval:
    """`Bc = m / (Cd*A)`, in kg/m^2. Monotone increasing in mass.

    Cd and area are treated as exact here. They are not -- Phase 8 gives `A` a
    range of its own -- but composing both uncertainties is a wider change than
    bounding the mass, and conflating them would hide which one dominates.
    """
    if cd <= 0 or area_m2 <= 0:
        raise ValueError("cd and area must be positive")
    denom = cd * area_m2
    return Interval(mass.lo / denom, mass.hi / denom)


def decay_time_interval(
    rho_of_h: Callable[[float, float], float],
    altitude_m: float,
    mass: Interval,
    cd: float,
    area_m2: float,
    t_max_s: float,
) -> tuple[Interval | None, dict]:
    """Time to reentry at both mass endpoints.

    Returns `(interval_or_None, detail)`. The interval is None when at least
    one endpoint does not reenter inside `t_max_s` -- an unbounded-above decay
    time is not a number, and returning `t_max` in its place would silently
    convert "never" into "just in time".

    Lighter decays faster, so the low-mass endpoint gives the lower bound.
    """
    t_light = decay_time_s(rho_of_h, altitude_m, mass.lo, cd, area_m2, t_max_s)
    t_heavy = decay_time_s(rho_of_h, altitude_m, mass.hi, cd, area_m2, t_max_s)

    detail = {
        "t_at_low_mass_s": t_light,
        "t_at_high_mass_s": t_heavy,
        "low_mass_kg": mass.lo,
        "high_mass_kg": mass.hi,
    }
    if t_light is None or t_heavy is None:
        detail["note"] = (
            "at least one mass endpoint does not reenter within the window, so "
            "the decay time is unbounded above and no interval is reported"
        )
        return None, detail

    # Monotonicity is a physical claim, so check it rather than trust it.
    if t_light > t_heavy:
        detail["note"] = (
            "decay time was not monotone in mass, which should be impossible "
            "for fixed Cd and area; endpoints sorted, treat with suspicion"
        )
    return Interval(min(t_light, t_heavy), max(t_light, t_heavy)), detail


@dataclass
class BoundedComplianceResult:
    """A compliance verdict that knows what it rests on.

    `renderable` is the only thing a UI needs to check. It is False whenever
    the verdict is a refusal rather than an answer.
    """

    verdict: Compliance
    mass_estimate: MassEstimate
    mass_interval_kg: Interval
    optimistic: ComplianceResult | None = None   # low mass, decays fastest
    pessimistic: ComplianceResult | None = None  # high mass, decays slowest
    ballistic_coefficient: Interval | None = None
    delta_v_interval_ms: Interval | None = None
    propellant_interval_kg: Interval | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def renderable(self) -> bool:
        return self.verdict in {
            Compliance.COMPLIANT_NATURAL,
            Compliance.COMPLIANT_WITH_DISPOSAL,
            Compliance.NON_COMPLIANT_INSUFFICIENT_PROPELLANT,
            Compliance.NON_COMPLIANT_NO_SOLUTION,
            Compliance.OUT_OF_SCOPE,
        }

    def as_dict(self) -> dict:
        return {
            "verdict": self.verdict.value,
            "renderable": self.renderable,
            "mass": self.mass_estimate.as_dict(),
            "mass_interval_kg": self.mass_interval_kg.as_list(),
            "ballistic_coefficient_kg_m2": (
                None if self.ballistic_coefficient is None
                else self.ballistic_coefficient.as_list()
            ),
            "delta_v_interval_ms": (
                None if self.delta_v_interval_ms is None
                else self.delta_v_interval_ms.as_list()
            ),
            "propellant_interval_kg": (
                None if self.propellant_interval_kg is None
                else self.propellant_interval_kg.as_list()
            ),
            "optimistic": None if self.optimistic is None else self.optimistic.as_dict(),
            "pessimistic": (
                None if self.pessimistic is None else self.pessimistic.as_dict()
            ),
            "notes": self.notes,
        }


def assess_compliance_bounded(
    rho_of_h: Callable[[float, float], float],
    operational_altitude_m: float,
    mass_estimate: MassEstimate,
    cd: float,
    area_m2: float,
    isp_s: float,
    propellant_available_kg: float,
    solar_activity: str,
    rule: DisposalRule | None = None,
) -> BoundedComplianceResult:
    """Compliance verdict carrying the mass uncertainty it rests on.

    **Refuses outright when the mass is unresolvable.** No verdict, no
    endpoints, no "most likely" -- `renderable` is False and the reasons the
    mass model gave are passed through verbatim. This is the requirement that
    a verdict built on an unbounded mass estimate must not render.
    """
    if rule is None:
        rule = default_rule()
    mass = Interval(*mass_estimate.interval_kg)

    if not mass_estimate.resolvable:
        return BoundedComplianceResult(
            verdict=Compliance.NOT_ASSESSABLE,
            mass_estimate=mass_estimate,
            mass_interval_kg=mass,
            ballistic_coefficient=ballistic_coefficient_interval(mass, cd, area_m2),
            notes=[
                "No compliance verdict: the dry mass could not be resolved to a "
                "single value, so any verdict would be an artefact of an "
                f"arbitrary choice within {mass}.",
                *mass_estimate.refusal_reasons,
                f"The mass is bounded to {mass} kg. That range is the honest "
                "output; narrowing it needs a better mass model (V2_BRIEF.md §5 "
                "stage 2), not a better propagator.",
            ],
        )

    ends = {}
    for label, m in (("optimistic", mass.lo), ("pessimistic", mass.hi)):
        ends[label] = assess_compliance(
            rho_of_h,
            operational_altitude_m,
            mass_at_eol_kg=m,
            cd=cd,
            area_m2=area_m2,
            isp_s=isp_s,
            propellant_available_kg=propellant_available_kg,
            solar_activity=solar_activity,
            rule=rule,
        )
    lo_res, hi_res = ends["optimistic"], ends["pessimistic"]

    dv = _interval_of(lo_res.delta_v_ms, hi_res.delta_v_ms)
    prop = _interval_of(lo_res.propellant_required_kg, hi_res.propellant_required_kg)
    bc = ballistic_coefficient_interval(mass, cd, area_m2)

    notes: list[str] = [
        f"Mass bounded to {mass} kg ({mass_estimate.provenance}); "
        f"ballistic coefficient {bc} kg/m2.",
    ]
    notes.extend(mass_estimate.warnings)

    if lo_res.verdict is hi_res.verdict:
        verdict = lo_res.verdict
    else:
        verdict = Compliance.AMBIGUOUS
        notes.insert(0, (
            f"No verdict: compliance flips across the mass range. At "
            f"{mass.lo:.0f} kg the design is {lo_res.verdict.value}; at "
            f"{mass.hi:.0f} kg it is {hi_res.verdict.value}. The answer is "
            "determined by mass uncertainty, not by the design."
        ))

    return BoundedComplianceResult(
        verdict=verdict,
        mass_estimate=mass_estimate,
        mass_interval_kg=mass,
        optimistic=lo_res,
        pessimistic=hi_res,
        ballistic_coefficient=bc,
        delta_v_interval_ms=dv,
        propellant_interval_kg=prop,
        notes=notes,
    )


def _interval_of(a: float | None, b: float | None) -> Interval | None:
    if a is None or b is None:
        return None
    return Interval(min(a, b), max(a, b))


def natural_decay_years_interval(
    rho_of_h: Callable[[float, float], float],
    altitude_m: float,
    mass_estimate: MassEstimate,
    cd: float,
    area_m2: float,
    horizon_years: float = 200.0,
) -> tuple[Interval | None, dict]:
    """Natural decay time in years, bounded by the mass interval.

    Available whether or not the mass resolves to a point -- a bounded decay
    time is useful even when a verdict is not, and it is the number a designer
    actually wants when the tool declines to rule.
    """
    interval, detail = decay_time_interval(
        rho_of_h,
        altitude_m,
        Interval(*mass_estimate.interval_kg),
        cd,
        area_m2,
        t_max_s=horizon_years * SECONDS_PER_YEAR,
    )
    detail["mass_resolvable"] = mass_estimate.resolvable
    if interval is None:
        return None, detail
    return (
        Interval(interval.lo / SECONDS_PER_YEAR, interval.hi / SECONDS_PER_YEAR),
        detail,
    )
