"""Post-mission disposal: delta-v, propellant, and the compliance verdict.

V2_BRIEF.md §3. Above roughly 550 km "will it decay?" stops being the
interesting question, because the answer is "not for centuries". The question
that replaces it is regulatory: can this design get itself out of orbit inside
the window the rule allows, and can it afford the propellant to do so?

Rules and their citations live in `data/disposal_rules.json`, not here. The
five-year figure is a legal threshold that has already changed once and will
change again; it is data.

---

**A limitation stated up front, because it shapes the whole module.**

`PHYSICS.md` §1 excludes orbital eccentricity: the state is `[a, m]` and every
orbit this simulator propagates is circular. Real post-mission disposal is
almost never circular -- an operator fires once near apogee, drops perigee into
the thermosphere, and lets drag at perigee do the rest over many elliptical
revolutions. That is both cheaper and faster than the circular equivalent.

So there are two delta-v numbers here and they are not interchangeable:

  * `hohmann_lower_circular` -- two burns to a lower *circular* orbit. More
    expensive, but the resulting orbit is one this engine can actually
    propagate, so the decay time behind the verdict is computed rather than
    asserted. **This is what `assess_compliance` uses.**

  * `perigee_lowering_dv` -- one burn to a target perigee. This is what an
    operator would really do, and the delta-v is exact. But the decay of the
    resulting ellipse is *not modelled by this simulator*, so it is reported
    as a reference figure with no decay time attached, never folded into a
    verdict.

The honest summary: the verdict is a conservative bound. A real mission
flying the elliptical profile would need less propellant than this module
demands, so a design this module calls compliant is compliant, while one it
calls non-compliant may still be recoverable with a profile the engine cannot
yet represent. `ComplianceResult.margin_note` carries that caveat with the
result rather than leaving it in a docstring.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable

import numpy as np

from .constants import G0, MU, R_E, REENTRY_ALTITUDE
from .dynamics import derivatives
from .integrator import propagate_adaptive
from .satellite import Outcome

ROOT = Path(__file__).resolve().parent.parent
RULES_PATH = ROOT / "data" / "disposal_rules.json"

SECONDS_PER_YEAR = 365.25 * 86400.0


# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class DisposalRule:
    """One post-mission disposal rule, with the provenance to cite it."""

    id: str
    label: str
    citation: str
    window_years: float
    clock_starts: str
    applies_below_altitude_km: float
    jurisdiction: str
    source: str
    quote: str = ""

    @property
    def window_s(self) -> float:
        return self.window_years * SECONDS_PER_YEAR

    def applies_at(self, altitude_km: float) -> bool:
        return altitude_km < self.applies_below_altitude_km


def load_rules(path: Path | str = RULES_PATH) -> dict[str, DisposalRule]:
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    return {
        r["id"]: DisposalRule(
            id=r["id"],
            label=r["label"],
            citation=r["citation"],
            window_years=float(r["window_years"]),
            clock_starts=r["clock_starts"],
            applies_below_altitude_km=float(r["applies_below_altitude_km"]),
            jurisdiction=r["jurisdiction"],
            source=r["source"],
            quote=r.get("quote", ""),
        )
        for r in payload["rules"]
    }


def default_rule(path: Path | str = RULES_PATH) -> DisposalRule:
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    return load_rules(path)[payload["default_rule"]]


# --------------------------------------------------------------------------
# Orbital mechanics for the disposal manoeuvre
# --------------------------------------------------------------------------

def circular_velocity(r_m: float) -> float:
    """Speed on a circular orbit of radius r. `v = sqrt(MU/r)`."""
    if r_m <= 0.0:
        raise ValueError(f"radius must be positive, got {r_m}")
    return math.sqrt(MU / r_m)


def vis_viva(r_m: float, a_m: float) -> float:
    """Speed at radius r on an orbit of semi-major axis a.

        v = sqrt(MU * (2/r - 1/a))

    The standard vis-viva relation. `a` here is a genuine semi-major axis of
    an ellipse, which is the one place in this codebase where `a` is not also
    the orbital radius -- everywhere else the circular assumption makes them
    equal (PHYSICS.md §2).
    """
    term = 2.0 / r_m - 1.0 / a_m
    if term <= 0.0:
        raise ValueError(
            f"unbound or invalid orbit: r={r_m:.1f} m, a={a_m:.1f} m"
        )
    return math.sqrt(MU * term)


def hohmann_lower_circular(r_from_m: float, r_to_m: float) -> tuple[float, float]:
    """Two-burn Hohmann transfer between circular orbits. Returns (dv1, dv2), m/s.

    Both returned values are magnitudes. Lowering an orbit means both burns
    are retrograde; raising it means both are prograde. The transfer ellipse
    has `a_t = (r_from + r_to)/2` and touches both circles.

    Textbook result, and the reason it is worth writing out rather than
    importing: for a small altitude change the total approaches
    `dv ~ (v/2) * |dr| / r`, which is the sanity check applied in
    `tests/test_disposal.py`.
    """
    if r_from_m <= 0.0 or r_to_m <= 0.0:
        raise ValueError("radii must be positive")
    a_t = 0.5 * (r_from_m + r_to_m)
    dv1 = abs(vis_viva(r_from_m, a_t) - circular_velocity(r_from_m))
    dv2 = abs(circular_velocity(r_to_m) - vis_viva(r_to_m, a_t))
    return dv1, dv2


def perigee_lowering_dv(r_circular_m: float, r_perigee_m: float) -> float:
    """Single retrograde burn lowering perigee from a circular orbit, m/s.

    What an operator actually flies. Reported for reference only -- the
    resulting ellipse is outside this simulator's circular state (see the
    module docstring), so no decay time is attached to it.
    """
    if r_perigee_m > r_circular_m:
        raise ValueError(
            f"target perigee {r_perigee_m:.1f} m is above the starting "
            f"circular radius {r_circular_m:.1f} m -- that is a raise, not a "
            "perigee lowering"
        )
    a_t = 0.5 * (r_circular_m + r_perigee_m)
    return circular_velocity(r_circular_m) - vis_viva(r_circular_m, a_t)


def propellant_mass(dv_ms: float, wet_mass_kg: float, isp_s: float) -> float:
    """Propellant burned for a given delta-v. Tsiolkovsky, solved for mass.

        dv = Isp*g0*ln(m0/m1)   =>   m_prop = m0 * (1 - exp(-dv/(Isp*g0)))

    `wet_mass_kg` is the mass at the *start* of the burn, i.e. at end of
    mission, not at launch.
    """
    if isp_s <= 0.0:
        raise ValueError(f"Isp must be positive, got {isp_s}")
    if wet_mass_kg <= 0.0:
        raise ValueError(f"mass must be positive, got {wet_mass_kg}")
    if dv_ms < 0.0:
        raise ValueError(f"delta-v must be non-negative, got {dv_ms}")
    return wet_mass_kg * (1.0 - math.exp(-dv_ms / (isp_s * G0)))


# --------------------------------------------------------------------------
# Decay time, and the altitude that complies
# --------------------------------------------------------------------------

def decay_time_s(
    rho_of_h: Callable[[float, float], float],
    altitude_m: float,
    mass_kg: float,
    cd: float,
    area_m2: float,
    t_max_s: float,
    tol: float = 1e-4,
    dt_max: float = 30 * 86400.0,
) -> float | None:
    """Time to fall from `altitude_m` to the reentry altitude, or None.

    Thrusters off: disposal decay is drag alone, so this is PHYSICS.md §3.2
    with `F = 0`, the same equation and the same integrator as everything
    else. Returns None if it has not reentered by `t_max_s`.
    """
    def deriv(t: float, y: np.ndarray) -> np.ndarray:
        return derivatives(
            y, rho_of_h(t, y[0] - R_E), thrust=0.0, cd=cd, area=area_m2, isp=None
        )

    traj = propagate_adaptive(
        deriv,
        np.array([R_E + altitude_m, mass_kg]),
        t_max=t_max_s,
        tol=tol,
        dt_max=dt_max,
    )
    if traj.outcome is not Outcome.REENTERED:
        return None
    return float(traj.outcome_time_s)


def highest_complying_altitude_m(
    rho_of_h: Callable[[float, float], float],
    mass_kg: float,
    cd: float,
    area_m2: float,
    window_s: float,
    h_lo_m: float = REENTRY_ALTITUDE,
    h_hi_m: float | None = None,
    tol_m: float = 1e3,
    max_iter: int = 40,
) -> float | None:
    """Highest circular altitude whose natural decay finishes inside the window.

    Bisection. Decay time increases monotonically with altitude for fixed
    ballistic coefficient, so the root is bracketed cleanly -- the same
    argument `critical_altitude` relies on in PHYSICS.md §4.

    Returns None if even `h_lo_m` fails to comply, which would mean the window
    is shorter than the fall from the reentry altitude itself.
    """
    if h_hi_m is None:
        h_hi_m = 2000e3

    def complies(h_m: float) -> bool:
        t = decay_time_s(rho_of_h, h_m, mass_kg, cd, area_m2, t_max_s=window_s)
        return t is not None and t <= window_s

    if not complies(h_lo_m):
        return None
    if complies(h_hi_m):
        return h_hi_m

    for _ in range(max_iter):
        if h_hi_m - h_lo_m < tol_m:
            break
        mid = 0.5 * (h_lo_m + h_hi_m)
        if complies(mid):
            h_lo_m = mid
        else:
            h_hi_m = mid
    return h_lo_m


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

class Compliance(Enum):
    COMPLIANT_NATURAL = "COMPLIANT_NATURAL"
    COMPLIANT_WITH_DISPOSAL = "COMPLIANT_WITH_DISPOSAL"
    NON_COMPLIANT_INSUFFICIENT_PROPELLANT = "NON_COMPLIANT_INSUFFICIENT_PROPELLANT"
    NON_COMPLIANT_NO_SOLUTION = "NON_COMPLIANT_NO_SOLUTION"
    OUT_OF_SCOPE = "OUT_OF_SCOPE"

    # The mass this verdict would rest on could not be resolved (Gate 9), so
    # there is no verdict. Distinct from every outcome above: those are answers,
    # this is a refusal to answer. `ComplianceResult.renderable` is False for it
    # so a UI cannot show it as though it were a result.
    NOT_ASSESSABLE = "NOT_ASSESSABLE"

    # Compliance flips between the ends of the input uncertainty. Also not a
    # verdict -- the honest output is "this design is not determined either way
    # by what we know".
    AMBIGUOUS = "AMBIGUOUS"


#: Verdicts that are answers. Anything else must not be presented as one.
DECIDED_VERDICTS = frozenset({
    Compliance.COMPLIANT_NATURAL,
    Compliance.COMPLIANT_WITH_DISPOSAL,
    Compliance.NON_COMPLIANT_INSUFFICIENT_PROPELLANT,
    Compliance.NON_COMPLIANT_NO_SOLUTION,
    Compliance.OUT_OF_SCOPE,
})


@dataclass
class ComplianceResult:
    """Everything behind the verdict, so it can be argued with rather than believed."""

    verdict: Compliance
    rule_id: str
    rule_label: str
    rule_citation: str
    window_years: float
    solar_activity: str

    operational_altitude_km: float
    natural_decay_years: float | None      # None => does not decay in the window
    disposal_altitude_km: float | None     # target circular altitude, if a burn is needed
    delta_v_ms: float | None
    propellant_required_kg: float | None
    propellant_available_kg: float | None

    # Reported for reference; not used in the verdict. See the module docstring.
    perigee_lowering_delta_v_ms: float | None = None

    margin_note: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def renderable(self) -> bool:
        """Whether a UI may present this as a verdict.

        False for NOT_ASSESSABLE and AMBIGUOUS. A compliance verdict resting on
        a mass the model could not resolve is not a weak verdict, it is not a
        verdict, and showing it next to the ones that are would be the single
        most misleading thing this tool could do.
        """
        return self.verdict in DECIDED_VERDICTS

    def as_dict(self) -> dict:
        return {
            "verdict": self.verdict.value,
            "renderable": self.renderable,
            "rule": {
                "id": self.rule_id,
                "label": self.rule_label,
                "citation": self.rule_citation,
                "window_years": self.window_years,
            },
            "solar_activity": self.solar_activity,
            "operational_altitude_km": self.operational_altitude_km,
            "natural_decay_years": self.natural_decay_years,
            "disposal_altitude_km": self.disposal_altitude_km,
            "delta_v_ms": self.delta_v_ms,
            "propellant_required_kg": self.propellant_required_kg,
            "propellant_available_kg": self.propellant_available_kg,
            "perigee_lowering_delta_v_ms": self.perigee_lowering_delta_v_ms,
            "margin_note": self.margin_note,
            "notes": self.notes,
        }


_CIRCULAR_CAVEAT = (
    "Delta-v is for a two-burn transfer to a lower circular orbit, because "
    "that is the orbit this engine can propagate (PHYSICS.md §1 excludes "
    "eccentricity). Real disposal lowers perigee with a single burn and costs "
    "less, so this is a conservative bound, not an estimate."
)


def assess_compliance(
    rho_of_h: Callable[[float, float], float],
    operational_altitude_m: float,
    mass_at_eol_kg: float,
    cd: float,
    area_m2: float,
    isp_s: float,
    propellant_available_kg: float,
    solar_activity: str,
    rule: DisposalRule | None = None,
    target_perigee_m: float = 200e3,
) -> ComplianceResult:
    """Does this design get itself out of orbit in time, and can it afford to?

    V2_BRIEF.md §3's compliance calculator. The sequence is:

      1. Is the rule even in scope at this altitude?
      2. Does natural decay comply? If so, no propellant is needed.
      3. If not, find the highest circular altitude that *does* comply, and
         price the transfer down to it.
      4. Compare against the propellant the design actually carries.

    Every number in the result is computed here or by the propagator. Nothing
    is looked up from a table of typical values.
    """
    if rule is None:
        rule = default_rule()

    h_op_km = operational_altitude_m / 1e3
    notes: list[str] = []

    if not rule.applies_at(h_op_km):
        return ComplianceResult(
            verdict=Compliance.OUT_OF_SCOPE,
            rule_id=rule.id, rule_label=rule.label, rule_citation=rule.citation,
            window_years=rule.window_years, solar_activity=solar_activity,
            operational_altitude_km=h_op_km,
            natural_decay_years=None, disposal_altitude_km=None,
            delta_v_ms=None, propellant_required_kg=None,
            propellant_available_kg=propellant_available_kg,
            notes=[
                f"{rule.label} applies below "
                f"{rule.applies_below_altitude_km:.0f} km; this orbit is at "
                f"{h_op_km:.1f} km. Disposal above the LEO region is a "
                "different regime (graveyard orbits) and is not modelled."
            ],
        )

    natural_s = decay_time_s(
        rho_of_h, operational_altitude_m, mass_at_eol_kg, cd, area_m2,
        t_max_s=rule.window_s,
    )
    if natural_s is not None and natural_s <= rule.window_s:
        return ComplianceResult(
            verdict=Compliance.COMPLIANT_NATURAL,
            rule_id=rule.id, rule_label=rule.label, rule_citation=rule.citation,
            window_years=rule.window_years, solar_activity=solar_activity,
            operational_altitude_km=h_op_km,
            natural_decay_years=natural_s / SECONDS_PER_YEAR,
            disposal_altitude_km=None, delta_v_ms=0.0,
            propellant_required_kg=0.0,
            propellant_available_kg=propellant_available_kg,
            notes=[
                f"Natural decay reenters in "
                f"{natural_s / SECONDS_PER_YEAR:.2f} years, inside the "
                f"{rule.window_years:.0f}-year window. No disposal burn needed."
            ],
        )

    h_target_m = highest_complying_altitude_m(
        rho_of_h, mass_at_eol_kg, cd, area_m2, rule.window_s,
        h_hi_m=operational_altitude_m,
    )
    if h_target_m is None:
        return ComplianceResult(
            verdict=Compliance.NON_COMPLIANT_NO_SOLUTION,
            rule_id=rule.id, rule_label=rule.label, rule_citation=rule.citation,
            window_years=rule.window_years, solar_activity=solar_activity,
            operational_altitude_km=h_op_km,
            natural_decay_years=None, disposal_altitude_km=None,
            delta_v_ms=None, propellant_required_kg=None,
            propellant_available_kg=propellant_available_kg,
            notes=[
                "No circular altitude down to the reentry altitude decays "
                f"inside {rule.window_years:.0f} years at this ballistic "
                "coefficient. Check the mass and area inputs."
            ],
        )

    r_op = R_E + operational_altitude_m
    r_target = R_E + h_target_m
    dv1, dv2 = hohmann_lower_circular(r_op, r_target)
    dv = dv1 + dv2
    required = propellant_mass(dv, mass_at_eol_kg, isp_s)

    dv_perigee = perigee_lowering_dv(r_op, R_E + target_perigee_m)

    notes.append(
        f"Natural decay does not reenter within {rule.window_years:.0f} years "
        f"at {h_op_km:.1f} km."
    )
    notes.append(
        f"Highest complying circular altitude is {h_target_m / 1e3:.1f} km; "
        f"transfer costs {dv1:.1f} + {dv2:.1f} = {dv:.1f} m/s."
    )
    notes.append(
        f"For reference, a single-burn drop to a {target_perigee_m / 1e3:.0f} km "
        f"perigee costs {dv_perigee:.1f} m/s, but the resulting ellipse is "
        "outside this engine's circular state and carries no computed decay time."
    )

    enough = required <= propellant_available_kg
    verdict = (
        Compliance.COMPLIANT_WITH_DISPOSAL if enough
        else Compliance.NON_COMPLIANT_INSUFFICIENT_PROPELLANT
    )
    if not enough:
        notes.append(
            f"Design carries {propellant_available_kg:.2f} kg but needs "
            f"{required:.2f} kg -- short by "
            f"{required - propellant_available_kg:.2f} kg."
        )

    return ComplianceResult(
        verdict=verdict,
        rule_id=rule.id, rule_label=rule.label, rule_citation=rule.citation,
        window_years=rule.window_years, solar_activity=solar_activity,
        operational_altitude_km=h_op_km,
        natural_decay_years=None,
        disposal_altitude_km=h_target_m / 1e3,
        delta_v_ms=dv,
        propellant_required_kg=required,
        propellant_available_kg=propellant_available_kg,
        perigee_lowering_delta_v_ms=dv_perigee,
        margin_note=_CIRCULAR_CAVEAT,
        notes=notes,
    )
