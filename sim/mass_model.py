"""Stage-1 mass closure: interpolation between real spacecraft. V2_BRIEF.md §5.

**What this is allowed to be.** V2_BRIEF.md §5 is explicit that real mass
estimating relationships are weeks of work and that stage 1 should not attempt
them. So this does exactly one thing: given a power requirement, find the two
real spacecraft in `data/reference_satellites.json` that bracket it and
interpolate between them. Nothing is fitted, nothing is regressed, and no
coefficient appears anywhere in this file.

**What it must never be.** An LLM producing a mass number. The prose-to-spec
path (Phase 10) may set the *requirement* -- "about 500 W, earth observation,
LEO" -- but the kilogram figure has to come from this table or it does not
exist.

---

**Why power is the predictor.** Spacecraft dry mass tracks electrical power
more tightly than it tracks anything else that is knowable at concept stage:
the power system (arrays, batteries, PMAD, harness, and the structure and
radiators to carry them) is a large and fairly constant fraction of the bus,
and payload power scales with payload mass. Bus dimensions are the obvious
alternative but are published as deployed spans about as often as stowed
envelopes, which makes them treacherous.

Interpolation is in log-log space. Mass spans 5 kg to 3650 kg and power 20 W
to 6400 W across this table -- nearly three orders of magnitude each -- so
linear interpolation between neighbours would be dominated by the large end.

**Matching rules, in order:**

  1. Same `payload_class`. Communications satellites are much lighter per watt
     than instrument-carrying ones (Intelsat 15 is 0.27 kg/W; CryoSat-2 is
     0.81 kg/W). Mixing them is the single biggest avoidable error.
  2. Never the same `family`. Predicting GOES-16 from GOES-17 would be
     interpolating a satellite from its own twin.
  3. Dry masses only. An entry whose source publishes launch mass is excluded
     from the interpolation, because for a GEO spacecraft that is roughly
     double the answer.

Every result carries the spacecraft it came from, the weights, and an
`estimated` tag. A mass with no provenance is not a product of this module.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFERENCE_PATH = ROOT / "data" / "reference_satellites.json"


@dataclass(frozen=True)
class ReferenceSatellite:
    """One row of the table. Immutable -- nothing here is ever adjusted."""

    id: str
    name: str
    family: str
    payload_class: str
    altitude_class: str
    mass_kg: float
    mass_kind: str          # "dry" | "launch"
    mass_confidence: str
    power_w: float | None
    source: str
    altitude_km: float | None = None
    solar_array_m2: float | None = None
    power_confidence: str | None = None
    held_out: bool = False

    @property
    def is_dry(self) -> bool:
        return self.mass_kind == "dry"

    @property
    def usable_for_power_interpolation(self) -> bool:
        return self.is_dry and self.power_w is not None and self.power_w > 0

    @property
    def kg_per_w(self) -> float | None:
        if not self.usable_for_power_interpolation:
            return None
        return self.mass_kg / self.power_w


def load_reference_satellites(
    path: Path | str = REFERENCE_PATH,
) -> list[ReferenceSatellite]:
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    out = []
    for row in payload["satellites"]:
        out.append(
            ReferenceSatellite(
                id=row["id"],
                name=row["name"],
                family=row["family"],
                payload_class=row["payload_class"],
                altitude_class=row["altitude_class"],
                mass_kg=float(row["mass_kg"]),
                mass_kind=row["mass_kind"],
                mass_confidence=row["mass_confidence"],
                power_w=None if row.get("power_w") is None else float(row["power_w"]),
                source=row["source"],
                altitude_km=row.get("altitude_km"),
                solar_array_m2=row.get("solar_array_m2"),
                power_confidence=row.get("power_confidence"),
                held_out=bool(row.get("held_out", False)),
            )
        )
    return out


def training_set(
    satellites: list[ReferenceSatellite] | None = None,
) -> list[ReferenceSatellite]:
    """The table the interpolator is allowed to see. Held-out entries removed."""
    sats = satellites if satellites is not None else load_reference_satellites()
    return [s for s in sats if not s.held_out]


def held_out_set(
    satellites: list[ReferenceSatellite] | None = None,
) -> list[ReferenceSatellite]:
    sats = satellites if satellites is not None else load_reference_satellites()
    return [s for s in sats if s.held_out]


@dataclass
class MassEstimate:
    """A dry mass, and everything needed to argue with it.

    `tag` is always "estimated". There is no path through this module that
    produces an untagged mass, because V2_BRIEF.md §5 requires a generated
    design to carry its provenance.
    """

    mass_kg: float
    tag: str
    method: str
    neighbours: list[dict] = field(default_factory=list)
    power_w: float | None = None
    payload_class: str | None = None
    extrapolated: bool = False
    bracket_ratio: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def provenance(self) -> str:
        """One sentence naming the spacecraft this was interpolated between."""
        if not self.neighbours:
            return "no reference spacecraft"
        names = " and ".join(
            f"{n['name']} ({n['mass_kg']:.0f} kg at {n['power_w']:.0f} W)"
            for n in self.neighbours
        )
        verb = "extrapolated beyond" if self.extrapolated else "interpolated between"
        return f"{verb} {names}"

    def as_dict(self) -> dict:
        return {
            "mass_kg": self.mass_kg,
            "tag": self.tag,
            "method": self.method,
            "provenance": self.provenance,
            "neighbours": self.neighbours,
            "power_w": self.power_w,
            "payload_class": self.payload_class,
            "extrapolated": self.extrapolated,
            "bracket_ratio": self.bracket_ratio,
            "warnings": self.warnings,
        }


class MassEstimationError(RuntimeError):
    """Raised when the table cannot support an estimate.

    Deliberately an error rather than a fallback guess. V2_BRIEF.md §6:
    nothing hardcoded, nothing faked. A design tool that cannot estimate a
    mass must say so, not invent one.
    """


def estimate_dry_mass(
    power_w: float,
    payload_class: str,
    satellites: list[ReferenceSatellite] | None = None,
    exclude_family: str | None = None,
    allow_class_fallback: bool = True,
) -> MassEstimate:
    """Estimate dry mass by interpolating between two real spacecraft.

    Args:
        power_w: required electrical power. The predictor.
        payload_class: "earth_observation", "communications", ... Matched
            first, because it is the largest single discriminator.
        exclude_family: skip this family, to prevent a satellite being
            predicted from its own twin.
        allow_class_fallback: if the payload class has too few usable entries
            to bracket, widen to the whole table and record a warning. Set
            False to make a thin class an error instead.

    Raises:
        MassEstimationError: if no two usable spacecraft can be found.
    """
    if power_w <= 0:
        raise ValueError(f"power must be positive, got {power_w}")

    pool = training_set(satellites)
    usable = [
        s for s in pool
        if s.usable_for_power_interpolation and s.family != exclude_family
    ]
    if not usable:
        raise MassEstimationError("no reference spacecraft with dry mass and power")

    warnings: list[str] = []
    in_class = [s for s in usable if s.payload_class == payload_class]

    if len(in_class) >= 2:
        candidates = in_class
    elif allow_class_fallback:
        warnings.append(
            f"only {len(in_class)} usable reference spacecraft in payload class "
            f"'{payload_class}'; widened to all classes. Communications and "
            "instrument-carrying spacecraft differ by ~3x in kg/W, so treat "
            "this estimate as weak."
        )
        candidates = usable
    else:
        raise MassEstimationError(
            f"payload class '{payload_class}' has {len(in_class)} usable "
            "reference spacecraft; need 2"
        )

    candidates = sorted(candidates, key=lambda s: s.power_w)

    below = [s for s in candidates if s.power_w <= power_w]
    above = [s for s in candidates if s.power_w > power_w]

    if below and above:
        lo, hi = below[-1], above[0]
        extrapolated = False
    else:
        # Outside the table. Use the two nearest and say so loudly -- a power
        # law extended past its last data point is the least reliable thing
        # this module can produce.
        extrapolated = True
        lo, hi = (candidates[0], candidates[1]) if not below else (
            candidates[-2], candidates[-1]
        )
        warnings.append(
            f"{power_w:.0f} W is outside the reference range "
            f"{candidates[0].power_w:.0f}-{candidates[-1].power_w:.0f} W for this "
            "class; the result is extrapolated, not interpolated."
        )

    mass = _loglog_interpolate(
        power_w, lo.power_w, lo.mass_kg, hi.power_w, hi.mass_kg
    )

    bracket_ratio = hi.power_w / lo.power_w if lo.power_w > 0 else math.inf
    if bracket_ratio > 4.0 and not extrapolated:
        warnings.append(
            f"the bracketing spacecraft are {bracket_ratio:.0f}x apart in power "
            f"({lo.name} at {lo.power_w:.0f} W, {hi.name} at {hi.power_w:.0f} W). "
            "A wide bracket spans different design regimes and the interpolation "
            "is correspondingly weak -- this is a gap in the table, not a "
            "property of the method."
        )

    return MassEstimate(
        mass_kg=mass,
        tag="estimated",
        method="log-log interpolation on power between two reference spacecraft",
        neighbours=[
            {
                "id": s.id, "name": s.name, "mass_kg": s.mass_kg,
                "power_w": s.power_w, "kg_per_w": s.kg_per_w, "source": s.source,
            }
            for s in (lo, hi)
        ],
        power_w=power_w,
        payload_class=payload_class,
        extrapolated=extrapolated,
        bracket_ratio=bracket_ratio,
        warnings=warnings,
    )


def _loglog_interpolate(
    x: float, x0: float, y0: float, x1: float, y1: float
) -> float:
    """Linear interpolation in log-log space -- a local power law through two points.

    Equivalent to `y = y0 * (x/x0)**alpha` with
    `alpha = log(y1/y0) / log(x1/x0)`, i.e. the power law that passes through
    both reference spacecraft exactly.
    """
    if x0 <= 0 or x1 <= 0 or y0 <= 0 or y1 <= 0:
        raise ValueError("log-log interpolation needs positive values")
    if x1 == x0:
        return 0.5 * (y0 + y1)
    t = (math.log(x) - math.log(x0)) / (math.log(x1) - math.log(x0))
    return math.exp(math.log(y0) + t * (math.log(y1) - math.log(y0)))


def score_held_out(
    satellites: list[ReferenceSatellite] | None = None,
) -> list[dict]:
    """Predict every held-out spacecraft and report the error. Gate 9.

    Each is predicted with its own family excluded, from a table that never
    contained it.
    """
    sats = satellites if satellites is not None else load_reference_satellites()
    rows = []
    for target in held_out_set(sats):
        est = estimate_dry_mass(
            power_w=target.power_w,
            payload_class=target.payload_class,
            satellites=sats,
            exclude_family=target.family,
        )
        error = (est.mass_kg - target.mass_kg) / target.mass_kg
        rows.append({
            "id": target.id,
            "name": target.name,
            "actual_kg": target.mass_kg,
            "predicted_kg": est.mass_kg,
            "error_frac": error,
            "within_25pct": abs(error) <= 0.25,
            "power_w": target.power_w,
            "payload_class": target.payload_class,
            "estimate": est.as_dict(),
        })
    return rows


def population_scatter(
    payload_class: str = "earth_observation",
    satellites: list[ReferenceSatellite] | None = None,
) -> dict:
    """How much do real spacecraft at similar power differ in mass?

    This bounds what *any* power-only predictor can achieve, independent of
    the interpolation scheme, so it is the honest way to read a miss: an error
    inside the population scatter is the predictor hitting its ceiling, not
    the method being wrong.

    Uses launch-mass entries too, deliberately. They are excluded from
    prediction because launch mass overstates dry mass, but for measuring
    spread at a given power they are evidence: two spacecraft of very
    different mass at the same power demonstrate the scatter whichever mass
    convention is quoted.
    """
    sats = satellites if satellites is not None else load_reference_satellites()
    members = [
        s for s in sats
        if s.payload_class == payload_class and s.power_w and s.power_w > 0
    ]
    ratios = sorted((s.mass_kg / s.power_w, s.name, s.power_w) for s in members)

    # Closest pair in power, and how far apart they are in mass.
    by_power = sorted(members, key=lambda s: s.power_w)
    worst_pair = None
    for a, b in zip(by_power, by_power[1:]):
        power_gap = b.power_w / a.power_w
        if power_gap > 1.25:
            continue  # not "the same power" in any useful sense
        mass_ratio = max(a.mass_kg, b.mass_kg) / min(a.mass_kg, b.mass_kg)
        if worst_pair is None or mass_ratio > worst_pair["mass_ratio"]:
            worst_pair = {
                "a": a.name, "b": b.name,
                "power_w": [a.power_w, b.power_w],
                "mass_kg": [a.mass_kg, b.mass_kg],
                "power_gap": power_gap,
                "mass_ratio": mass_ratio,
            }

    return {
        "payload_class": payload_class,
        "n": len(members),
        "kg_per_w_min": ratios[0][0] if ratios else None,
        "kg_per_w_max": ratios[-1][0] if ratios else None,
        "kg_per_w_spread": (ratios[-1][0] / ratios[0][0]) if ratios else None,
        "closest_power_pair": worst_pair,
    }
