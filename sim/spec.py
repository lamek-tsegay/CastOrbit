"""Mission spec: the contract between prose and the engine. V2_BRIEF.md §7, Phase 10.

A design tool with a natural-language front door has one failure mode that
matters more than all the others: **the language model inventing a number the
engine was supposed to compute.** A plausible-sounding "dry mass: 240 kg" in a
generated spec would flow straight into the ballistic coefficient and out the
other end as a compliance verdict, and nothing downstream would know it was
fiction. That is precisely the Lovable-mockup failure `V2_BRIEF.md` §8 exists
to prevent, arriving through the front door instead of the back.

So this module is a boundary, not a parser. It defines:

  * `MissionSpec` -- the *requirements* a design states. Altitude, inclination,
    payload class, power, duration. Things a person asks for.
  * `FORBIDDEN_FIELDS` -- the quantities the engine derives. Mass, area,
    ballistic coefficient, decay time, delta-v, verdicts. Things a person
    receives. Any of these appearing in a spec is rejected outright, whatever
    produced it.
  * `Question` -- what comes back instead of a value when the prose does not
    determine one.

**The guarantee lives in `validate_spec_payload`, not in the extractor.**
Whether a spec came from a language model, a regex, a web form or a file, it
goes through the same validation, and the same rejection. Nothing about the
safety of this boundary depends on the model behaving.

`sim/prose.py` holds the extraction; this file holds the contract.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone

# --------------------------------------------------------------------------
# What the engine computes, and therefore what a spec may never contain
# --------------------------------------------------------------------------

#: Quantities derived by the engine. A spec asserting any of these is rejected.
#: Keyed by field name; the value explains which component owns it, so a
#: rejection message can point at the real source rather than just refusing.
FORBIDDEN_FIELDS: dict[str, str] = {
    "mass_kg": "sim/mass_model.py estimates dry mass from power and payload class",
    "dry_mass_kg": "sim/mass_model.py estimates dry mass",
    "wet_mass_kg": "derived from dry mass plus the propellant budget",
    "bus_mass_kg": "sim/mass_model.py estimates dry mass",
    "area_m2": "sim/geometry.py solves projected area from geometry and attitude",
    "ram_area_m2": "sim/geometry.py solves projected area",
    "cross_section_m2": "sim/geometry.py solves projected area",
    "cd_times_area_m2": "the product of a solved area and the drag convention",
    "ballistic_coefficient": "m/(Cd*A), derived from mass and area",
    "decay_time_years": "sim/disposal.py propagates decay; it is an output",
    "decay_time_s": "sim/disposal.py propagates decay",
    "reentry_date": "sim/disposal.py propagates decay",
    "delta_v_ms": "sim/disposal.py computes disposal delta-v",
    "propellant_required_kg": "sim/disposal.py computes it from delta-v",
    "compliance": "sim/disposal.py returns the verdict",
    "verdict": "sim/disposal.py returns the verdict",
    "compliant": "sim/disposal.py returns the verdict",
    "density": "sim/atmosphere.py, from NRLMSIS",
    "rho": "sim/atmosphere.py, from NRLMSIS",
    "cd": "a modelling convention from data/satellite_specs.json, not a requirement",
    "drag_coefficient": "a modelling convention from data/satellite_specs.json",
}

VALID_PAYLOAD_CLASSES = frozenset({
    "earth_observation", "communications", "technology", "navigation",
})

#: Definitional inclinations only -- ones that are what the name *means*, not
#: ones that require a calculation. "Sun-synchronous" is deliberately absent:
#: its inclination depends on altitude through J2, which PHYSICS.md §1
#: excludes, so this tool asks for it rather than deriving it from physics it
#: does not model.
NAMED_INCLINATIONS_DEG: dict[str, float] = {
    "equatorial": 0.0,
    "polar": 90.0,
}

# Bounds are sanity limits, not opinions about good design. They exist to
# catch a unit error or a hallucinated order of magnitude, not to police
# choices.
BOUNDS: dict[str, tuple[float, float]] = {
    "altitude_km": (150.0, 2000.0),      # below 150 km nothing is an orbit;
                                          # 2000 km is the LEO edge in 47 CFR 25.283(e)
    "inclination_deg": (0.0, 90.0),      # retrograde is passed as 180 - i
    "power_w": (1.0, 30000.0),
    "mission_duration_years": (0.01, 30.0),
    "propellant_available_kg": (0.0, 5000.0),
    "isp_s": (50.0, 5000.0),
    # Geometry. These are design *inputs* -- a chassis is something you choose,
    # like power. The engine computes projected area FROM them (sim/geometry.py);
    # it is the area, not the dimensions, that a spec may never assert.
    "bus_length_m": (0.05, 20.0),
    "bus_width_m": (0.05, 20.0),
    "bus_thickness_m": (0.01, 10.0),
    "solar_array_span_m": (0.05, 60.0),
    "solar_array_chord_m": (0.05, 20.0),
}

#: Optional geometry inputs. All-or-nothing for the bus: a partial chassis
#: cannot be solved, and filling the gap with a typical value is exactly the
#: fabrication this schema exists to prevent.
BUS_FIELDS = ("bus_length_m", "bus_width_m", "bus_thickness_m")
ARRAY_FIELDS = ("solar_array_span_m", "solar_array_chord_m")

REQUIRED_FIELDS = (
    "altitude_km", "inclination_deg", "payload_class", "power_w",
    "mission_duration_years",
)


class SpecError(ValueError):
    """A spec that cannot be trusted. Raised rather than repaired."""


@dataclass(frozen=True)
class Question:
    """What the front door returns instead of inventing a value.

    `why` is not decoration. A user asked for a satellite and got a question
    back; they are owed the reason, and "the tool would otherwise have to make
    it up" is the reason.
    """

    field: str
    question: str
    why: str

    def as_dict(self) -> dict:
        return {"field": self.field, "question": self.question, "why": self.why}


@dataclass
class MissionSpec:
    """Requirements a design states. Nothing the engine derives.

    `provenance` records, per field, where the value came from -- `"stated"`
    for something the prose gave explicitly, `"definitional"` for a named
    orbit whose inclination is what the name means. There is no `"inferred"`
    or `"typical"`: if a value is not in one of those two categories it is not
    in the spec, and a `Question` was returned instead.
    """

    altitude_km: float
    inclination_deg: float
    payload_class: str
    power_w: float
    mission_duration_years: float
    epoch: str | None = None
    propellant_available_kg: float | None = None
    isp_s: float | None = None
    target_altitude_km: float | None = None
    bus_length_m: float | None = None
    bus_width_m: float | None = None
    bus_thickness_m: float | None = None
    solar_array_span_m: float | None = None
    solar_array_chord_m: float | None = None
    label: str = "unnamed"
    provenance: dict[str, str] = field(default_factory=dict)

    @property
    def has_geometry(self) -> bool:
        """Whether the spec determines a chassis the area solver can use.

        False means `sim/geometry.py` has nothing to solve and ram area is
        refused rather than assumed -- the same posture the mass model takes
        when the table cannot support an estimate.
        """
        return all(getattr(self, f) is not None for f in BUS_FIELDS)

    def __post_init__(self):
        validate_spec_payload(self.as_payload())

    def as_payload(self) -> dict:
        """The plain dict form -- what an LLM would emit and what gets validated."""
        out = {
            "altitude_km": self.altitude_km,
            "inclination_deg": self.inclination_deg,
            "payload_class": self.payload_class,
            "power_w": self.power_w,
            "mission_duration_years": self.mission_duration_years,
            "label": self.label,
        }
        for k in ("epoch", "propellant_available_kg", "isp_s", "target_altitude_km",
                  *BUS_FIELDS, *ARRAY_FIELDS):
            v = getattr(self, k)
            if v is not None:
                out[k] = v
        return out

    def as_dict(self) -> dict:
        d = self.as_payload()
        d["provenance"] = dict(self.provenance)
        return d

    def epoch_datetime(self) -> datetime:
        """Epoch as a UTC datetime, defaulting to now.

        The default is a genuine default, not an invented value: a mission with
        no stated start date starts when you run it, and the compliance window
        is measured in years so the choice barely moves the answer.
        """
        if self.epoch is None:
            return datetime.now(timezone.utc)
        d = datetime.fromisoformat(self.epoch)
        return d if d.tzinfo else d.replace(tzinfo=timezone.utc)


def validate_spec_payload(payload: dict) -> None:
    """Reject anything that is not a legitimate requirement. Raises `SpecError`.

    Runs on every spec regardless of origin. The forbidden-field check comes
    first and is the reason this function exists -- a spec that smuggles in a
    mass is not a slightly-wrong spec, it is a fabricated result wearing the
    shape of an input.
    """
    if not isinstance(payload, dict):
        raise SpecError(f"spec must be an object, got {type(payload).__name__}")

    smuggled = sorted(set(payload) & set(FORBIDDEN_FIELDS))
    if smuggled:
        lines = "\n".join(
            f"  - {k}: {FORBIDDEN_FIELDS[k]}" for k in smuggled
        )
        raise SpecError(
            "spec contains quantities the engine computes; these are outputs, "
            f"not requirements:\n{lines}\n"
            "A spec states what is wanted, never what the answer is."
        )

    missing = [k for k in REQUIRED_FIELDS if payload.get(k) is None]
    if missing:
        raise SpecError(f"spec is missing required fields: {', '.join(missing)}")

    if payload["payload_class"] not in VALID_PAYLOAD_CLASSES:
        raise SpecError(
            f"unknown payload_class {payload['payload_class']!r}; "
            f"expected one of {sorted(VALID_PAYLOAD_CLASSES)}"
        )

    for key, (lo, hi) in BOUNDS.items():
        if payload.get(key) is None:
            continue
        v = payload[key]
        if not isinstance(v, (int, float)) or isinstance(v, bool):
            raise SpecError(f"{key} must be a number, got {v!r}")
        if not math.isfinite(v):
            raise SpecError(f"{key} must be finite, got {v!r}")
        if not lo <= v <= hi:
            raise SpecError(f"{key} = {v} is outside the plausible range {lo}-{hi}")

    target = payload.get("target_altitude_km")
    if target is not None and target <= payload["altitude_km"]:
        raise SpecError(
            f"target_altitude_km ({target}) must be above the insertion "
            f"altitude ({payload['altitude_km']}); lowering is disposal, "
            "not a mission"
        )

    given_bus = [f for f in BUS_FIELDS if payload.get(f) is not None]
    if given_bus and len(given_bus) != len(BUS_FIELDS):
        raise SpecError(
            f"partial bus geometry: got {sorted(given_bus)}, need all of "
            f"{list(BUS_FIELDS)}. A missing dimension cannot be filled with a "
            "typical value -- state all three or none, and the area will be "
            "refused rather than guessed."
        )
    given_array = [f for f in ARRAY_FIELDS if payload.get(f) is not None]
    if given_array and len(given_array) != len(ARRAY_FIELDS):
        raise SpecError(
            f"partial array geometry: got {sorted(given_array)}, need both of "
            f"{list(ARRAY_FIELDS)}."
        )

    if payload.get("epoch") is not None:
        try:
            datetime.fromisoformat(payload["epoch"])
        except (TypeError, ValueError) as exc:
            raise SpecError(f"epoch must be ISO-8601, got {payload['epoch']!r}") from exc


def spec_from_payload(payload: dict, provenance: dict | None = None) -> MissionSpec:
    """Build a `MissionSpec` from a plain dict, validating on the way in.

    This is the single entry point for anything a language model produces.
    """
    validate_spec_payload(payload)
    known = {
        "altitude_km", "inclination_deg", "payload_class", "power_w",
        "mission_duration_years", "epoch", "propellant_available_kg",
        "isp_s", "target_altitude_km", "label",
        *BUS_FIELDS, *ARRAY_FIELDS,
    }
    unknown = sorted(set(payload) - known - {"provenance"})
    if unknown:
        raise SpecError(
            f"spec contains unrecognised fields: {', '.join(unknown)}. "
            "The schema is closed -- an unknown field is either a typo or an "
            "attempt to pass the engine something it should be computing."
        )
    kwargs = {k: payload[k] for k in known if k in payload}
    return MissionSpec(provenance=dict(provenance or payload.get("provenance") or {}),
                       **kwargs)
