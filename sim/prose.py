"""Prose to spec. V2_BRIEF.md §7, Phase 10 -- the smallest phase in the project.

Turns "a 400 kg imaging satellite at 600 km, 5 year mission, 1.2 kW" into a
`MissionSpec` the sizing code consumes -- or into questions, when the prose
does not determine one.

---

**The rule, stated once.** This module extracts *requirements*. It never
produces a quantity the engine computes. If prose says "a 400 kg satellite",
the 400 kg is **discarded with a note**, not passed through: dry mass is
`sim/mass_model.py`'s output, and accepting a stated mass would let a sentence
overwrite a computation. The same goes for ram area, ballistic coefficient,
decay time and verdicts (`sim/spec.FORBIDDEN_FIELDS`).

That is deliberately more aggressive than it needs to be for a well-behaved
user, and the reason is the badly-behaved case. Once a language model is in the
loop, "the prose said so" and "the model made it up" are indistinguishable
downstream.

---

**Two front doors, one boundary.**

`extract()` is a deliberately conservative deterministic parser. It reads
explicit numbers with units and a small set of definitional orbit names, and
asks about everything else. It exists because the guarantee this phase makes
should not depend on a network call, and because it is testable offline.

`llm_prompt()` / `spec_from_llm_response()` are the language-model path. The
prompt constrains the model to the same schema; the response goes through
`sim.spec.validate_spec_payload` exactly as the deterministic path does. **No
part of the safety argument rests on the model complying** -- if it emits a
mass, validation rejects the spec.

Neither door can widen the schema. That is the whole design.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field

from .spec import (
    FORBIDDEN_FIELDS,
    NAMED_INCLINATIONS_DEG,
    REQUIRED_FIELDS,
    MissionSpec,
    Question,
    SpecError,
    spec_from_payload,
)


@dataclass
class ExtractionResult:
    """A spec, or the questions standing between the prose and one."""

    spec: MissionSpec | None
    questions: list[Question] = field(default_factory=list)
    discarded: list[str] = field(default_factory=list)
    raw: dict = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return self.spec is not None

    def as_dict(self) -> dict:
        return {
            "complete": self.complete,
            "spec": None if self.spec is None else self.spec.as_dict(),
            "questions": [q.as_dict() for q in self.questions],
            "discarded": self.discarded,
        }


# --------------------------------------------------------------------------
# Unit-aware number reading
# --------------------------------------------------------------------------

_NUM = r"(\d+(?:\.\d+)?)"

_WORD_NUMBERS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
    "seven": 7, "eight": 8, "nine": 10 - 1, "ten": 10, "twelve": 12,
    "fifteen": 15, "twenty": 20, "twentyfive": 25,
}


def _search_float(pattern: str, text: str, scale: float = 1.0) -> float | None:
    m = re.search(pattern, text, re.I)
    return float(m.group(1)) * scale if m else None


def _find_altitude_km(text: str) -> float | None:
    for pat, scale in (
        (rf"{_NUM}\s*(?:km|kilometres|kilometers)\b", 1.0),
        (rf"{_NUM}\s*(?:nmi|nautical miles)\b", 1.852),
    ):
        # "5 year mission at 600 km" -- require the number to read as an
        # altitude, not as any other km-quantity in the sentence.
        for m in re.finditer(pat, text, re.I):
            window = text[max(0, m.start() - 40):m.end() + 20].lower()
            if any(w in window for w in ("altitude", "orbit", "at ", "circular", "shell")):
                return float(m.group(1)) * scale
        v = _search_float(pat, text, scale)
        if v is not None:
            return v
    return None


def _find_power_w(text: str) -> float | None:
    v = _search_float(rf"{_NUM}\s*(?:kw|kilowatts?)\b", text, 1000.0)
    if v is not None:
        return v
    return _search_float(rf"{_NUM}\s*(?:w|watts?)\b(?!\w)", text)


def _find_years(text: str) -> float | None:
    v = _search_float(rf"{_NUM}\s*(?:-|\s)?\s*year", text)
    if v is not None:
        return v
    m = re.search(r"\b(" + "|".join(_WORD_NUMBERS) + r")[- ]year", text, re.I)
    return float(_WORD_NUMBERS[m.group(1).lower()]) if m else None


def normalise_inclination(deg: float) -> tuple[float, str | None]:
    """Express an inclination in the engine's 0-90 sampling-latitude convention.

    A retrograde orbit at 98.6 deg samples density at the same latitudes as a
    prograde one at 81.4 deg, and latitude is all this engine takes from
    inclination (PHYSICS.md §10.2). The conversion is therefore lossless *for
    this model*, and it is a change of convention rather than a guess -- but it
    is still a transformation of the user's input, so it returns a note and the
    caller records it.

    Shared by both front doors deliberately. A conversion applied on one path
    and not the other is the kind of inconsistency that produces two different
    answers for the same mission.
    """
    if 90.0 < deg <= 180.0:
        return 180.0 - deg, (
            f"read inclination {deg:g} deg as retrograde and converted to the "
            f"engine's sampling-latitude convention, {180.0 - deg:g} deg "
            "(PHYSICS.md §10.2 samples density at one latitude, so the two are "
            "equivalent here)"
        )
    return deg, None


def _find_inclination(text: str) -> tuple[float | None, str | None]:
    """Returns `(degrees, provenance)`; provenance is 'stated' or 'definitional'."""
    m = re.search(rf"{_NUM}\s*(?:deg|degree|°)", text, re.I)
    if m:
        v, _ = normalise_inclination(float(m.group(1)))
        return v, "stated"
    low = text.lower()
    for name, deg in NAMED_INCLINATIONS_DEG.items():
        if name in low:
            return deg, "definitional"
    if "iss" in low.split() or "space station" in low:
        return 51.6, "definitional"
    return None, None


_PAYLOAD_KEYWORDS: list[tuple[str, tuple[str, ...]]] = [
    ("earth_observation", (
        "imaging", "imagery", "earth observation", "observation", "remote sensing",
        "sar", "radar", "weather", "meteorolog", "optical", "multispectral",
        "hyperspectral", "camera", "surveillance", "monitoring",
    )),
    ("communications", (
        "communication", "comsat", "broadband", "internet", "relay", "telecom",
        "backhaul", "connectivity", "data relay", "iot",
    )),
    ("navigation", ("navigation", "gnss", "positioning", "timing signal")),
    ("technology", (
        "technology demonstration", "tech demo", "demonstrator", "experimental",
        "pathfinder", "in-orbit demonstration",
    )),
]


def _find_payload_class(text: str) -> str | None:
    low = text.lower()
    hits = [cls for cls, words in _PAYLOAD_KEYWORDS if any(w in low for w in words)]
    # Ambiguity is underspecification, not a tie to be broken.
    return hits[0] if len(hits) == 1 else None


#: Phrases that name a quantity the engine computes. Matching one means the
#: prose tried to assert an output; the value is dropped and recorded.
_FORBIDDEN_PHRASES: list[tuple[str, str]] = [
    (rf"{_NUM}\s*(?:kg|kilograms?)\b", "mass"),
    (rf"{_NUM}\s*(?:m2|m\^2|square metres|square meters|m²)", "area"),
    (r"ballistic coefficient", "ballistic coefficient"),
    (r"\bcd\s*=\s*\d", "drag coefficient"),
    (r"drag coefficient", "drag coefficient"),
]


def _find_discards(text: str) -> list[str]:
    out = []
    for pat, what in _FORBIDDEN_PHRASES:
        m = re.search(pat, text, re.I)
        if m:
            owner = {
                "mass": FORBIDDEN_FIELDS["dry_mass_kg"],
                "area": FORBIDDEN_FIELDS["ram_area_m2"],
                "ballistic coefficient": FORBIDDEN_FIELDS["ballistic_coefficient"],
                "drag coefficient": FORBIDDEN_FIELDS["cd"],
            }[what]
            out.append(
                f"ignored a stated {what} ({m.group(0).strip()!r}): {owner}"
            )
    return out


# --------------------------------------------------------------------------
# Questions
# --------------------------------------------------------------------------

_QUESTIONS: dict[str, tuple[str, str]] = {
    "altitude_km": (
        "What circular altitude should the satellite operate at, in km?",
        "Altitude drives atmospheric density over about four orders of "
        "magnitude across LEO; nothing downstream is meaningful without it.",
    ),
    "inclination_deg": (
        "What orbital inclination, in degrees?",
        "Inclination sets the latitude at which density is sampled "
        "(PHYSICS.md §10.2). Sun-synchronous is not enough on its own: its "
        "inclination depends on altitude through J2, which this engine does "
        "not model (PHYSICS.md §1).",
    ),
    "payload_class": (
        "What is the payload -- earth observation, communications, navigation, "
        "or a technology demonstration?",
        "Communications satellites are roughly 3x lighter per watt than "
        "instrument-carrying ones, so the class picks which reference "
        "spacecraft the mass is interpolated between.",
    ),
    "power_w": (
        "What is the spacecraft's electrical power, in watts?",
        "Power is the predictor the mass model interpolates on "
        "(sim/mass_model.py). Without it there is no mass, and without a mass "
        "there is no ballistic coefficient and no verdict.",
    ),
    "mission_duration_years": (
        "How long is the operational mission, in years?",
        "The disposal clock starts at end of mission (47 CFR 25.283(e)), so "
        "the duration sets when the compliance window opens.",
    ),
}


def _question_for(field_name: str) -> Question:
    q, why = _QUESTIONS[field_name]
    return Question(field=field_name, question=q, why=why)


# --------------------------------------------------------------------------
# The deterministic door
# --------------------------------------------------------------------------

def extract(text: str, label: str = "from-prose") -> ExtractionResult:
    """Read requirements out of prose. Conservative by design.

    Only explicit values and definitional orbit names are accepted. Anything
    else produces a `Question`. There is no "reasonable default" path, because
    a reasonable default is exactly what a fabricated number looks like.
    """
    if not text or not text.strip():
        return ExtractionResult(
            spec=None,
            questions=[_question_for(f) for f in REQUIRED_FIELDS],
            discarded=[],
        )

    provenance: dict[str, str] = {}
    payload: dict = {"label": label}

    alt = _find_altitude_km(text)
    if alt is not None:
        payload["altitude_km"] = alt
        provenance["altitude_km"] = "stated"

    inc, inc_prov = _find_inclination(text)
    if inc is not None:
        payload["inclination_deg"] = inc
        provenance["inclination_deg"] = inc_prov

    cls = _find_payload_class(text)
    if cls is not None:
        payload["payload_class"] = cls
        provenance["payload_class"] = "stated"

    pw = _find_power_w(text)
    if pw is not None:
        payload["power_w"] = pw
        provenance["power_w"] = "stated"

    yrs = _find_years(text)
    if yrs is not None:
        payload["mission_duration_years"] = yrs
        provenance["mission_duration_years"] = "stated"

    discarded = _find_discards(text)
    questions = [_question_for(f) for f in REQUIRED_FIELDS if f not in payload]

    if questions:
        return ExtractionResult(None, questions, discarded, raw=payload)

    try:
        spec = spec_from_payload(payload, provenance)
    except SpecError as exc:
        # A value was read but is not usable -- out of range, say. That is
        # still a question, not a silent correction.
        return ExtractionResult(
            None,
            [Question(field="spec", question=f"That spec cannot be used: {exc}",
                      why="Values are range-checked before the engine sees them.")],
            discarded,
            raw=payload,
        )
    return ExtractionResult(spec, [], discarded, raw=payload)


# --------------------------------------------------------------------------
# The language-model door
# --------------------------------------------------------------------------

LLM_SYSTEM_PROMPT = """\
You convert a spacecraft mission description into a JSON specification.

Emit ONLY a JSON object. No prose, no code fences.

You may emit exactly these fields:
  altitude_km               number   circular operating altitude, 150-2000
  inclination_deg           number   0-90 (for retrograde, emit 180 - i)
  payload_class             string   one of: earth_observation, communications,
                                     navigation, technology
  power_w                   number   spacecraft electrical power, 1-30000
  mission_duration_years    number   operational life before disposal, 0.01-30
  epoch                     string   optional, ISO-8601
  propellant_available_kg   number   optional, only if the description states it
  isp_s                     number   optional, only if the description states it
  target_altitude_km        number   optional, only if raising to a shell
  label                     string   optional short name

You MUST NOT emit any of the following. They are computed by the simulator,
and asserting one would replace a calculation with a guess:
  mass, dry_mass_kg, wet_mass_kg, area_m2, ram_area_m2, cd, drag_coefficient,
  ballistic_coefficient, decay_time, delta_v_ms, propellant_required_kg,
  compliance, verdict, density

If the description states a mass or an area, IGNORE IT. Those are outputs.

If any REQUIRED field (altitude_km, inclination_deg, payload_class, power_w,
mission_duration_years) is not determined by the description, do not guess.
Instead emit:
  {"questions": [{"field": "<name>", "question": "<what to ask>"}]}

Do not infer inclination from "sun-synchronous": it depends on altitude
through J2, which this simulator does not model. Ask instead.
"""


def llm_prompt(text: str) -> list[dict]:
    """Messages for the constrained-JSON call. No network access happens here."""
    return [
        {"role": "system", "content": LLM_SYSTEM_PROMPT},
        {"role": "user", "content": text},
    ]


def spec_from_llm_response(
    response: str | dict, label: str = "from-llm"
) -> ExtractionResult:
    """Validate whatever the model returned. Trusts nothing.

    Handles the three things a model actually does: return a spec, return
    questions, or return something malformed. All three are ordinary outcomes
    here; only the first produces a spec, and only after full validation.
    """
    if isinstance(response, str):
        try:
            payload = json.loads(response)
        except json.JSONDecodeError as exc:
            return ExtractionResult(
                None,
                [Question("response", f"The model did not return valid JSON: {exc}",
                          "The schema is enforced on this side of the call.")],
            )
    else:
        payload = response

    if not isinstance(payload, dict):
        return ExtractionResult(
            None,
            [Question("response", "The model returned a non-object payload.",
                      "The schema is enforced on this side of the call.")],
        )

    if "questions" in payload:
        qs = []
        for q in payload.get("questions") or []:
            if isinstance(q, dict) and q.get("field") in _QUESTIONS:
                qs.append(_question_for(q["field"]))
            elif isinstance(q, dict):
                qs.append(Question(
                    field=str(q.get("field", "unknown")),
                    question=str(q.get("question", "Please clarify.")),
                    why="The description did not determine this.",
                ))
        return ExtractionResult(None, qs or [_question_for(f) for f in REQUIRED_FIELDS])

    payload.setdefault("label", label)
    discarded = [
        f"model emitted {k!r}, which the engine computes: {FORBIDDEN_FIELDS[k]}"
        for k in sorted(set(payload) & set(FORBIDDEN_FIELDS))
    ]

    # Same convention handling as the deterministic door. A model told "emit
    # 180 - i" will sometimes emit the raw inclination anyway; that is a
    # convention slip, not a fabricated value, so it is converted and noted
    # rather than rejected.
    inc = payload.get("inclination_deg")
    if isinstance(inc, (int, float)) and not isinstance(inc, bool):
        converted, note = normalise_inclination(float(inc))
        if note:
            payload["inclination_deg"] = converted
            discarded.append(note)
    try:
        spec = spec_from_payload(payload, {k: "stated" for k in payload})
    except SpecError as exc:
        return ExtractionResult(
            None,
            [Question("spec", f"The model's spec was rejected: {exc}",
                      "Every spec is validated regardless of what produced it.")],
            discarded,
        )
    return ExtractionResult(spec, [], discarded)
