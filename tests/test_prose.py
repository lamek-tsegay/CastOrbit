"""Prose to spec. V2_BRIEF.md §7, Phase 10.

Gate 10 is two claims:
  1. Ten varied descriptions produce schema-valid specs.
  2. Underspecified input asks a question rather than inventing values.

The second matters more. A front door that fills gaps with plausible numbers
is worse than no front door, because the fabrication arrives wearing the shape
of a requirement and nothing downstream can tell.
"""

from __future__ import annotations

import json

import pytest

from sim.prose import (
    LLM_SYSTEM_PROMPT,
    extract,
    llm_prompt,
    spec_from_llm_response,
)
from sim.spec import (
    FORBIDDEN_FIELDS,
    REQUIRED_FIELDS,
    VALID_PAYLOAD_CLASSES,
    MissionSpec,
    SpecError,
    spec_from_payload,
    validate_spec_payload,
)

# --------------------------------------------------------------------------
# Ten varied descriptions -- Gate 10, part 1
#
# Varied on purpose: different phrasings for each field, different orderings,
# both word and numeric durations, both W and kW, definitional and explicit
# inclinations, and two that state a mass or area that must be discarded.
# --------------------------------------------------------------------------

TEN_DESCRIPTIONS: list[tuple[str, str]] = [
    ("starlink_like",
     "A broadband communications satellite at 550 km, 53 degrees, 4 kW, "
     "5 year mission."),
    ("sso_imager_explicit",
     "An optical imaging spacecraft in a 786 km orbit at 98.6 degrees "
     "inclination, drawing 1700 W, operating for 7 years."),
    ("polar_weather",
     "Polar weather monitoring satellite, altitude 820 km, 320 watts, "
     "five year mission."),
    ("iss_tech_demo",
     "A technology demonstration flying in the ISS orbit at 400 km with "
     "250 W available, planned for a two year mission."),
    ("equatorial_iot",
     "An equatorial IoT relay at 600 km, 150 watts, 3 year mission."),
    ("sar_heavy",
     "A radar earth observation platform at 693 km, 65 degrees, 4.8 kW, "
     "seven year mission."),
    ("nav_meo_edge",
     "A navigation satellite at 1900 km, 55 degrees, 2 kW, ten year mission."),
    ("cubesat_low",
     "A small multispectral imaging cubesat at 500 km, 40 degrees, 20 W, "
     "two year mission."),
    ("mass_stated_must_be_ignored",
     "A 400 kg communications satellite at 700 km, 45 degrees, 2.5 kW, "
     "6 year mission."),
    ("area_stated_must_be_ignored",
     "An earth observation satellite with 3.5 m2 of ram area at 620 km, "
     "80 degrees, 900 W, 4 year mission."),
]


@pytest.mark.parametrize("name,text", TEN_DESCRIPTIONS, ids=[n for n, _ in TEN_DESCRIPTIONS])
def test_gate_10_descriptions_produce_valid_specs(name, text, capsys):
    """**Phase 10 gate, part 1:** ten varied descriptions -> schema-valid specs."""
    r = extract(text, label=name)

    assert r.complete, (
        f"{name} did not produce a spec; questions: "
        f"{[q.field for q in r.questions]}"
    )
    # Validating again from the plain payload proves the spec would survive
    # a round trip through JSON, which is how it will actually travel.
    validate_spec_payload(r.spec.as_payload())
    spec_from_payload(json.loads(json.dumps(r.spec.as_payload())))

    with capsys.disabled():
        s = r.spec
        print(f"    {name:32s} {s.altitude_km:6.0f} km  {s.inclination_deg:5.1f} deg  "
              f"{s.power_w:6.0f} W  {s.mission_duration_years:4.1f} yr  "
              f"{s.payload_class}")


def test_gate_10_all_ten_are_distinct_and_cover_the_schema():
    """A gate passed by ten near-identical sentences would prove nothing."""
    specs = [extract(t, label=n).spec for n, t in TEN_DESCRIPTIONS]
    assert all(s is not None for s in specs)

    assert len({s.altitude_km for s in specs}) >= 8
    assert len({s.power_w for s in specs}) >= 8
    assert len({s.payload_class for s in specs}) >= 3
    assert len({s.mission_duration_years for s in specs}) >= 5
    # Both unit conventions and both inclination provenances are exercised.
    assert any(s.power_w >= 2000 for s in specs)   # kW phrasing
    assert any(s.power_w < 100 for s in specs)     # W phrasing
    provs = {p for s in specs for p in s.provenance.values()}
    assert provs <= {"stated", "definitional"}, (
        f"a value entered a spec by some route other than being stated: {provs}"
    )


# --------------------------------------------------------------------------
# Gate 10, part 2 -- underspecified input asks rather than invents
# --------------------------------------------------------------------------

UNDERSPECIFIED: list[tuple[str, str, str]] = [
    ("no_power", "An imaging satellite at 600 km, 45 degrees, 5 year mission.",
     "power_w"),
    ("no_altitude", "A communications satellite at 53 degrees, 3 kW, 5 years.",
     "altitude_km"),
    ("no_inclination", "An imaging satellite at 600 km, 800 W, 5 year mission.",
     "inclination_deg"),
    ("no_duration", "An imaging satellite at 600 km, 45 degrees, 800 W.",
     "mission_duration_years"),
    ("no_payload", "A satellite at 600 km, 45 degrees, 800 W, 5 year mission.",
     "payload_class"),
    ("sun_synchronous_only",
     "A sun-synchronous imaging satellite at 700 km, 1 kW, 5 year mission.",
     "inclination_deg"),
    ("vague", "Something like a Starlink but a bit bigger.", "altitude_km"),
    ("empty", "   ", "altitude_km"),
]


@pytest.mark.parametrize("name,text,expected_field", UNDERSPECIFIED,
                         ids=[n for n, _, _ in UNDERSPECIFIED])
def test_gate_10_underspecified_input_asks(name, text, expected_field, capsys):
    """**Phase 10 gate, part 2:** ask, never invent."""
    r = extract(text, label=name)

    assert not r.complete, f"{name} produced a spec from underspecified prose"
    assert r.spec is None
    fields = {q.field for q in r.questions}
    assert expected_field in fields, (
        f"{name}: expected a question about {expected_field}, got {fields}"
    )
    for q in r.questions:
        assert q.question.strip().endswith("?") or "cannot be used" in q.question
        assert q.why, f"question about {q.field} gives no reason"

    with capsys.disabled():
        print(f"    {name:22s} asks: {sorted(fields)}")


def test_sun_synchronous_is_asked_about_not_derived():
    """SSO inclination depends on altitude through J2, which is not modelled.

    Deriving it would mean using physics PHYSICS.md §1 excludes. Asking is the
    honest option, and the question says why.
    """
    r = extract("A sun-synchronous imaging satellite at 700 km, 1 kW, 5 years.")
    assert not r.complete
    q = next(q for q in r.questions if q.field == "inclination_deg")
    assert "J2" in q.why and "sun-synchronous" in q.why.lower()


def test_ambiguous_payload_class_is_a_question_not_a_coin_flip():
    """Two classes matched is underspecification, not a tie to break."""
    r = extract(
        "An imaging and communications relay at 600 km, 45 degrees, 1 kW, 5 years."
    )
    assert not r.complete
    assert "payload_class" in {q.field for q in r.questions}


def test_out_of_range_values_are_questioned_not_clamped():
    r = extract("An imaging satellite at 50 km, 45 degrees, 800 W, 5 year mission.")
    assert not r.complete
    assert any("cannot be used" in q.question for q in r.questions)


# --------------------------------------------------------------------------
# The boundary: no engine-computed quantity may enter a spec
# --------------------------------------------------------------------------

def test_stated_mass_is_discarded_with_a_note():
    """The prose says 400 kg. The engine computes mass. The 400 kg is dropped."""
    r = extract(
        "A 400 kg communications satellite at 700 km, 45 degrees, 2.5 kW, 6 years."
    )
    assert r.complete
    assert not hasattr(r.spec, "mass_kg")
    assert "mass_kg" not in r.spec.as_payload()
    assert any("mass" in d for d in r.discarded), r.discarded
    assert any("mass_model" in d for d in r.discarded), (
        "the note should point at what actually owns the number"
    )


def test_stated_area_is_discarded_with_a_note():
    r = extract(
        "An earth observation satellite with 3.5 m2 of ram area at 620 km, "
        "80 degrees, 900 W, 4 years."
    )
    assert r.complete
    assert "ram_area_m2" not in r.spec.as_payload()
    assert any("area" in d for d in r.discarded), r.discarded


@pytest.mark.parametrize("forbidden", sorted(FORBIDDEN_FIELDS))
def test_every_forbidden_field_is_rejected(forbidden):
    """Exhaustive over the list, so adding a field to it cannot be cosmetic."""
    payload = {
        "altitude_km": 600.0, "inclination_deg": 45.0,
        "payload_class": "earth_observation", "power_w": 800.0,
        "mission_duration_years": 5.0, forbidden: 1.0,
    }
    with pytest.raises(SpecError, match="quantities the engine computes"):
        validate_spec_payload(payload)


def test_rejection_names_the_component_that_owns_the_number():
    payload = {
        "altitude_km": 600.0, "inclination_deg": 45.0,
        "payload_class": "earth_observation", "power_w": 800.0,
        "mission_duration_years": 5.0, "dry_mass_kg": 240.0,
    }
    with pytest.raises(SpecError) as exc:
        validate_spec_payload(payload)
    assert "mass_model" in str(exc.value)


def test_schema_is_closed():
    """An unrecognised field is a typo or a smuggling attempt. Either way, no."""
    payload = {
        "altitude_km": 600.0, "inclination_deg": 45.0,
        "payload_class": "earth_observation", "power_w": 800.0,
        "mission_duration_years": 5.0, "solar_flux_assumption": 150.0,
    }
    with pytest.raises(SpecError, match="unrecognised fields"):
        spec_from_payload(payload)


def test_bad_payload_class_is_rejected():
    payload = {
        "altitude_km": 600.0, "inclination_deg": 45.0,
        "payload_class": "space_laser", "power_w": 800.0,
        "mission_duration_years": 5.0,
    }
    with pytest.raises(SpecError, match="unknown payload_class"):
        validate_spec_payload(payload)


def test_target_shell_must_be_above_insertion():
    with pytest.raises(SpecError, match="lowering is disposal"):
        MissionSpec(altitude_km=600.0, inclination_deg=45.0,
                    payload_class="earth_observation", power_w=800.0,
                    mission_duration_years=5.0, target_altitude_km=500.0)


# --------------------------------------------------------------------------
# The language-model door goes through the same boundary
# --------------------------------------------------------------------------

def test_llm_prompt_forbids_the_computed_quantities_by_name():
    p = llm_prompt("a satellite")
    assert p[0]["role"] == "system" and p[1]["content"] == "a satellite"
    for name in ("dry_mass_kg", "ram_area_m2", "ballistic_coefficient", "verdict"):
        assert name in LLM_SYSTEM_PROMPT, f"{name} not forbidden in the prompt"
    assert "do not guess" in LLM_SYSTEM_PROMPT.lower()
    for f in REQUIRED_FIELDS:
        assert f in LLM_SYSTEM_PROMPT
    for c in VALID_PAYLOAD_CLASSES:
        assert c in LLM_SYSTEM_PROMPT


def test_llm_output_is_validated_not_trusted():
    """A model that ignores the prompt must still not get a mass through.

    This is the load-bearing test of the whole phase: the guarantee cannot
    depend on the model complying.
    """
    rogue = json.dumps({
        "altitude_km": 700.0, "inclination_deg": 45.0,
        "payload_class": "communications", "power_w": 2500.0,
        "mission_duration_years": 6.0,
        "dry_mass_kg": 240.0,           # the model made this up
    })
    r = spec_from_llm_response(rogue)
    assert not r.complete
    assert r.spec is None
    assert any("dry_mass_kg" in d for d in r.discarded)


def test_llm_questions_are_passed_through():
    r = spec_from_llm_response(json.dumps({
        "questions": [{"field": "power_w", "question": "How much power?"}]
    }))
    assert not r.complete
    assert {q.field for q in r.questions} == {"power_w"}


def test_llm_malformed_output_is_a_question_not_a_crash():
    for bad in ("not json at all", "[1,2,3]", '{"altitude_km": '):
        r = spec_from_llm_response(bad)
        assert not r.complete
        assert r.questions


def test_llm_valid_output_produces_a_spec():
    r = spec_from_llm_response(json.dumps({
        "altitude_km": 705.0, "inclination_deg": 81.8,
        "payload_class": "earth_observation", "power_w": 1200.0,
        "mission_duration_years": 7.0, "label": "llm-spec",
    }))
    assert r.complete
    assert r.spec.altitude_km == 705.0
    assert r.spec.label == "llm-spec"


def test_both_doors_handle_retrograde_the_same_way():
    """A convention slip must not give two answers for the same mission.

    The prompt tells the model to emit `180 - i`, and a model will sometimes
    emit the raw inclination anyway. That is a convention slip, not a
    fabricated value, so both doors convert it and record the conversion --
    and they must agree, or the same satellite gets different physics
    depending on which door it came through.
    """
    from sim.prose import normalise_inclination

    converted, note = normalise_inclination(98.6)
    assert converted == pytest.approx(81.4)
    assert note and "retrograde" in note

    prose = extract(
        "An imaging satellite at 700 km, 98.6 degrees, 1 kW, 5 year mission."
    )
    llm = spec_from_llm_response(json.dumps({
        "altitude_km": 700.0, "inclination_deg": 98.6,
        "payload_class": "earth_observation", "power_w": 1000.0,
        "mission_duration_years": 5.0,
    }))
    assert prose.complete and llm.complete
    assert prose.spec.inclination_deg == pytest.approx(llm.spec.inclination_deg)
    assert prose.spec.inclination_deg == pytest.approx(81.4)
    # The LLM path records the conversion; nothing is silently rewritten.
    assert any("retrograde" in d for d in llm.discarded)


def test_inclination_above_180_is_rejected_not_wrapped():
    """Conversion applies to retrograde orbits, not to nonsense."""
    r = spec_from_llm_response(json.dumps({
        "altitude_km": 700.0, "inclination_deg": 260.0,
        "payload_class": "earth_observation", "power_w": 1000.0,
        "mission_duration_years": 5.0,
    }))
    assert not r.complete


# --------------------------------------------------------------------------
# The spec actually drives the engine
# --------------------------------------------------------------------------

def test_spec_serialises_and_round_trips():
    r = extract(TEN_DESCRIPTIONS[0][1], label="rt")
    payload = json.loads(json.dumps(r.spec.as_payload()))
    again = spec_from_payload(payload)
    assert again.as_payload() == r.spec.as_payload()


def test_extraction_result_serialises_for_a_ui():
    r = extract("An imaging satellite at 600 km, 45 degrees, 5 year mission.")
    d = r.as_dict()
    assert d["complete"] is False
    assert d["spec"] is None
    assert d["questions"] and all("why" in q for q in d["questions"])


def test_spec_feeds_a_mission_without_carrying_physical_quantities():
    """The spec's whole job: hand the engine requirements and nothing else.

    Mass and area arrive from the mass model and the geometry solver, not from
    the sentence that started this.
    """
    from sim.geometry import geometry_from_spec, ram_area_for_mode
    from sim.mass_model import estimate_dry_mass
    from sim.satellite import ThrusterMode
    import json as _json
    from pathlib import Path

    r = extract(TEN_DESCRIPTIONS[1][1], label="pipeline")
    assert r.complete
    spec = r.spec

    mass = estimate_dry_mass(spec.power_w, spec.payload_class)
    specs_json = _json.loads(
        (Path(__file__).resolve().parent.parent / "data" / "satellite_specs.json")
        .read_text()
    )
    area = ram_area_for_mode(geometry_from_spec(specs_json), ThrusterMode.SAFE_MODE)

    # Both came from the engine; neither appears anywhere in the spec.
    assert mass.tag == "estimated"
    assert area > 0
    assert not (set(spec.as_payload()) & set(FORBIDDEN_FIELDS))
