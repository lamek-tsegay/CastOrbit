"""Studio payload. V2_BRIEF.md §7, Phase 11.

Gate 11 has three claims. Two of them are checkable here; the third
("verify in a real browser") is not something pytest can assert, and is done
with the Playwright harness instead.

  1. Every number traces to an engine field.
  2. Provenance is visible -- carried on every field, not implied.
  3. Refusals render as refusals -- guaranteed here by making a refusal a
     distinct, value-carrying state in the payload rather than a null.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim.atmosphere import SpaceWeather
from sim.spec import FORBIDDEN_FIELDS
from sim.studio import CD_BASELINE, GALLERY, build_studio_payload

DATA = Path(__file__).resolve().parent.parent / "data" / "SW-All.csv"
VALID_KINDS = {"stated", "computed", "estimated", "refused"}


@pytest.fixture(scope="module")
def payload():
    if not DATA.exists():
        pytest.skip(f"{DATA} not present")
    return build_studio_payload(SpaceWeather.load(DATA))


def test_payload_is_json_serialisable(payload):
    """It travels as a static asset; anything non-serialisable breaks silently."""
    text = json.dumps(payload)
    assert json.loads(text)["meta"]["cd_baseline"] == CD_BASELINE


def test_every_field_declares_its_provenance(payload):
    """**Gate 11:** provenance present on every field, and a recognised kind."""
    seen = 0
    for design in payload["designs"]:
        for name, f in design["fields"].items():
            assert f["kind"] in VALID_KINDS, f"{design['label']}.{name}: {f['kind']}"
            assert "detail" in f or f["kind"] == "computed", (
                f"{design['label']}.{name} has no explanation of where it came from"
            )
            seen += 1
    assert seen > 30, "too few fields to be a meaningful check"


def test_all_four_provenance_kinds_are_exercised(payload, capsys):
    """A gallery that never produces a refusal leaves the refusal path unrendered.

    The designs are chosen to hit every kind, so a regression that quietly
    stopped emitting refusals would fail here rather than just looking tidier.
    """
    kinds = {
        f["kind"]
        for d in payload["designs"]
        for f in d["fields"].values()
    }
    with capsys.disabled():
        print(f"\n  provenance kinds present: {sorted(kinds)}")
    assert kinds == VALID_KINDS, f"missing: {VALID_KINDS - kinds}"


def test_refusals_carry_a_reason_and_never_a_value(payload):
    """**Gate 11:** a refusal is a state with content, not a null.

    A refused field with `value: null` and no `detail` would render as a blank
    cell, which reads as "small" or "not applicable" rather than "this tool
    declined to answer".
    """
    refusals = [
        (d["label"], n, f)
        for d in payload["designs"]
        for n, f in d["fields"].items()
        if f["kind"] == "refused"
    ]
    assert refusals, "no refusals in the gallery -- the path is untested"
    for label, name, f in refusals:
        assert f["value"] is None, f"{label}.{name} is refused but carries a value"
        assert f.get("detail"), f"{label}.{name} is refused with no reason given"
        assert len(f["detail"]) > 40, f"{label}.{name} reason is too thin to be useful"


def test_estimated_fields_always_carry_an_interval(payload):
    """An estimate without its interval is the interval-vs-point failure again."""
    for d in payload["designs"]:
        for name, f in d["fields"].items():
            if f["kind"] == "estimated":
                assert f.get("interval"), f"{d['label']}.{name} estimated, no interval"
                lo, hi = f["interval"]
                assert lo <= f["value"] <= hi


def test_mass_names_its_two_neighbours(payload):
    """**Gate 11:** provenance visible, not just present.

    "1019 kg (estimated)" is a number with a disclaimer. Naming CryoSat-2 and
    Himawari-8, with their masses, powers and citations, is a number a reader
    can go and check.
    """
    checked = 0
    for d in payload["designs"]:
        m = d["fields"].get("dry_mass_kg")
        if not m:
            continue
        assert len(m["sources"]) == 2, f"{d['label']}: mass has no named neighbours"
        for s in m["sources"]:
            assert s["name"] and s["source"]
            assert s["mass_kg"] > 0 and s["power_w"] > 0
        checked += 1
    assert checked >= 4


def test_refused_mass_still_names_its_neighbours(payload):
    """Declining to interpolate does not mean declining to show the working."""
    d = next(x for x in payload["designs"] if x["label"] == "microsat-unresolvable-mass")
    m = d["fields"]["dry_mass_kg"]
    assert m["kind"] == "refused"
    assert m["value"] is None
    assert len(m["sources"]) == 2
    assert m["interval"][0] < m["interval"][1]


def test_no_spec_ever_carries_an_engine_computed_quantity(payload):
    """The Phase 10 boundary, still holding at the Phase 11 output."""
    for d in payload["designs"]:
        if d["spec"] is None:
            continue
        smuggled = set(d["spec"]) & set(FORBIDDEN_FIELDS)
        assert not smuggled, f"{d['label']} spec contains {smuggled}"


def test_unrenderable_verdicts_are_marked_unrenderable(payload):
    """**Gate 11:** refusals must be flagged so a UI cannot show them as results."""
    seen_refusal = False
    for d in payload["designs"]:
        c = d.get("compliance")
        if not c:
            continue
        assert "renderable" in c
        if c["verdict"] in {"NOT_ASSESSABLE", "AMBIGUOUS"}:
            assert c["renderable"] is False
            assert c["notes"], "a refusal with no notes gives the UI nothing to show"
            seen_refusal = True
        else:
            assert c["renderable"] is True
    assert seen_refusal, "gallery contains no non-renderable verdict"


def test_gallery_covers_both_refusal_verdicts(payload, capsys):
    verdicts = {
        d["label"]: (d.get("compliance") or {}).get("verdict")
        for d in payload["designs"]
    }
    with capsys.disabled():
        print("\n  verdicts:")
        for k, v in verdicts.items():
            print(f"    {k:30s} {v}")
    assert "NOT_ASSESSABLE" in verdicts.values()
    assert "AMBIGUOUS" in verdicts.values()


def test_solar_band_carries_the_caveat_verbatim(payload):
    """The caveat must come from Python, so the UI cannot drift from it."""
    from sim.mission import SOLAR_BAND_CAVEAT

    bands = [d["solar_band"] for d in payload["designs"] if d["solar_band"]]
    assert bands
    for b in bands:
        assert b["caveat"] == SOLAR_BAND_CAVEAT
        assert b["caveat"].startswith("Bounds, not scenarios")
        assert set(b["levels"]) == {"low", "mean", "high"}


def test_underspecified_design_asks_and_computes_nothing(payload):
    d = next(x for x in payload["designs"] if x["label"] == "underspecified")
    assert d["blocked"]
    assert d["spec"] is None
    assert d["fields"] == {}
    assert d["compliance"] is None
    assert d["orbit"] is None
    assert len(d["extraction"]["questions"]) == 5
    for q in d["extraction"]["questions"]:
        assert q["why"]


def test_missing_chassis_refuses_area_rather_than_assuming_one(payload):
    d = next(x for x in payload["designs"] if x["label"] == "no-chassis-stated")
    area = d["fields"]["ram_area_knife_edge_m2"]
    assert area["kind"] == "refused"
    assert "not assumed" in area["detail"]
    assert d["compliance"] is None, "cannot fly a design with no area"


def test_orbit_track_is_segmented_and_labelled_as_display_geometry(payload):
    for d in payload["designs"]:
        o = d.get("orbit")
        if not o:
            continue
        assert o["segments"] and all(len(s) > 1 for s in o["segments"])
        for seg in o["segments"]:
            for lat, lon in seg:
                assert -90.0 <= lat <= 90.0
                assert -180.0 <= lon <= 180.0
        # No segment may straddle the seam; that is what the splitting is for.
        for seg in o["segments"]:
            for (_, a), (_, b) in zip(seg, seg[1:]):
                assert abs(b - a) <= 180.0
        assert "outside the validated physics" in o["note"]


def test_gallery_designs_are_all_built(payload):
    assert len(payload["designs"]) == len(GALLERY)
    assert {d["label"] for d in payload["designs"]} == {g[0] for g in GALLERY}


def test_disposal_rule_is_cited_in_meta(payload):
    r = payload["meta"]["disposal_rule"]
    assert r["citation"] == "47 CFR § 25.283(e)"
    assert r["clock_starts"] == "end_of_mission"
    assert r["window_years"] == 5.0
