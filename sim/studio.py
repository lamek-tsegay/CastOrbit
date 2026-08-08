"""Studio payload: the whole V2 pipeline, with provenance attached. Phase 11.

`ARCHITECTURE.md` §2 is unchanged — the frontend never simulates anything. This
module runs prose → spec → mass → geometry → flight → compliance in Python and
emits one JSON document the studio plays back.

**Provenance is the product here, not a decoration.** A studio that shows
"412 kg" next to "compliant" and lets a reader assume both are equally solid is
the failure `V2_BRIEF.md` §8 names. So every number leaves this module wrapped
in a `Field` carrying how it was obtained:

| kind | meaning | how the UI must show it |
|---|---|---|
| `stated` | the spec said so | plain |
| `computed` | derived exactly from stated inputs | plain |
| `estimated` | derived with real uncertainty; carries an interval | visibly marked, interval shown |
| `refused` | the engine declined; carries the reason | as a refusal, never as a blank or a hedge |

`refused` is a value-carrying state, not an error. A design whose mass cannot
be resolved still has a legitimate answer — "between 83 and 301 kg, and here is
why that cannot be narrowed" — and the studio's job is to render that as the
answer rather than hiding the row.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field as dc_field
from datetime import datetime, timezone
from pathlib import Path

from .atmosphere import SOLAR_ACTIVITY_PERCENTILES, ClimatologyDensity, SpaceWeather
from .bounded import Interval, assess_compliance_bounded, natural_decay_years_interval
from .disposal import default_rule
from .geometry import Box, Geometry, Panel
from .mass_model import MassEstimationError, estimate_dry_mass
from .prose import extract
from .spec import MissionSpec

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "out"

#: Drag convention. A modelling choice from data/satellite_specs.json, not a
#: spec input -- which is why `cd` is in FORBIDDEN_FIELDS.
CD_BASELINE = 2.2


@dataclass
class Field:
    """A number and how it was obtained. The unit of provenance in this payload."""

    value: float | str | None
    unit: str = ""
    kind: str = "computed"           # stated | computed | estimated | refused
    detail: str = ""                 # why, or where from
    interval: list[float] | None = None
    sources: list[dict] = dc_field(default_factory=list)

    def as_dict(self) -> dict:
        d = asdict(self)
        if not d["sources"]:
            d.pop("sources")
        if d["interval"] is None:
            d.pop("interval")
        if not d["detail"]:
            d.pop("detail")
        if not d["unit"]:
            d.pop("unit")
        return d


def _cd_field() -> Field:
    return Field(
        CD_BASELINE, "", "stated",
        "Free-molecular convention from data/satellite_specs.json "
        "(swept 2.0-2.4 elsewhere). Not a spec input.",
    )


def _mass_field(est) -> Field:
    """Mass, with its interval and the two spacecraft it came from.

    The neighbours are the point of this. "412 kg (estimated)" is a number with
    a disclaimer; "between CryoSat-2 at 684 kg and Himawari-8 at 1300 kg" is a
    number a reader can check.
    """
    sources = [
        {
            "name": n["name"],
            "mass_kg": n["mass_kg"],
            "power_w": n["power_w"],
            "kg_per_w": round(n["kg_per_w"], 4) if n["kg_per_w"] else None,
            "source": n["source"],
        }
        for n in est.neighbours
    ]
    if est.resolvable:
        return Field(
            round(est.mass_kg, 1), "kg", "estimated",
            est.provenance, [round(v, 1) for v in est.interval_kg], sources,
        )
    return Field(
        None, "kg", "refused",
        " ".join(est.refusal_reasons),
        [round(v, 1) for v in est.interval_kg], sources,
    )


def _geometry_fields(spec: MissionSpec) -> tuple[Field, Field, Geometry | None]:
    """Ram area from the stated chassis, or a refusal when none was stated."""
    if not spec.has_geometry:
        refusal = Field(
            None, "m^2", "refused",
            "No chassis dimensions were stated, so sim/geometry.py has nothing "
            "to solve. Ram area is not assumed from a typical bus: that would "
            "put a fabricated number into the drag term.",
        )
        return refusal, refusal, None

    parts = [Box(spec.bus_length_m, spec.bus_width_m, spec.bus_thickness_m)]
    if spec.solar_array_span_m and spec.solar_array_chord_m:
        parts.append(Panel(
            span_m=spec.solar_array_span_m,
            chord_m=spec.solar_array_chord_m,
            center_m=(spec.bus_length_m / 2 + spec.solar_array_span_m / 2, 0.0, 0.0),
        ))
    g = Geometry(parts, label=spec.label)
    knife, _ = g.knife_edge_area()
    broad, _ = g.broadside_area()
    detail = (
        f"sim/geometry.py projected-area solver over the stated "
        f"{spec.bus_length_m:g} x {spec.bus_width_m:g} x {spec.bus_thickness_m:g} m bus"
        + (f" plus a {spec.solar_array_span_m:g} x {spec.solar_array_chord_m:g} m array"
           if len(parts) > 1 else "")
    )
    return (
        Field(round(knife, 4), "m^2", "computed", detail + " (minimum-area attitude)"),
        Field(round(broad, 4), "m^2", "computed", detail + " (maximum-area attitude)"),
        g,
    )


def _solar_band_field(band: dict) -> dict:
    from .mission import SOLAR_BAND_CAVEAT
    return {
        "caveat": SOLAR_BAND_CAVEAT,
        "percentiles": SOLAR_ACTIVITY_PERCENTILES,
        "levels": band,
    }


def build_design(
    prose: str,
    sw: SpaceWeather,
    label: str,
    spec_override: MissionSpec | None = None,
) -> dict:
    """Run one design end to end and return the studio payload for it.

    `spec_override` supplies a spec directly for cases the deterministic prose
    parser cannot express yet (chassis dimensions in particular). The prose is
    still recorded and still parsed, so the console shows the real extraction
    result rather than a staged one.
    """
    extraction = extract(prose, label=label)
    spec = spec_override or extraction.spec

    payload: dict = {
        "label": label,
        "prose": prose,
        "extraction": extraction.as_dict(),
        "spec": None if spec is None else spec.as_dict(),
        "fields": {},
        "orbit": None,
        "compliance": None,
        "solar_band": None,
        "blocked": None,
    }

    if spec is None:
        payload["blocked"] = (
            "The description does not determine a spec. The engine was not run; "
            "these are the questions that stand between the prose and a design."
        )
        return payload

    fields: dict[str, Field] = {
        "altitude_km": Field(spec.altitude_km, "km", "stated", "from the description"),
        "inclination_deg": Field(
            spec.inclination_deg, "deg",
            spec.provenance.get("inclination_deg", "stated"),
            "sampled latitude for density (PHYSICS.md §10.2)",
        ),
        "power_w": Field(spec.power_w, "W", "stated", "from the description"),
        "mission_duration_years": Field(
            spec.mission_duration_years, "yr", "stated",
            "disposal clock starts at end of mission (47 CFR 25.283(e))",
        ),
        "payload_class": Field(spec.payload_class, "", "stated",
                               "selects the mass-model reference class"),
        "cd": _cd_field(),
    }

    try:
        est = estimate_dry_mass(spec.power_w, spec.payload_class)
    except MassEstimationError as exc:
        payload["blocked"] = f"Mass could not be estimated: {exc}"
        payload["fields"] = {k: v.as_dict() for k, v in fields.items()}
        return payload
    fields["dry_mass_kg"] = _mass_field(est)

    knife_f, broad_f, geom = _geometry_fields(spec)
    fields["ram_area_knife_edge_m2"] = knife_f
    fields["ram_area_broadside_m2"] = broad_f

    mass_iv = Interval(*est.interval_kg)
    if geom is not None:
        area = knife_f.value
        bc = Interval(mass_iv.lo / (CD_BASELINE * area), mass_iv.hi / (CD_BASELINE * area))
        fields["ballistic_coefficient_kg_m2"] = Field(
            None if not est.resolvable else round(est.mass_kg / (CD_BASELINE * area), 1),
            "kg/m^2",
            "estimated" if est.resolvable else "refused",
            "m/(Cd*A); inherits the mass interval (sim/bounded.py)",
            [round(bc.lo, 1), round(bc.hi, 1)],
        )
        fields["cd_times_area_m2"] = Field(
            round(CD_BASELINE * area, 4), "m^2", "computed",
            "the only separately observable combination (README central finding)",
        )
    else:
        fields["ballistic_coefficient_kg_m2"] = Field(
            None, "kg/m^2", "refused",
            "needs a ram area, which needs a stated chassis.",
        )

    # ---- flight and compliance, across the solar band -------------------
    if geom is not None:
        band: dict = {}
        for level in SOLAR_ACTIVITY_PERCENTILES:
            atmos = ClimatologyDensity.for_level(
                level, sw, lat_deg=spec.inclination_deg
            )
            decay, _ = natural_decay_years_interval(
                atmos, spec.altitude_km * 1e3, est, CD_BASELINE, knife_f.value
            )
            res = assess_compliance_bounded(
                atmos,
                spec.altitude_km * 1e3,
                est,
                cd=CD_BASELINE,
                area_m2=knife_f.value,
                isp_s=spec.isp_s or 1666.0,
                propellant_available_kg=spec.propellant_available_kg or 0.0,
                solar_activity=level,
                rule=default_rule(),
            )
            band[level] = {
                "natural_decay_years": None if decay is None else
                    [round(decay.lo, 2), round(decay.hi, 2)],
                "compliance": res.as_dict(),
            }
        payload["solar_band"] = _solar_band_field(band)
        payload["compliance"] = band["mean"]["compliance"]

    # ---- orbit track, for the centre panel ------------------------------
    payload["orbit"] = _orbit_track(spec)

    payload["fields"] = {k: v.as_dict() for k, v in fields.items()}
    return payload


def _orbit_track(spec: MissionSpec, n: int = 240) -> dict:
    """Lat/lon/altitude samples for the centre panel.

    Display geometry only, exactly as `sim/groundtrack.py` documents: the
    inclination and altitude are the spec's, the phase is illustrative, and
    lat/lon modelling is outside the validated physics (PHYSICS.md §1).
    """
    import math

    from .constants import MU, R_E

    a = R_E + spec.altitude_km * 1e3
    period_s = 2.0 * math.pi * math.sqrt(a ** 3 / MU)
    inc = math.radians(spec.inclination_deg)
    earth_rot = 360.9856 / 86400.0

    # Split at the +/-180 seam. A polyline that runs straight through the
    # wrap draws a line across the whole globe, which looks like an orbit
    # doing something it is not. Segmenting is display bookkeeping, done here
    # so the frontend stays a renderer (ARCHITECTURE.md §2).
    segments: list[list[list[float]]] = [[]]
    prev_lon = None
    for i in range(n):
        t = period_s * i / (n - 1)
        u = 2.0 * math.pi * t / period_s
        lat = math.degrees(math.asin(math.sin(inc) * math.sin(u)))
        raw = math.degrees(math.atan2(math.cos(inc) * math.sin(u), math.cos(u)))
        lon = ((raw - earth_rot * t + 180.0) % 360.0) - 180.0
        if prev_lon is not None and abs(lon - prev_lon) > 180.0:
            segments.append([])
        segments[-1].append([round(lat, 3), round(lon, 3)])
        prev_lon = lon

    return {
        "segments": [s for s in segments if len(s) > 1],
        "altitude_km": spec.altitude_km,
        "period_minutes": round(period_s / 60.0, 2),
        "note": (
            "Ground track is illustrative display geometry. Inclination and "
            "altitude are the spec's and the period is computed from them, but "
            "lat/lon modelling is outside the validated physics (PHYSICS.md §1)."
        ),
    }


# --------------------------------------------------------------------------
# The gallery
# --------------------------------------------------------------------------

def _spec(**kw) -> MissionSpec:
    base = dict(bus_length_m=3.0, bus_width_m=1.5, bus_thickness_m=0.3,
                solar_array_span_m=8.0, solar_array_chord_m=1.4,
                isp_s=1666.0, propellant_available_kg=40.0)
    base.update(kw)
    return MissionSpec(**base)


#: Designs the studio ships with. Chosen to exercise every provenance state --
#: a clean verdict, an unresolvable mass, an ambiguous verdict, a missing
#: chassis, and prose that does not determine a spec at all. A gallery of five
#: compliant designs would leave the refusal paths unrendered and untested.
GALLERY: list[tuple[str, str, MissionSpec | None]] = [
    (
        "sso-imager",
        "An optical imaging spacecraft in a 786 km orbit at 81.4 degrees "
        "inclination, drawing 1700 W, operating for 7 years.",
        _spec(altitude_km=786.0, inclination_deg=81.4,
              payload_class="earth_observation", power_w=1700.0,
              mission_duration_years=7.0, label="sso-imager"),
    ),
    (
        "leo-comsat",
        "A broadband communications satellite at 550 km, 53 degrees, 4 kW, "
        "5 year mission.",
        _spec(altitude_km=550.0, inclination_deg=53.0,
              payload_class="communications", power_w=4000.0,
              mission_duration_years=5.0, label="leo-comsat"),
    ),
    (
        "microsat-unresolvable-mass",
        "A small imaging microsatellite at 620 km, 80 degrees, 320 W, "
        "4 year mission.",
        _spec(altitude_km=620.0, inclination_deg=80.0,
              payload_class="earth_observation", power_w=320.0,
              mission_duration_years=4.0, label="microsat-unresolvable-mass"),
    ),
    (
        "knife-edge-ambiguous",
        "An earth observation platform at 390 km, 45 degrees, 1700 W, "
        "5 year mission.",
        _spec(altitude_km=390.0, inclination_deg=45.0,
              payload_class="earth_observation", power_w=1700.0,
              mission_duration_years=5.0, bus_length_m=3.4, bus_width_m=2.0,
              bus_thickness_m=0.6, solar_array_span_m=9.0,
              solar_array_chord_m=2.0, label="knife-edge-ambiguous"),
    ),
    (
        "no-chassis-stated",
        "A radar earth observation platform at 693 km, 65 degrees, 4.8 kW, "
        "seven year mission.",
        None,   # prose alone: no chassis, so ram area is refused
    ),
    (
        "underspecified",
        "Something like a Starlink but a bit bigger.",
        None,
    ),
]


def build_studio_payload(sw: SpaceWeather | None = None) -> dict:
    if sw is None:
        sw = SpaceWeather.load(DATA / "SW-All.csv")
    designs = []
    for label, prose, override in GALLERY:
        designs.append(build_design(prose, sw, label, spec_override=override))
    rule = default_rule()
    return {
        "meta": {
            "generated": datetime.now(timezone.utc).isoformat(),
            "cd_baseline": CD_BASELINE,
            "atmosphere_model": "NRLMSIS 2.1 via pymsis, climatology levels",
            "disposal_rule": {
                "id": rule.id, "label": rule.label, "citation": rule.citation,
                "window_years": rule.window_years, "clock_starts": rule.clock_starts,
                "source": rule.source,
            },
            "provenance_kinds": {
                "stated": "the spec said so",
                "computed": "derived exactly from stated inputs",
                "estimated": "derived with real uncertainty; carries an interval",
                "refused": "the engine declined; carries the reason",
            },
        },
        "designs": designs,
    }


def main() -> None:
    OUT.mkdir(exist_ok=True)
    payload = build_studio_payload()
    path = OUT / "studio.json"
    path.write_text(json.dumps(payload, indent=1))
    print(f"wrote {path}  ({path.stat().st_size / 1024:.0f} kB, "
          f"{len(payload['designs'])} designs)")
    for d in payload["designs"]:
        if d["blocked"]:
            print(f"  {d['label']:28s} BLOCKED: {d['blocked'][:60]}")
        else:
            v = (d["compliance"] or {}).get("verdict", "-")
            m = d["fields"]["dry_mass_kg"]
            print(f"  {d['label']:28s} mass={m['kind']:9s} verdict={v}")


if __name__ == "__main__":
    main()
