// Three-panel studio. Phase 11.
//
// Layout only is borrowed from reference/orbital-forge: console left, 3D
// centre, spec panel right (360px | 1fr | 340px). Nothing else is. That
// reference computes mass in JavaScript from a made-up formula, which is
// precisely the failure V2_BRIEF.md §8 exists to prevent; here every number
// comes from out/studio.json and was computed by sim/studio.py.
//
// There is no physics in this file. It selects a precomputed design and
// renders it, per ARCHITECTURE.md §2.

import { useEffect, useMemo, useRef, useState } from "react";
import Globe from "react-globe.gl";
import { FieldRow, SolarBand, Verdict } from "../components/Provenance";

const EARTH_RADIUS_KM = 6371;

// The order fields appear in the spec panel: requirements first, then what the
// engine made of them. Reading top to bottom should walk the pipeline.
const FIELD_ORDER = [
  ["Altitude", "altitude_km"],
  ["Inclination", "inclination_deg"],
  ["Power", "power_w"],
  ["Mission duration", "mission_duration_years"],
  ["Payload class", "payload_class"],
  ["Drag coefficient", "cd"],
  ["Dry mass", "dry_mass_kg"],
  ["Ram area (knife-edge)", "ram_area_knife_edge_m2"],
  ["Ram area (broadside)", "ram_area_broadside_m2"],
  ["Cd·A", "cd_times_area_m2"],
  ["Ballistic coefficient", "ballistic_coefficient_kg_m2"],
];

export default function StudioView({ studio }) {
  const [selected, setSelected] = useState(studio.designs[0]?.label);
  const design = useMemo(
    () => studio.designs.find((d) => d.label === selected) ?? studio.designs[0],
    [studio, selected],
  );

  return (
    <div className="studio">
      <ConsolePanel
        designs={studio.designs}
        selected={design.label}
        onSelect={setSelected}
        design={design}
      />
      <OrbitPanel design={design} />
      <SpecPanel design={design} meta={studio.meta} />
    </div>
  );
}

// ---------------------------------------------------------------- console

function ConsolePanel({ designs, selected, onSelect, design }) {
  return (
    <section className="studio-panel studio-console">
      <div className="panel-title">Console</div>
      <p className="view-note">
        Descriptions are parsed to a spec in Python (<code>sim/prose.py</code>),
        then run through the engine. Select one to load its result.
      </p>

      <div className="design-list">
        {designs.map((d) => (
          <button
            key={d.label}
            className={`design-item${d.label === selected ? " active" : ""}`}
            onClick={() => onSelect(d.label)}
          >
            <span className="design-name">{d.label}</span>
            <StatusChip design={d} />
          </button>
        ))}
      </div>

      <div className="console-thread">
        <div className="console-msg console-user">{design.prose}</div>
        <ExtractionReply design={design} />
      </div>
    </section>
  );
}

function StatusChip({ design }) {
  if (design.blocked) return <span className="chip chip-ask">asks</span>;
  const kind = design.fields?.dry_mass_kg?.kind;
  const v = design.compliance?.verdict;
  if (kind === "refused" || (v && !design.compliance.renderable)) {
    return <span className="chip chip-refused">refusal</span>;
  }
  return <span className="chip chip-ok">verdict</span>;
}

function ExtractionReply({ design }) {
  const ex = design.extraction;
  if (design.blocked) {
    return (
      <div className="console-msg console-system console-ask">
        <div className="refusal-title">Not enough to build a spec</div>
        <p className="refusal-why">{design.blocked}</p>
        <ul className="question-list">
          {ex.questions.map((q) => (
            <li key={q.field}>
              <span className="q-field">{q.field}</span>
              <span className="q-text">{q.question}</span>
              <span className="q-why">{q.why}</span>
            </li>
          ))}
        </ul>
      </div>
    );
  }
  return (
    <div className="console-msg console-system">
      <div>Spec accepted. The engine computed the panel on the right.</div>
      {ex.discarded?.length > 0 && (
        <ul className="discard-list">
          {ex.discarded.map((d, i) => (
            <li key={i}>{d}</li>
          ))}
        </ul>
      )}
      {ex.questions?.length > 0 && (
        <p className="refusal-why">
          The prose alone left {ex.questions.length} field
          {ex.questions.length > 1 ? "s" : ""} undetermined; this design supplies
          them explicitly.
        </p>
      )}
    </div>
  );
}

// ------------------------------------------------------------------ orbit

function OrbitPanel({ design }) {
  const wrapRef = useRef(null);
  const [dims, setDims] = useState({ width: 600, height: 520 });

  useEffect(() => {
    if (!wrapRef.current) return undefined;
    const ro = new ResizeObserver((entries) => {
      const { width, height } = entries[0].contentRect;
      if (width > 0 && height > 0) setDims({ width, height });
    });
    ro.observe(wrapRef.current);
    return () => ro.disconnect();
  }, []);

  // The track is drawn as paths, not points. Extruded points at an orbital
  // altitude render as tall bars, which reads as a spike rather than a
  // trajectory. Segments arrive pre-split at the +/-180 seam from Python.
  const paths = useMemo(() => {
    const o = design.orbit;
    if (!o?.segments) return [];
    const alt = o.altitude_km / EARTH_RADIUS_KM;
    return o.segments.map((seg) => seg.map(([lat, lng]) => [lat, lng, alt]));
  }, [design]);

  return (
    <section className="studio-panel studio-orbit" ref={wrapRef}>
      {design.orbit ? (
        <>
          <Globe
            width={dims.width}
            height={dims.height}
            backgroundColor="#000000"
            globeImageUrl="/textures/earth-dark.jpg"
            showAtmosphere
            atmosphereColor="#4fa3ff"
            atmosphereAltitude={0.12}
            pathsData={paths}
            pathPointLat={(p) => p[0]}
            pathPointLng={(p) => p[1]}
            pathPointAlt={(p) => p[2]}
            pathColor={() => "#4fa3ff"}
            pathStroke={1.6}
            pathTransitionDuration={0}
          />
          <div className="orbit-overlay">
            <div className="orbit-stats">
              <Stat label="Altitude" value={`${design.orbit.altitude_km} km`} />
              <Stat label="Period" value={`${design.orbit.period_minutes} min`} />
              <Stat
                label="Inclination"
                value={`${design.fields.inclination_deg?.value}°`}
              />
            </div>
            <p className="orbit-note">{design.orbit.note}</p>
          </div>
        </>
      ) : (
        <div className="orbit-empty">
          <div className="refusal-title">No orbit to draw</div>
          <p className="refusal-why">
            {design.blocked ??
              "The design was not flown, so there is no trajectory to show."}
          </p>
        </div>
      )}
    </section>
  );
}

function Stat({ label, value }) {
  return (
    <div className="orbit-stat">
      <div className="orbit-stat-label">{label}</div>
      <div className="orbit-stat-value">{value}</div>
    </div>
  );
}

// ------------------------------------------------------------------- spec

function SpecPanel({ design, meta }) {
  return (
    <section className="studio-panel studio-spec">
      <div className="panel-title">Design</div>

      {design.blocked ? (
        <div className="field-refusal">
          <div className="refusal-title">Nothing computed</div>
          <p className="refusal-why">{design.blocked}</p>
        </div>
      ) : (
        <>
          <Verdict compliance={design.compliance} />

          <div className="panel-title" style={{ marginTop: 18 }}>
            Fields
          </div>
          <p className="view-note" style={{ marginBottom: 10 }}>
            Every value below is a field of <code>out/studio.json</code>,
            computed by <code>sim/studio.py</code>. Badges are the provenance
            kinds in <code>meta.provenance_kinds</code>.
          </p>
          {FIELD_ORDER.map(([label, key]) => (
            <FieldRow key={key} name={label} field={design.fields[key]} />
          ))}

          <SolarBand band={design.solar_band} />

          <div className="rule-cite">
            Rule applied: {meta.disposal_rule.label} ({meta.disposal_rule.citation}),{" "}
            {meta.disposal_rule.window_years} yr from{" "}
            {meta.disposal_rule.clock_starts.replace(/_/g, " ")}.
          </div>
        </>
      )}
    </section>
  );
}
