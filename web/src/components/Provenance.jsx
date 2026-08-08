// Rendering the four provenance kinds from out/studio.json.
//
// This file is where Gate 11 lives. The requirement is not that provenance
// exists in the payload -- sim/studio.py guarantees that -- but that a reader
// can SEE which numbers are solid and which are not, without hovering,
// expanding, or reading the docs.
//
// So the four kinds get four visibly different treatments:
//
//   stated     plain            the spec said so
//   computed   plain            derived exactly from stated inputs
//   estimated  marked + range   real uncertainty, interval always shown
//   refused    refusal block    the engine declined, and says why
//
// `refused` is the one that matters most. A refusal must never render as a
// blank cell, an em-dash, or a hedged number -- all three read as "small
// value" or "not important" rather than "this tool declined to answer".

const KIND_LABEL = {
  stated: "stated",
  computed: "computed",
  estimated: "estimated",
  refused: "refused",
};

export function KindBadge({ kind }) {
  return <span className={`kind-badge kind-${kind}`}>{KIND_LABEL[kind] ?? kind}</span>;
}

function fmt(value, unit) {
  if (value === null || value === undefined) return "—";
  if (typeof value === "string") return value;
  const abs = Math.abs(value);
  const digits = abs >= 100 ? 0 : abs >= 1 ? 2 : 4;
  return `${value.toFixed(digits)}${unit ? ` ${unit}` : ""}`;
}

// One row of the spec panel. The whole point of the component is that the
// `refused` branch is structurally different, not just styled differently --
// there is no number to render, so it renders the reason instead.
export function FieldRow({ name, field }) {
  if (!field) return null;
  const refused = field.kind === "refused";

  return (
    <div className={`field-row${refused ? " field-refused" : ""}`}>
      <div className="field-head">
        <span className="field-name">{name}</span>
        <KindBadge kind={field.kind} />
      </div>

      {refused ? (
        <div className="field-refusal">
          <div className="refusal-title">No single value — the engine declined</div>
          {field.interval && (
            <div className="field-value">
              {fmt(field.interval[0], "")}–{fmt(field.interval[1], field.unit)}
              <span className="field-range-note"> (bounded range, not an estimate)</span>
            </div>
          )}
          <p className="refusal-why">{field.detail}</p>
        </div>
      ) : (
        <>
          <div className="field-value">
            {fmt(field.value, field.unit)}
            {field.interval && (
              <span className="field-interval">
                {" "}
                [{fmt(field.interval[0], "")}–{fmt(field.interval[1], field.unit)}]
              </span>
            )}
          </div>
          {field.detail && <p className="field-detail">{field.detail}</p>}
        </>
      )}

      {field.sources?.length > 0 && <SourceList sources={field.sources} />}
    </div>
  );
}

// The named neighbours behind an interpolated mass. Required by Gate 11:
// "412 kg (estimated)" is a number with a disclaimer; naming the two real
// spacecraft it sits between is a number a reader can go and check.
function SourceList({ sources }) {
  return (
    <div className="source-list">
      <div className="source-label">Interpolated between</div>
      {sources.map((s) => (
        <div key={s.name} className="source-row">
          <span className="source-name">{s.name}</span>
          <span className="source-nums">
            {s.mass_kg} kg @ {s.power_w} W
            {s.kg_per_w != null && ` · ${s.kg_per_w.toFixed(3)} kg/W`}
          </span>
          <span className="source-cite">{s.source}</span>
        </div>
      ))}
    </div>
  );
}

// Verdicts. NOT_ASSESSABLE and AMBIGUOUS are refusals and must read as
// refusals -- `renderable` comes straight from sim/disposal.py, so the UI
// never has to decide which verdicts are real.
const VERDICT_TEXT = {
  COMPLIANT_NATURAL: "Complies — natural decay",
  COMPLIANT_WITH_DISPOSAL: "Complies — with a disposal burn",
  NON_COMPLIANT_INSUFFICIENT_PROPELLANT: "Does not comply — not enough propellant",
  NON_COMPLIANT_NO_SOLUTION: "Does not comply — no solution found",
  OUT_OF_SCOPE: "Out of scope for this rule",
  NOT_ASSESSABLE: "No verdict — the mass could not be resolved",
  AMBIGUOUS: "No verdict — compliance flips across the mass range",
};

export function Verdict({ compliance }) {
  if (!compliance) {
    return (
      <div className="verdict verdict-refused">
        <div className="verdict-line">No verdict — the design was not flown</div>
        <p className="refusal-why">
          Compliance needs a ram area, which needs a stated chassis. Nothing was
          assumed in its place.
        </p>
      </div>
    );
  }

  const ok = compliance.renderable;
  const v = compliance.verdict;
  return (
    <div className={`verdict ${ok ? "verdict-decided" : "verdict-refused"}`}>
      <div className="verdict-line">{VERDICT_TEXT[v] ?? v}</div>
      {!ok && (
        <div className="verdict-refusal-tag">
          This is a refusal, not a result. It is shown because the alternative —
          a blank, or a hedged guess — would read as an answer.
        </div>
      )}
      {compliance.notes?.slice(0, 3).map((n, i) => (
        <p key={i} className="refusal-why">{n}</p>
      ))}
      {ok && compliance.delta_v_interval_ms && (
        <div className="verdict-nums">
          <span>
            Δv {compliance.delta_v_interval_ms[0].toFixed(1)}–
            {compliance.delta_v_interval_ms[1].toFixed(1)} m/s
          </span>
          {compliance.propellant_interval_kg && (
            <span>
              propellant {compliance.propellant_interval_kg[0].toFixed(2)}–
              {compliance.propellant_interval_kg[1].toFixed(2)} kg
            </span>
          )}
        </div>
      )}
      {compliance.mass?.resolvable === false && (
        <p className="refusal-why">{compliance.mass.summary}</p>
      )}
    </div>
  );
}

// The solar band. The caveat text comes from sim/mission.py's
// SOLAR_BAND_CAVEAT so it cannot drift from the Python side, and it renders
// next to the numbers rather than in a footnote.
export function SolarBand({ band }) {
  if (!band) return null;
  const levels = ["low", "mean", "high"];
  return (
    <div className="solar-band">
      <div className="panel-title">Solar activity band</div>
      {/* The caveat text is SOLAR_BAND_CAVEAT from sim/mission.py, verbatim, so
          it cannot drift from the Python side. Its own first sentence is the
          headline -- bolding it here rather than prepending a copy, which is
          how this rendered the phrase twice. */}
      <div className="band-caveat-note">
        <strong>{band.caveat.split(". ")[0]}.</strong>{" "}
        {band.caveat.split(". ").slice(1).join(". ")}
      </div>
      <table className="data-table">
        <thead>
          <tr>
            <th>Level</th>
            <th>Natural decay</th>
            <th>Verdict</th>
          </tr>
        </thead>
        <tbody>
          {levels.map((l) => {
            const e = band.levels?.[l];
            if (!e) return null;
            const c = e.compliance;
            return (
              <tr key={l}>
                <td className="label">
                  {l} <span className="pct">p{band.percentiles[l]}</span>
                </td>
                <td>
                  {e.natural_decay_years
                    ? `${e.natural_decay_years[0]}–${e.natural_decay_years[1]} yr`
                    : "does not reenter in the window"}
                </td>
                <td className={c?.renderable ? "" : "fail"}>
                  {VERDICT_TEXT[c?.verdict] ?? "—"}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}
