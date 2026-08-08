import { OUTCOME_COLORS, OUTCOME_LABELS } from "../lib/outcomeColors";

export default function OutcomeLegend({ counts }) {
  return (
    <div className="legend">
      {Object.entries(OUTCOME_COLORS).map(([key, color]) => (
        <span className="legend-item" key={key}>
          <span className="legend-swatch" style={{ background: color }} />
          {OUTCOME_LABELS[key]}
          {counts ? <span className="legend-count">{counts[key] ?? 0}</span> : null}
        </span>
      ))}
    </div>
  );
}
