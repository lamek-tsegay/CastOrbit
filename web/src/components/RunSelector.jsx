import { formatCd } from "../lib/format";

// Both Cd conventions in batch.json must stay visible and reachable -- the
// hard rule is that the UI never defaults to one silently. This renders both
// options with their own outcome counts always in view (never a plain
// dropdown hiding the alternative), and `selected` can start as `null` so a
// caller can force an explicit choice before rendering anything else.
export default function RunSelector({ runs, selected, onSelect, nSatellites }) {
  return (
    <div className="run-selector">
      {runs.map((run) => {
        const lost = run.outcome_counts.REENTERED ?? 0;
        const survived = nSatellites - lost;
        const isActive = selected === run.label;
        return (
          <button
            key={run.label}
            className={`run-selector-option${isActive ? " active" : ""}`}
            onClick={() => onSelect(run.label)}
          >
            <span className="run-selector-title">Cd = {formatCd(run.config.cd)}</span>
            <span className="run-selector-sub">{run.description}</span>
            <span className="run-selector-counts">
              {lost} lost / {survived} survived
            </span>
            <span className="run-selector-drag">
              Cd·A ∈ [{run.config.effective_drag_range_m2[0].toFixed(2)},{" "}
              {run.config.effective_drag_range_m2[1].toFixed(2)}] m²
            </span>
          </button>
        );
      })}
    </div>
  );
}
