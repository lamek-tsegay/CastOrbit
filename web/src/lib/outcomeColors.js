// Outcome color mapping, per ARCHITECTURE.md §7. The same four colours,
// used identically in every view, are what let a user carry an association
// ("red = lost") from the globe to the altitude chart without re-reading a
// legend each time.
export const OUTCOME_COLORS = {
  REACHED_SHELL: "#2a9d3f",
  REENTERED: "#d1495b",
  PROPELLANT_EXHAUSTED: "#e8a33d",
  INDETERMINATE: "#8d99ae",
};

export const OUTCOME_LABELS = {
  REACHED_SHELL: "Reached shell",
  REENTERED: "Reentered",
  PROPELLANT_EXHAUSTED: "Propellant exhausted",
  INDETERMINATE: "Indeterminate",
};

export function outcomeColor(outcome) {
  return OUTCOME_COLORS[outcome] ?? "#8d99ae";
}
