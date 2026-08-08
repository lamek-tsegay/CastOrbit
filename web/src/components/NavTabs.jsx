const VIEWS = [
  { id: "studio", label: "Studio" },
  { id: "globe", label: "Globe" },
  { id: "altitude", label: "Altitude" },
  { id: "sweeps", label: "Sweeps" },
  { id: "validation", label: "Validation" },
];

export default function NavTabs({ active, onChange }) {
  return (
    <nav className="nav-tabs">
      {VIEWS.map((v) => (
        <button
          key={v.id}
          className={`nav-tab${active === v.id ? " active" : ""}`}
          onClick={() => onChange(v.id)}
        >
          {v.label}
        </button>
      ))}
    </nav>
  );
}
