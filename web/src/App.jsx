import { useMemo, useState } from "react";
import "./App.css";
import NavTabs from "./components/NavTabs";
import Scrubber from "./components/Scrubber";
import { DataProvider, useData } from "./context/DataContext";
import { PlaybackProvider } from "./context/PlaybackContext";
import AltitudeView from "./views/AltitudeView";
import GlobeView from "./views/GlobeView";
import SweepsView from "./views/SweepsView";
import ValidationView from "./views/ValidationView";

export default function App() {
  return (
    <DataProvider>
      <Shell />
    </DataProvider>
  );
}

function Shell() {
  const { status, batch, sweeps, error } = useData();
  const [view, setView] = useState("globe");

  const maxSeconds = useMemo(() => {
    if (!batch) return 0;
    const epoch = new Date(batch.meta.epoch);
    const end = new Date(batch.meta.window_end);
    return (end - epoch) / 1000;
  }, [batch]);

  if (status === "loading") {
    return <CenterMessage>Loading out/batch.json and out/sweeps.json…</CenterMessage>;
  }
  if (status === "error") {
    return (
      <CenterMessage>
        Failed to load data: {String(error)}.
        <br />
        Run <code>npm run sync-data</code> after generating{" "}
        <code>out/batch.json</code> and <code>out/sweeps.json</code>{" "}
        (<code>python -m sim.export</code>, <code>python -m sim.sweeps</code>).
      </CenterMessage>
    );
  }

  return (
    <PlaybackProvider epochIso={batch.meta.epoch} maxSeconds={maxSeconds}>
      <div className="app">
        <header className="app-header">
          <div className="app-title">
            <span className="app-name">CastOrbit</span>
            <span className="app-subtitle">
              {batch.meta.n_satellites} satellites · safe mode ·{" "}
              {batch.meta.epoch.slice(0, 16).replace("T", " ")} UT
            </span>
          </div>
          <NavTabs active={view} onChange={setView} />
        </header>

        <Scrubber />

        <main className="app-main">
          {view === "globe" && <GlobeView batch={batch} />}
          {view === "altitude" && <AltitudeView batch={batch} />}
          {view === "sweeps" && <SweepsView sweeps={sweeps} />}
          {view === "validation" && <ValidationView sweeps={sweeps} batch={batch} />}
        </main>
      </div>
    </PlaybackProvider>
  );
}

function CenterMessage({ children }) {
  return <div className="center-message">{children}</div>;
}
