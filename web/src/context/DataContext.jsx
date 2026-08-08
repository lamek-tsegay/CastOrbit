import { createContext, useContext, useEffect, useState } from "react";

// Loads out/batch.json and out/sweeps.json (synced into public/data/ by
// `npm run sync-data`) exactly once and hands them to every view via
// context. No view fetches its own copy or transforms the shape -- every
// number a component renders traces back to one of these two objects.
const DataContext = createContext(null);

export function DataProvider({ children }) {
  const [state, setState] = useState({
    status: "loading",
    batch: null,
    sweeps: null,
    error: null,
  });

  useEffect(() => {
    let cancelled = false;
    Promise.all([
      fetch(`${import.meta.env.BASE_URL}data/batch.json`).then((r) => {
        if (!r.ok) throw new Error(`batch.json: HTTP ${r.status}`);
        return r.json();
      }),
      fetch(`${import.meta.env.BASE_URL}data/sweeps.json`).then((r) => {
        if (!r.ok) throw new Error(`sweeps.json: HTTP ${r.status}`);
        return r.json();
      }),
    ])
      .then(([batch, sweeps]) => {
        if (!cancelled) setState({ status: "ready", batch, sweeps, error: null });
      })
      .catch((error) => {
        if (!cancelled) setState({ status: "error", batch: null, sweeps: null, error });
      });
    return () => {
      cancelled = true;
    };
  }, []);

  return <DataContext.Provider value={state}>{children}</DataContext.Provider>;
}

export function useData() {
  const ctx = useContext(DataContext);
  if (!ctx) throw new Error("useData must be used within a DataProvider");
  return ctx;
}
