import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState } from "react";
import { isoPlusSeconds } from "../lib/format";

// The shared scrubber (ARCHITECTURE.md §7: "Playback scrubber with a UT
// clock, shared across views"). One clock, one current time, read by
// whichever view cares (Globe and Altitude animate against it; Sweeps and
// Validation are static and simply ignore it).
const PlaybackContext = createContext(null);

export function PlaybackProvider({ epochIso, maxSeconds, children }) {
  const [currentT, setCurrentT] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(1800); // simulated seconds per real second
  const rafRef = useRef(null);
  const lastRef = useRef(null);

  useEffect(() => {
    if (!playing) {
      lastRef.current = null;
      return undefined;
    }
    const tick = (now) => {
      if (lastRef.current == null) lastRef.current = now;
      const realDeltaS = (now - lastRef.current) / 1000;
      lastRef.current = now;
      setCurrentT((prev) => {
        const next = prev + realDeltaS * speed;
        if (next >= maxSeconds) {
          setPlaying(false);
          return maxSeconds;
        }
        return next;
      });
      rafRef.current = requestAnimationFrame(tick);
    };
    rafRef.current = requestAnimationFrame(tick);
    return () => {
      if (rafRef.current) cancelAnimationFrame(rafRef.current);
    };
  }, [playing, speed, maxSeconds]);

  const seek = useCallback(
    (t) => setCurrentT(Math.max(0, Math.min(maxSeconds, t))),
    [maxSeconds],
  );
  const togglePlay = useCallback(() => {
    setPlaying((p) => {
      if (!p && currentT >= maxSeconds) setCurrentT(0); // replay from the start
      return !p;
    });
  }, [currentT, maxSeconds]);

  const utClock = useMemo(() => isoPlusSeconds(epochIso, currentT), [epochIso, currentT]);

  const value = useMemo(
    () => ({
      currentT, maxSeconds, playing, speed, utClock, epochIso,
      seek, togglePlay, setSpeed,
    }),
    [currentT, maxSeconds, playing, speed, utClock, epochIso, seek, togglePlay],
  );

  return <PlaybackContext.Provider value={value}>{children}</PlaybackContext.Provider>;
}

export function usePlayback() {
  const ctx = useContext(PlaybackContext);
  if (!ctx) throw new Error("usePlayback must be used within a PlaybackProvider");
  return ctx;
}
