import { usePlayback } from "../context/PlaybackContext";
import { formatHours, formatUTC } from "../lib/format";

const SPEEDS = [
  { label: "15 min/s", value: 900 },
  { label: "30 min/s", value: 1800 },
  { label: "1 h/s", value: 3600 },
  { label: "4 h/s", value: 14400 },
];

export default function Scrubber() {
  const { currentT, maxSeconds, playing, speed, utClock, togglePlay, seek, setSpeed } =
    usePlayback();

  return (
    <div className="scrubber">
      <button
        className="scrubber-play"
        onClick={togglePlay}
        aria-label={playing ? "Pause" : "Play"}
      >
        {playing ? "⏸" : "▶"}
      </button>

      <input
        type="range"
        className="scrubber-range"
        min={0}
        max={maxSeconds}
        step={60}
        value={currentT}
        onChange={(e) => seek(Number(e.target.value))}
        aria-label="Playback time"
      />

      <select
        className="scrubber-speed"
        value={speed}
        onChange={(e) => setSpeed(Number(e.target.value))}
        aria-label="Playback speed"
      >
        {SPEEDS.map((s) => (
          <option key={s.value} value={s.value}>{s.label}</option>
        ))}
      </select>

      <div className="scrubber-clock">
        <span className="scrubber-utc">{formatUTC(utClock)}</span>
        <span className="scrubber-elapsed">{formatHours(currentT)}</span>
      </div>
    </div>
  );
}
