// Copies out/batch.json and out/sweeps.json into public/data/ so the app can
// `fetch()` them as static assets. Run automatically before `dev` and
// `build` -- see package.json -- and can be run standalone after
// regenerating either file from Python (`sim.export`, `sim.sweeps`).
//
// Deliberately a copy, not a symlink: `vite build` needs real files under
// public/ to include in dist/, and a plain copy keeps that working the same
// way in dev and in a production build.
import { copyFileSync, existsSync, mkdirSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const repoRoot = join(here, "..", "..");
const outDir = join(repoRoot, "out");
const destDir = join(here, "..", "public", "data");

const files = ["batch.json", "sweeps.json"];

mkdirSync(destDir, { recursive: true });

let missing = [];
for (const name of files) {
  const src = join(outDir, name);
  if (!existsSync(src)) {
    missing.push(src);
    continue;
  }
  copyFileSync(src, join(destDir, name));
  console.log(`synced ${name}`);
}

if (missing.length) {
  console.error(
    "Missing source file(s):\n  " +
      missing.join("\n  ") +
      "\nGenerate them first:\n" +
      "  .venv/bin/python -m sim.export   (writes out/batch.json)\n" +
      "  .venv/bin/python -m sim.sweeps   (writes out/sweeps.json)\n"
  );
  process.exit(1);
}
