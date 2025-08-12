# ez-convphase

ConvPhase wrapper for ragged FASTA alignments. Cleans headers, masks problematic columns, runs **ConvPhase** in parallel, then rebuilds the original alignment around the phased positions (column-for-column). Supports regular, Haploview & MolD header styles. TSV functions not supported.

> Not affiliated with iTaxoTools/ConvPhase/PHASE.

---

## Why?

- Alignments can contain tri-allelic, missing and gap characters, which are not supported by ConvPhase.
- ConvPhase requires strict inputs, this wrappers removes columns with prohibited characters and round-trips back to the **original alignment length**.
- `ez-convphase` does: *sanitize → phase → rebuild* — and preserves your IDs.

---

## Features

- Header-safe: works with `sample|Species` and `sample.Species`.
- Hap suffix mapping on the **left token only**: `_a/_b`, `-a/-b`, `.a/.b`, `/a,/b`, space-`a/b`, trailing `a/b`.
- Exact reconstruction via JSON mask (`kept_indices` + per-column chars).
- Parallel ConvPhase runs (GNU `parallel`).
- Clear logs; per-stage outputs.

---

## Requirements

- macOS or Linux
- **Python 3.9+**
  - `pip install biopython`
  - (GUI) `tkinter` — e.g., `conda install -c conda-forge tk`
- **GNU parallel**
  - macOS: `brew install parallel`
  - Debian/Ubuntu: `sudo apt-get update && sudo apt-get install -y parallel`
- **ConvPhase** (from iTaxoTools)
  - Get the ConvPhase CLI/binary from the iTaxoTools release for your OS.
  - Ensure the binary (often `itaxotools-convphase` or `convphase`) is on your `PATH`, **or** set its path in the script via `CONVPHASE_BIN`.

**Self-check:**
```bash
# does ConvPhase run?
itaxotools-convphase -h 2>/dev/null || echo "ConvPhase not on PATH; set CONVPHASE_BIN in the script"
