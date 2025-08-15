# ez-convphase

**ez-convphase** is a wrapper for phasing ragged FASTA alignments using itaxotools ConvPhase (implements SeqPhase and PHASE2.1). Cleans headers, masks any problematic columns, runs **ConvPhase** in **parallel**, then rebuilds the original alignment around the phased positions (column-for-column). Supports standard, Haploview & MolD header styles. TSV functions not supported at this time.

> Not affiliated with iTaxoTools/ConvPhase/PHASE.

------------
**Added a GUI!**
------------

---

## Why?

- Alignments can contain tri-allelic, missing and gap characters, which are not supported by ConvPhase's conservative input requirements.
- The scripts in this wrapper remove columns with prohibited characters and replaces them after phasing, restoring the **original alignment length**.
- `ez-convphase` does: *sanitize → phase → rebuild* — and preserves your IDs.

---

## GUI Interface

- EZ-ConvPhase now includes a user-friendly GUI for easier pipeline management.

### Installation

# Install Python dependencies
pip install -r requirements.txt

- Ensure system dependencies
- macOS: brew install parallel
- Ubuntu: sudo apt-get install parallel

---

## Features

- Header-safe: works with `sample|Species` and `sample.Species`.
- Hap suffix mapping on the **left token only**: `_a/_b`, `-a/-b`, `.a/.b`, `/a,/b`, space-`a/b`, trailing `a/b`.
- Exact reconstruction via JSON mask (`kept_indices` + per-column chars).
- Parallel ConvPhase runs (GNU `parallel`).
- Clear logs; per-stage outputs.

---

**GNU parallel**

- macOS: `brew install parallel`
- Debian/Ubuntu: `sudo apt-get update && sudo apt-get install -y parallel`

**ConvPhase** (from iTaxoTools)

- Get the ConvPhase CLI/binary from the iTaxoTools release for your OS.
- Ensure the binary (often `itaxotools-convphase` or `convphase`) is on your `PATH`, **or** set its path in the script via `CONVPHASE_BIN`.
     
**ConvPhase (CLI)**

- Install from PyPI: `python -m pip install itaxotools-convphase`
- Command: convphase (or python -m itaxotools.convphase as a fallback)
- If you keep a local folder/binary (e.g., ConvPhase CLI/convphase), set CONVPHASE_BIN to that path

---
 
**Self-check:**
# does ConvPhase run? (package installs the 'convphase' command)
"${CONVPHASE_BIN:-convphase}" -h >/dev/null 2>&1 \
  || echo "ConvPhase CLI not found. Install with 'python -m pip install itaxotools-convphase' or set CONVPHASE_BIN to your binary."
