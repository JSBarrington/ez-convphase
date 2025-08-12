#!/usr/bin/env bash
#
# This script runs a complete bioinformatics pipeline in parallel.
# It takes one argument: the path to the folder containing your FASTA files.
#
# --- Output directory handling ---
# Stop the script if any command fails
set -euo pipefail
IN_DIR="${1:-}"
if [[ -z "${IN_DIR}" ]]; then
  echo "Usage: $0 <input_folder_of_fastas> [output_folder]" >&2
  exit 2
fi

# Default OUT if not provided: a sibling folder named 'ezconvphase_out'
OUT_DIR="${2:-"$(cd "${IN_DIR}/.." && pwd)/ezconvphase_out"}"

# Back-compat: if the rest of the script uses $OUT, set it here:
OUT="$OUT_DIR"

# Standard subfolders
mkdir -p "$OUT_DIR"/{clean,mask,phased,rebuilt,logs}
echo "[Output] $OUT_DIR"

# -------------------------------------------------------------------
# ## --- USER CONFIGURATION --- ##
# -------------------------------------------------------------------
JOBS=4             # Number of parallel jobs to run (0 for all cores)
ITER_VAL=10000     # ConvPhase MCMC Iterations
THIN_VAL=10        # Thinning interval
BURN_VAL=1000      # Burn-in
# -------------------------------------------------------------------


# --- Setup & Path Definitions (Automatic) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Use the first command-line argument ($1) as the input directory.
INDIR="$1"
OUTDIR="$ROOT_DIR/${OUT_DIR}"

SANITIZER="$SCRIPT_DIR/sanitize_for_convphase.py"
REBUILDER="$SCRIPT_DIR/rebuild_from_mask.py"

# --- Sanity Checks ---
[[ -f "$SANITIZER" ]] || { echo "ERROR: Missing sanitizer script: $SANITIZER"; exit 1; }
[[ -f "$REBUILDER" ]] || { echo "ERROR: Missing rebuilder script: $REBUILDER"; exit 1; }
[[ -d "$INDIR" ]] || { echo "ERROR: Input directory not found: $INDIR"; exit 1; }


# --- Main Processing Function ---
process_one_locus() {
  local fa_path="$1"
  local base; base="$(basename "$fa_path")"
  local stem="${base%.*}"
  local log_file="$OUTDIR/logs/${stem}.log"
  local clean_fa="$OUTDIR/clean/${stem}.clean.fasta"
  local mask_json="$OUTDIR/mask/${stem}.mask.json"
  local phased_prefix="$OUTDIR/phased/${stem}"
  local phased_fa_final="$OUTDIR/phased/${stem}.phased.fasta"
  local rebuilt_fa="$OUTDIR/rebuilt/${stem}.rebuilt.fasta"

  if [[ -f "$rebuilt_fa" ]]; then
      echo "Output for $stem already exists. Skipping."
      return
  fi

  {
    echo "[$(date)] START: $stem"
    echo "  Sanitizing..."
    python "$SANITIZER" "$fa_path" "$clean_fa" "$mask_json"
    echo "  Phasing..."
    convphase "$clean_fa" "$phased_prefix" -n "$ITER_VAL" -t "$THIN_VAL" -b "$BURN_VAL"
    local src_phased=""
    if   [[ -s "${phased_prefix}.fas" ]]; then src_phased="${phased_prefix}.fas";
    elif [[ -s "${phased_prefix}.fasta" ]]; then src_phased="${phased_prefix}.fasta";
    elif [[ -s "${phased_prefix}" ]]; then src_phased="${phased_prefix}";
    fi
    if [[ -z "$src_phased" ]]; then
        echo "  ERROR: Phased output file not found for $phased_prefix"
        exit 24
    fi
    mv -f "$src_phased" "$phased_fa_final"
    echo "  Rebuilding..."
    python "$REBUILDER" "$phased_fa_final" "$mask_json" "$rebuilt_fa"
    echo "[$(date)] DONE: $stem"
  } >"$log_file" 2>&1
}

export -f process_one_locus
export OUTDIR SANITIZER REBUILDER ITER_VAL THIN_VAL BURN_VAL

# --- Parallel Execution ---
mkdir -p "$OUTDIR/clean" "$OUTDIR/mask" "$OUTDIR/phased" "$OUTDIR/rebuilt" "$OUTDIR/logs"

echo "--- Starting ConvPhase Pipeline ---"
echo "Input Directory:  $INDIR"
echo "Output Directory: $OUTDIR"
echo "Running up to $JOBS jobs in parallel..."
echo "-----------------------------------"

find "$INDIR" -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fas" \) -print0 | \
  parallel -0 --jobs "$JOBS" --bar 'process_one_locus {}'

echo "-----------------------------------"
echo "All jobs complete. Check '$OUTDIR' for results."
echo "-----------------------------------"