#!/usr/bin/env python3
import sys, json, argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

ALLOWED = set("ACGTRYSWKM")     # ConvPhase-acceptable (ACGT + biallelic codes)
FORBIDDEN = set("N?-BDHV")      # missing + tri-allelic + gaps

def read_rectangular(path: str):
    aln = AlignIO.read(path, "fasta")
    L = aln.get_alignment_length()
    lens = [len(rec.seq) for rec in aln]
    if len(set(lens)) != 1:
        sys.stderr.write(f"[ERROR] Ragged input (lengths={sorted(set(lens))}). Provide a proper aligned file.\n")
        sys.exit(2)
    return aln, L

def find_mask(aln, L):
    n = len(aln)
    kept, removed = [], []
    # record removed chars per column per sequence (by index order)
    removed_cols = []  # list of {"index": i, "chars": [char_for_seq0, ...]}
    for i in range(L):
        col = [aln[j, i].upper() for j in range(n)]
        if all(c in ALLOWED for c in col):
            kept.append(i)
        else:
            removed.append(i)
            removed_cols.append({"index": i, "chars": col})
    return kept, removed, removed_cols

def build_clean(aln, kept):
    clean = []
    for rec in aln:
        seq = "".join(rec.seq[i] for i in kept)
        r = rec[:0]         # copy metadata (id, description)
        r.seq = Seq(str(seq))  # ensure Seq object
        clean.append(r)
    return MultipleSeqAlignment(clean)

def write_fasta_unwrapped(path, aln):
    # Write one-line-per-sequence FASTA (avoids wrap-length confusion)
    with open(path, "w") as f:
        for rec in aln:
            f.write(f">{rec.id}\n{str(rec.seq)}\n")

def main():
    ap = argparse.ArgumentParser(
        description="ConvPhase-safe sanitizer: global-column removal of forbidden chars; stores reinsertion mask."
    )
    ap.add_argument("in_fa")
    ap.add_argument("out_fa")
    ap.add_argument("out_mask_json")
    ap.add_argument("--allow", default="ACGTRYSWKM", help="Allowed symbols (default ACGTRYSWKM)")
    ap.add_argument("--forbid", default="N?-BDHV", help="Forbidden symbols (default N?-BDHV)")
    args = ap.parse_args()

    global ALLOWED, FORBIDDEN
    ALLOWED = set(args.allow.upper())
    FORBIDDEN = set(args.forbid.upper())

    aln, L = read_rectangular(args.in_fa)
    seq_ids = [rec.id for rec in aln]

    kept, removed, removed_cols = find_mask(aln, L)
    if not kept:
        sys.stderr.write("[ERROR] All columns removed; nothing left to analyze.\n")
        sys.exit(3)

    clean = build_clean(aln, kept)
    write_fasta_unwrapped(args.out_fa, clean)

    mask = {
        "original_length": L,
        "seq_ids": seq_ids,               # in order matching rows in removed_cols["chars"]
        "kept_indices": kept,
        "removed_indices": removed,
        "removed_columns": removed_cols,  # list of {index, chars: [per-seq char]}
        "allowed": sorted(list(ALLOWED)),
        "forbidden": sorted(list(FORBIDDEN)),
    }
    with open(args.out_mask_json, "w") as f:
        json.dump(mask, f, indent=2)

    kept_pct = 100.0 * len(kept) / L
    sys.stderr.write(f"[OK] Kept {len(kept)}/{L} sites ({kept_pct:.1f}%). Output length={clean.get_alignment_length()} cols, records={len(clean)}.\n")

if __name__ == "__main__":
    main()
