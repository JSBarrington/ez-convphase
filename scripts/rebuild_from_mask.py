#!/usr/bin/env python3
import sys, json, argparse
from Bio import SeqIO
from Bio.Seq import Seq

# --- map phased IDs -> original IDs (strip hap suffix only from the LEFT token) ---
def guess_parent_id(phased_id: str, known_ids: set[str]) -> str | None:
    if phased_id in known_ids:
        return phased_id

    hap_sfx = ("_a","_b","-a","-b",".a",".b"," a"," b","/a","/b")
    cands = set()

    # Trim whole-ID hap suffixes (common for plain sample IDs)
    for s in hap_sfx:
        if phased_id.endswith(s):
            cands.add(phased_id[:-len(s)])
    if phased_id.endswith(("a","b")):
        cands.add(phased_id[:-1])

    # If there's a delimiter, only strip suffixes from the LEFT token (sample)
    for delim in ("|", "."):
        if delim in phased_id:
            left, right = phased_id.split(delim, 1)
            left_opts = {left}
            for s in hap_sfx:
                if left.endswith(s):
                    left_opts.add(left[:-len(s)])
            if left.endswith(("a","b")):
                left_opts.add(left[:-1])
            for lo in left_opts:
                cands.add(lo + delim + right)

    for c in cands:
        if c in known_ids:
            return c
    return None

def load_mask(mask_path: str):
    m = json.load(open(mask_path))
    # accept either schema
    orig_len = m.get("original_length", m.get("orig_length"))
    original_ids = m.get("seq_ids", m.get("original_ids"))
    kept = m.get("kept_indices")
    removed_cols = m.get("removed_columns")  # [{index:int, chars:[...]}]

    if orig_len is None or original_ids is None or kept is None:
        sys.stderr.write("[ERROR] Mask JSON missing required keys; need original_length/orig_length,"
                         " seq_ids/original_ids, and kept_indices.\n")
        sys.exit(2)
    if removed_cols is None:
        sys.stderr.write("[ERROR] Mask JSON lacks removed_columns (per-column chars). "
                         "Cannot rebuild without it.\n")
        sys.exit(2)

    return orig_len, original_ids, kept, removed_cols

def main():
    ap = argparse.ArgumentParser(description="Rebuild original-length alignment using sanitize mask JSON.")
    ap.add_argument("phased_fa")
    ap.add_argument("mask_json")
    ap.add_argument("out_fa")
    args = ap.parse_args()

    phased = list(SeqIO.parse(args.phased_fa, "fasta"))
    if not phased:
        sys.stderr.write(f"[ERROR] No sequences in {args.phased_fa}\n")
        sys.exit(2)

    orig_len, original_ids, kept, removed_cols = load_mask(args.mask_json)

    kept = list(kept)
    kept_set = set(kept)
    removed_map = {rc["index"]: rc["chars"] for rc in removed_cols}

    if len(kept) + len(removed_map) != orig_len:
        sys.stderr.write("[WARN] kept + removed != original length; continuing.\n")

    id_to_row = {id_: i for i, id_ in enumerate(original_ids)}
    known_ids = set(original_ids)

    rebuilt = []
    for rec in phased:
        parent_id = guess_parent_id(rec.id, known_ids)
        if parent_id is None:
            preview = ", ".join(list(sorted(known_ids))[:5])
            sys.stderr.write(
                f"[ERROR] Phased sequence ID '{rec.id}' could not be mapped to any original ID.\n"
                f"       First few known IDs: {preview} ...\n"
            )
            sys.exit(2)

        row_idx = id_to_row[parent_id]
        pseq = str(rec.seq)
        if len(pseq) != len(kept):
            sys.stderr.write(f"[ERROR] Length mismatch for {rec.id}: phased={len(pseq)} vs kept={len(kept)}.\n")
            sys.exit(2)

        out = []
        kptr = 0
        for i in range(orig_len):
            if i in kept_set:
                out.append(pseq[kptr]); kptr += 1
            else:
                out.append(removed_map[i][row_idx])

        new_rec = rec[:0]
        new_rec.id = rec.id
        new_rec.description = ""
        new_rec.seq = Seq("".join(out))
        rebuilt.append(new_rec)

    with open(args.out_fa, "w") as f:
        for r in rebuilt:
            f.write(f">{r.id}\n{str(r.seq)}\n")

    sys.stderr.write(f"[OK] Rebuilt alignment length={orig_len} for {len(rebuilt)} sequences → {args.out_fa}\n")

if __name__ == "__main__":
    main()
