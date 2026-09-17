#!/usr/bin/env python3

import csv
import json
import os
import re
import sys
from pathlib import Path

import ijson
import requests


BULK_JSON = Path("addgene_bulk/plasmids_with_sequences_download")
BEST_FULL_TSV = Path("results/itr_positive_best_full_sequences.tsv")
DELTA11_PLASMID_TSV = Path("results/delta11_exact_plasmid_summary.tsv")
GENBANK_DIR = Path("genbank_files")
OUT_TSV = Path("results/test_bulk_vs_genbank_sequence_consistency.tsv")

API_BASE = "https://api.developers.addgene.org"
TOKEN_ENV = "ADDGENE_API_TOKEN"

DELTA11_FORWARD = "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG"
DELTA11_REVERSE = "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA"
MOTIF_LEN = len(DELTA11_FORWARD)


def load_delta11_positive_plasmids(path):
    ids = set()

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

        if "plasmid_id" not in fieldnames:
            raise ValueError(f"{path} missing plasmid_id column")

        count_col = None
        for candidate in [
            "max_exact_delta11_count_per_sequence",
            "max_delta11_exact_total_count",
            "delta11_exact_total_count",
            "max_exact_delta11_count",
            "delta11_exact_count",
        ]:
            if candidate in fieldnames:
                count_col = candidate
                break

        if count_col is None:
            raise ValueError(
                f"{path} missing recognizable delta11 count column. Columns: {fieldnames}"
            )

        for row in reader:
            try:
                count = int(row[count_col])
            except Exception:
                continue

            if count >= 1:
                ids.add(str(row["plasmid_id"]))

    return ids


def load_first_n_delta11_positive_best_full_sequences(best_path, delta11_positive_ids, n=20):
    rows = []

    with best_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

        required = ["plasmid_id", "sequence_id"]
        for col in required:
            if col not in fieldnames:
                raise ValueError(f"{best_path} missing required column: {col}")

        for row in reader:
            pid = str(row["plasmid_id"])
            if pid not in delta11_positive_ids:
                continue

            rows.append(row)

            if len(rows) >= n:
                break

    return rows


def extract_bulk_sequences_for_targets(bulk_json_path, targets):
    """
    Return dict keyed by sequence_id containing bulk sequence information.

    The bulk JSON is streamed with ijson.
    """
    target_sequence_ids = {str(r["sequence_id"]) for r in targets}
    found = {}

    with bulk_json_path.open("rb") as f:
        for plasmid in ijson.items(f, "plasmids.item"):
            plasmid_id = str(plasmid.get("id", ""))

            sequences = plasmid.get("sequences") or []

            if isinstance(sequences, dict):
                sequence_groups = sequences.values()
            else:
                sequence_groups = [sequences]

            for group in sequence_groups:
                if isinstance(group, dict):
                    group = [group]
                if not isinstance(group, list):
                    continue

                for seqrec in group:
                    if not isinstance(seqrec, dict):
                        continue
                    sid = str(seqrec.get("id", ""))

                    if not sid or sid == "None":
                        url = str(seqrec.get("genbank_api_url") or seqrec.get("genbank_url") or "")
                        m = re.search(r"/genbank/(\d+)/?", url)
                        if m:
                            sid = m.group(1)

                    if sid not in target_sequence_ids:
                        continue

                    seq = (
                        seqrec.get("sequence")
                        or seqrec.get("bases")
                        or seqrec.get("nucleotides")
                        or ""
                    )

                    seq = clean_sequence(seq)

                    found[sid] = {
                        "plasmid_id": plasmid_id,
                        "sequence_id": sid,
                        "sequence": seq,
                        "bulk_length": len(seq),
                        "sequence_type": str(seqrec.get("sequence_type", "")),
                        "sequence_source": str(seqrec.get("source", "")),
                        "sequence_name": str(seqrec.get("name", "")),
                        "is_full_sequence": str(seqrec.get("is_full_sequence", "")),
                    }

            if len(found) == len(target_sequence_ids):
                break

    return found


def clean_sequence(seq):
    if seq is None:
        return ""
    seq = str(seq)
    seq = re.sub(r"[^A-Za-z]", "", seq)
    return seq.upper()


def genbank_path_for_sequence_id(sequence_id):
    return GENBANK_DIR / f"sequence_{sequence_id}.gb"


def download_genbank_if_needed(sequence_id):
    GENBANK_DIR.mkdir(exist_ok=True)

    path = genbank_path_for_sequence_id(sequence_id)
    if path.exists() and path.stat().st_size > 0:
        return path, "cached"

    token = os.environ.get(TOKEN_ENV)
    if not token:
        raise RuntimeError(
            f"{TOKEN_ENV} is not set, and GenBank file is missing: {path}"
        )

    url = f"{API_BASE}/download/genbank/{sequence_id}/"
    headers = {"Authorization": f"Token {token}"}

    r = requests.get(url, headers=headers, timeout=120)
    if r.status_code != 200:
        raise RuntimeError(
            f"Failed to download GenBank for sequence_id={sequence_id}: "
            f"HTTP {r.status_code}; response starts: {r.text[:300]}"
        )

    path.write_text(r.text)
    return path, "downloaded"


def parse_genbank_origin_sequence(gb_path):
    """
    Minimal GenBank ORIGIN parser.

    Returns uppercase A/C/G/T/N sequence from ORIGIN to //.
    """
    in_origin = False
    chunks = []

    with gb_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue

            if in_origin:
                if line.startswith("//"):
                    break
                chunks.append(line)

    seq = "".join(chunks)
    seq = clean_sequence(seq)
    return seq


def find_motif_positions_circular(seq, motif):
    """
    Find exact motif positions on circular sequence.

    Returns sorted unique 0-based start positions in original sequence coordinates.
    Uses seq + seq[:motif_len-1], not seq + seq.
    """
    seq = clean_sequence(seq)
    motif = clean_sequence(motif)

    if not seq or not motif:
        return []

    n = len(seq)
    extended = seq + seq[: len(motif) - 1]
    positions = []

    start = 0
    while True:
        idx = extended.find(motif, start)
        if idx == -1:
            break

        if idx < n:
            positions.append(idx)

        start = idx + 1

    return sorted(set(positions))


def delta11_positions(seq):
    fwd = find_motif_positions_circular(seq, DELTA11_FORWARD)
    rev = find_motif_positions_circular(seq, DELTA11_REVERSE)

    all_hits = []
    for p in fwd:
        all_hits.append(("forward", p))
    for p in rev:
        all_hits.append(("reverse", p))

    all_hits.sort(key=lambda x: (x[1], x[0]))
    return all_hits


def format_hits(hits):
    if not hits:
        return ""
    return ";".join([f"{orient}:{pos0 + 1}" for orient, pos0 in hits])


def main():
    if not BULK_JSON.exists():
        raise FileNotFoundError(BULK_JSON)
    if not BEST_FULL_TSV.exists():
        raise FileNotFoundError(BEST_FULL_TSV)
    if not DELTA11_PLASMID_TSV.exists():
        raise FileNotFoundError(DELTA11_PLASMID_TSV)

    OUT_TSV.parent.mkdir(exist_ok=True)
    GENBANK_DIR.mkdir(exist_ok=True)

    delta11_positive_ids = load_delta11_positive_plasmids(DELTA11_PLASMID_TSV)
    print(f"Delta11-positive plasmids loaded: {len(delta11_positive_ids)}")

    targets = load_first_n_delta11_positive_best_full_sequences(
        BEST_FULL_TSV,
        delta11_positive_ids,
        n=20,
    )
    print(f"Testing first {len(targets)} delta11-positive best full sequence records")

    if not targets:
        raise RuntimeError("No target sequence records found")

    print("Streaming bulk JSON to recover target sequences...")
    bulk_by_sid = extract_bulk_sequences_for_targets(BULK_JSON, targets)
    print(f"Bulk sequences found: {len(bulk_by_sid)} / {len(targets)}")

    rows = []

    for i, target in enumerate(targets, start=1):
        plasmid_id = str(target["plasmid_id"])
        sequence_id = str(target["sequence_id"])

        print(f"[{i}/{len(targets)}] sequence_id={sequence_id} plasmid_id={plasmid_id}")

        bulk_info = bulk_by_sid.get(sequence_id)
        bulk_seq = bulk_info["sequence"] if bulk_info else ""

        gb_status = ""
        gb_path = ""
        gb_seq = ""
        gb_error = ""

        try:
            gb_file, gb_status = download_genbank_if_needed(sequence_id)
            gb_path = str(gb_file)
            gb_seq = parse_genbank_origin_sequence(gb_file)
        except Exception as e:
            gb_error = repr(e)

        bulk_hits = delta11_positions(bulk_seq)
        gb_hits = delta11_positions(gb_seq)

        length_match = bool(bulk_seq and gb_seq and len(bulk_seq) == len(gb_seq))
        sequence_match = bool(bulk_seq and gb_seq and bulk_seq == gb_seq)
        delta11_positions_match = bulk_hits == gb_hits

        rows.append(
            {
                "plasmid_id": plasmid_id,
                "sequence_id": sequence_id,
                "plasmid_name": target.get("plasmid_name", ""),
                "bulk_found": "yes" if bulk_info else "no",
                "genbank_status": gb_status,
                "genbank_path": gb_path,
                "genbank_error": gb_error,
                "bulk_length": len(bulk_seq),
                "genbank_length": len(gb_seq),
                "length_match": "yes" if length_match else "no",
                "sequence_match": "yes" if sequence_match else "no",
                "bulk_delta11_count": len(bulk_hits),
                "genbank_delta11_count": len(gb_hits),
                "delta11_positions_match": "yes" if delta11_positions_match else "no",
                "bulk_delta11_hits_1based": format_hits(bulk_hits),
                "genbank_delta11_hits_1based": format_hits(gb_hits),
            }
        )

    fieldnames = [
        "plasmid_id",
        "sequence_id",
        "plasmid_name",
        "bulk_found",
        "genbank_status",
        "genbank_path",
        "genbank_error",
        "bulk_length",
        "genbank_length",
        "length_match",
        "sequence_match",
        "bulk_delta11_count",
        "genbank_delta11_count",
        "delta11_positions_match",
        "bulk_delta11_hits_1based",
        "genbank_delta11_hits_1based",
    ]

    with OUT_TSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote: {OUT_TSV}")

    n = len(rows)
    length_ok = sum(1 for r in rows if r["length_match"] == "yes")
    seq_ok = sum(1 for r in rows if r["sequence_match"] == "yes")
    pos_ok = sum(1 for r in rows if r["delta11_positions_match"] == "yes")

    print(f"Records tested: {n}")
    print(f"Length matches: {length_ok} / {n}")
    print(f"Exact sequence matches: {seq_ok} / {n}")
    print(f"Delta11 position matches: {pos_ok} / {n}")

    if seq_ok == n and pos_ok == n:
        print("PASS: bulk and GenBank coordinates appear directly comparable for this test set.")
    elif length_ok == n and pos_ok == n:
        print("PARTIAL PASS: lengths and delta11 positions match, but full sequences differ.")
    else:
        print("WARNING: coordinate consistency is not fully established.")


if __name__ == "__main__":
    main()
