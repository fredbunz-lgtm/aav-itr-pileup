#!/usr/bin/env python3

import csv

import ijson

BULK_JSON = "addgene_bulk/plasmids_with_sequences_download"
ITR_PLASMID_SUMMARY = "results/addgene_candidate_itr_plasmid_summary.tsv"
OUT_TSV = "results/itr_positive_best_full_sequences.tsv"

PREFERRED_GROUPS = [
    "public_addgene_full_sequences",
    "public_user_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_partial_sequences",
]

FULL_GROUPS = {
    "public_addgene_full_sequences",
    "public_user_full_sequences",
}

def load_itr_positive_ids(path):
    ids = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ids.add(str(row["plasmid_id"]))
    return ids

itr_ids = load_itr_positive_ids(ITR_PLASMID_SUMMARY)
print(f"ITR-positive plasmid IDs loaded: {len(itr_ids)}")

rows = []
seen = 0
with_full = 0
with_any_sequence = 0
missing_sequence_id = 0

with open(BULK_JSON, "rb") as f:
    for plasmid in ijson.items(f, "plasmids.item"):
        pid = str(plasmid.get("id", ""))
        if pid not in itr_ids:
            continue

        seen += 1

        seqs_by_group = plasmid.get("sequences") or {}

        chosen_group = ""
        chosen_seq = None

        for group in PREFERRED_GROUPS:
            seqs = seqs_by_group.get(group) or []
            if not seqs:
                continue

            # If multiple records in same preferred group, choose longest sequence.
            chosen_seq = max(seqs, key=lambda s: len(s.get("sequence") or ""))
            chosen_group = group
            break

        if chosen_seq is None:
            rows.append({
                "plasmid_id": pid,
                "plasmid_name": plasmid.get("name", ""),
                "plasmid_url": plasmid.get("url", ""),
                "sequence_group": "",
                "sequence_id": "",
                "sequence_description": "",
                "sequence_length": "",
                "is_full_sequence_group": 0,
            })
            continue

        with_any_sequence += 1
        if chosen_group in FULL_GROUPS:
            with_full += 1

        sequence_id = chosen_seq.get("sequence_id", "")
        if not sequence_id:
            missing_sequence_id += 1

        rows.append({
            "plasmid_id": pid,
            "plasmid_name": plasmid.get("name", ""),
            "plasmid_url": plasmid.get("url", ""),
            "sequence_group": chosen_group,
            "sequence_id": sequence_id,
            "sequence_description": chosen_seq.get("sequence_description", ""),
            "sequence_length": len(chosen_seq.get("sequence") or ""),
            "is_full_sequence_group": int(chosen_group in FULL_GROUPS),
        })

out_fields = [
    "plasmid_id",
    "plasmid_name",
    "plasmid_url",
    "sequence_group",
    "sequence_id",
    "sequence_description",
    "sequence_length",
    "is_full_sequence_group",
]

with open(OUT_TSV, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_fields)
    writer.writeheader()
    writer.writerows(rows)

print()
print(f"ITR-positive plasmids expected: {len(itr_ids)}")
print(f"ITR-positive plasmids seen in bulk: {seen}")
print(f"With any chosen sequence: {with_any_sequence}")
print(f"With chosen full sequence: {with_full}")
print(f"Chosen records missing sequence_id: {missing_sequence_id}")
print(f"Wrote: {OUT_TSV}")
