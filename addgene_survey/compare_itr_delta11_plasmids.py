#!/usr/bin/env python3

import csv

ITR_SUMMARY = "results/addgene_candidate_itr_plasmid_summary.tsv"
DELTA11_FILE = "results/delta11_exact_motif_validation.tsv"

itr_positive = set()
itr_2plus = set()

with open(ITR_SUMMARY, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        pid = row["plasmid_id"]
        itr_positive.add(pid)

        if int(row["inferred_plasmid_itr_count"]) >= 2:
            itr_2plus.add(pid)

delta11_positive = set()

with open(DELTA11_FILE, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if int(row["delta11_exact_total_count"]) > 0:
            delta11_positive.add(row["plasmid_id"])

both_itr_delta = itr_positive & delta11_positive
both_2plus_delta = itr_2plus & delta11_positive

print(f"ITR-positive plasmids: {len(itr_positive):,}")
print(f"2+ ITR plasmids: {len(itr_2plus):,}")
print(f"Exact delta11-positive plasmids: {len(delta11_positive):,}")
print()

print(
    f"ITR-positive with exact delta11: "
    f"{len(both_itr_delta):,} / {len(itr_positive):,} = "
    f"{len(both_itr_delta) / len(itr_positive):.1%}"
)

print(
    f"2+ ITR with exact delta11: "
    f"{len(both_2plus_delta):,} / {len(itr_2plus):,} = "
    f"{len(both_2plus_delta) / len(itr_2plus):.1%}"
)

print()
print(f"ITR-positive but no exact delta11: {len(itr_positive - delta11_positive):,}")
print(f"Exact delta11 but not ITR-positive: {len(delta11_positive - itr_positive):,}")
