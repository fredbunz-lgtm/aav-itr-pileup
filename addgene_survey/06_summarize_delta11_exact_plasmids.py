#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict

IN_TSV = "results/delta11_exact_motif_validation.tsv"
OUT_TSV = "results/delta11_exact_plasmid_summary.tsv"

def pct(n, d):
    if d == 0:
        return "NA"
    return f"{100.0 * n / d:.1f}%"

plasmid_to_counts = defaultdict(list)
plasmid_to_names = {}

with open(IN_TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    fieldnames = reader.fieldnames or []
    print("Input columns:")
    for name in fieldnames:
        print(f"  {name}")

    required = {"plasmid_id", "delta11_exact_total_count"}
    missing = required - set(fieldnames)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    name_col = None
    for candidate in ["plasmid_name", "name"]:
        if candidate in fieldnames:
            name_col = candidate
            break

    for row in reader:
        pid = row["plasmid_id"]
        try:
            count = int(row["delta11_exact_total_count"])
        except ValueError:
            count = 0

        plasmid_to_counts[pid].append(count)

        if name_col:
            plasmid_to_names[pid] = row.get(name_col, "")

summary_rows = []

for pid, counts in plasmid_to_counts.items():
    max_count = max(counts)
    total_across_sequence_records = sum(counts)
    sequence_records = len(counts)
    sequence_records_with_motif = sum(1 for x in counts if x > 0)

    summary_rows.append({
        "plasmid_id": pid,
        "plasmid_name": plasmid_to_names.get(pid, ""),
        "sequence_records_examined": sequence_records,
        "sequence_records_with_exact_delta11": sequence_records_with_motif,
        "max_exact_delta11_count_per_sequence": max_count,
        "sum_exact_delta11_counts_across_sequence_records": total_across_sequence_records,
        "has_exact_delta11": int(max_count >= 1),
        "has_2plus_exact_delta11": int(max_count >= 2),
        "has_4plus_exact_delta11": int(max_count >= 4),
    })

summary_rows.sort(
    key=lambda r: (
        -int(r["max_exact_delta11_count_per_sequence"]),
        r["plasmid_id"],
    )
)

out_fields = [
    "plasmid_id",
    "plasmid_name",
    "sequence_records_examined",
    "sequence_records_with_exact_delta11",
    "max_exact_delta11_count_per_sequence",
    "sum_exact_delta11_counts_across_sequence_records",
    "has_exact_delta11",
    "has_2plus_exact_delta11",
    "has_4plus_exact_delta11",
]

with open(OUT_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=out_fields)
    writer.writeheader()
    writer.writerows(summary_rows)

total_plasmids = len(summary_rows)
max_count_counter = Counter(
    int(r["max_exact_delta11_count_per_sequence"]) for r in summary_rows
)

n_ge1 = sum(1 for r in summary_rows if int(r["max_exact_delta11_count_per_sequence"]) >= 1)
n_ge2 = sum(1 for r in summary_rows if int(r["max_exact_delta11_count_per_sequence"]) >= 2)
n_ge4 = sum(1 for r in summary_rows if int(r["max_exact_delta11_count_per_sequence"]) >= 4)

print()
print(f"Plasmids summarized: {total_plasmids}")
print()
print("Max exact delta11 motif count per plasmid:")
for count, n in sorted(max_count_counter.items()):
    print(f"  {count}\t{n}")

print()
print(f"Plasmids with >=1 exact delta11 motif: {n_ge1} / {total_plasmids} = {pct(n_ge1, total_plasmids)}")
print(f"Plasmids with >=2 exact delta11 motifs: {n_ge2} / {total_plasmids} = {pct(n_ge2, total_plasmids)}")
print(f"Plasmids with >=4 exact delta11 motifs: {n_ge4} / {total_plasmids} = {pct(n_ge4, total_plasmids)}")

print()
print(f"Wrote: {OUT_TSV}")
