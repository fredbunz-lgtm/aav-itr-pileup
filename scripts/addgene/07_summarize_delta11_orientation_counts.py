#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict

IN_TSV = "results/delta11_exact_motif_validation.tsv"

seq_pattern_counts = Counter()
plasmid_patterns = defaultdict(list)

with open(IN_TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        pid = row["plasmid_id"]
        sid = row["sequence_id"]

        fwd = int(row["delta11_forward_exact_count"])
        rev = int(row["delta11_reverse_exact_count"])
        total = int(row["delta11_exact_total_count"])

        pattern = (fwd, rev, total)
        seq_pattern_counts[pattern] += 1
        plasmid_patterns[pid].append(pattern)

print("Sequence-record-level forward/reverse/total exact motif count patterns:")
for (fwd, rev, total), n in sorted(seq_pattern_counts.items(), key=lambda x: (-x[1], x[0])):
    print(f"  forward={fwd}\treverse={rev}\ttotal={total}\tsequence_records={n}")

# For each plasmid, choose the sequence record with the highest total exact motif count.
# If tied, choose the one with the highest forward count, then highest reverse count.
plasmid_best_pattern_counts = Counter()

for pid, patterns in plasmid_patterns.items():
    best = sorted(patterns, key=lambda x: (x[2], x[0], x[1]), reverse=True)[0]
    plasmid_best_pattern_counts[best] += 1

print()
print("Plasmid-level best sequence forward/reverse/total exact motif count patterns:")
for (fwd, rev, total), n in sorted(plasmid_best_pattern_counts.items(), key=lambda x: (-x[1], x[0])):
    print(f"  forward={fwd}\treverse={rev}\ttotal={total}\tplasmids={n}")
