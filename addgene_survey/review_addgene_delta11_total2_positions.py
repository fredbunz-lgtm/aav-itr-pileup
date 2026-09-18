#!/usr/bin/env python3

import csv
import ijson
from collections import Counter

BULK_JSON = "addgene_bulk/plasmids_with_sequences_download"
CANDIDATE_IDS = "results/refined_aav_candidate_ids.txt"
OUT_TSV = "results/review_addgene_delta11_total2_positions.tsv"

DELTA11_FORWARD = "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG"
DELTA11_REVERSE = "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA"

SEQ_GROUPS = [
    "public_addgene_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_full_sequences",
    "public_user_partial_sequences",
]

def find_all(seq, motif):
    positions = []
    start = 0
    while True:
        i = seq.find(motif, start)
        if i == -1:
            break
        positions.append(i)
        start = i + 1
    return positions

def distance_bin(d):
    if d == "":
        return "NA"
    d = int(d)
    if d < 100:
        return "<100"
    if d < 200:
        return "100-199"
    if d < 500:
        return "200-499"
    if d < 1000:
        return "500-999"
    if d < 2000:
        return "1000-1999"
    if d < 5000:
        return "2000-4999"
    return "5000+"

with open(CANDIDATE_IDS) as f:
    candidate_ids = {line.strip() for line in f if line.strip()}

out_fields = [
    "plasmid_id",
    "plasmid_name",
    "plasmid_url",
    "sequence_group",
    "sequence_id",
    "sequence_description",
    "sequence_length",
    "forward_count",
    "reverse_count",
    "total_count",
    "forward_positions_0based",
    "reverse_positions_0based",
    "distance_between_two_hits",
    "distance_bin",
]

rows_written = 0
pattern_counter = Counter()
distance_counter = Counter()

with open(OUT_TSV, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_fields)
    writer.writeheader()

    with open(BULK_JSON, "rb") as f:
        for plasmid in ijson.items(f, "plasmids.item"):
            pid = str(plasmid.get("id", ""))
            if pid not in candidate_ids:
                continue

            pname = plasmid.get("name", "")
            purl = plasmid.get("url", "")

            seqs_by_group = plasmid.get("sequences") or {}

            for group in SEQ_GROUPS:
                seqs = seqs_by_group.get(group) or []

                for seqrec in seqs:
                    seq = seqrec.get("sequence") or ""
                    if not seq:
                        continue

                    seq = seq.upper()
                    fwd = find_all(seq, DELTA11_FORWARD)
                    rev = find_all(seq, DELTA11_REVERSE)
                    total = len(fwd) + len(rev)

                    if total != 2:
                        continue

                    all_positions = sorted(fwd + rev)
                    if len(all_positions) == 2:
                        dist = str(all_positions[1] - all_positions[0])
                    else:
                        dist = ""

                    b = distance_bin(dist)

                    pattern = (len(fwd), len(rev), total)
                    pattern_counter[pattern] += 1
                    distance_counter[b] += 1

                    writer.writerow({
                        "plasmid_id": pid,
                        "plasmid_name": pname,
                        "plasmid_url": purl,
                        "sequence_group": group,
                        "sequence_id": seqrec.get("sequence_id", ""),
                        "sequence_description": seqrec.get("sequence_description", ""),
                        "sequence_length": len(seq),
                        "forward_count": len(fwd),
                        "reverse_count": len(rev),
                        "total_count": total,
                        "forward_positions_0based": ",".join(map(str, fwd)),
                        "reverse_positions_0based": ",".join(map(str, rev)),
                        "distance_between_two_hits": dist,
                        "distance_bin": b,
                    })

                    rows_written += 1

print(f"Rows written: {rows_written}")
print(f"Wrote: {OUT_TSV}")

print()
print("Patterns among total=2 sequence records:")
for pattern, n in sorted(pattern_counter.items(), key=lambda x: (-x[1], x[0])):
    fwd, rev, total = pattern
    print(f"  forward={fwd}\treverse={rev}\ttotal={total}\tsequence_records={n}")

print()
print("Distance bins between the two exact motif hits:")
for b, n in sorted(distance_counter.items()):
    print(f"  {b}\t{n}")
