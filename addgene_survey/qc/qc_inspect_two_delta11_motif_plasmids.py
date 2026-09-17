#!/usr/bin/env python3

import csv
from collections import Counter

import ijson

BULK_JSON = "addgene_bulk/plasmids_with_sequences_download"
TARGETS_TSV = "results/delta11_exact_2motif_plasmids.tsv"
OUT_TSV = "results/delta11_exact_2motif_contexts.tsv"

DELTA11_FORWARD = "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG"
DELTA11_REVERSE = "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA"
MOTIF_LEN = len(DELTA11_FORWARD)

CONTEXT = 100

SEQ_GROUPS = [
    "public_addgene_full_sequences",
    "public_user_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_partial_sequences",
]

def find_all_circular(seq, motif):
    """
    Return motif start positions in circular sequence coordinates.

    Uses seq + seq[:motif_len-1] to detect origin-spanning motifs,
    but only reports starts in the original sequence coordinate range.
    """
    seq2 = seq + seq[:len(motif) - 1]
    positions = []
    start = 0
    n = len(seq)

    while True:
        i = seq2.find(motif, start)
        if i == -1:
            break
        if i < n:
            positions.append(i)
        start = i + 1

    return positions

def circular_context(seq, pos, motif_len, context):
    n = len(seq)
    start = pos - context
    end = pos + motif_len + context

    chars = []
    for i in range(start, end):
        chars.append(seq[i % n])

    left = "".join(chars[:context])
    hit = "".join(chars[context:context + motif_len])
    right = "".join(chars[context + motif_len:])

    return left, hit, right

def circular_distance(a, b, n):
    """
    Shortest circular distance between two motif starts.
    """
    d = abs(a - b)
    return min(d, n - d)

def linear_distance(a, b):
    return abs(a - b)

def classify_orientation(fwd_positions, rev_positions):
    if len(fwd_positions) > 0 and len(rev_positions) == 0:
        return "forward_only"
    if len(rev_positions) > 0 and len(fwd_positions) == 0:
        return "reverse_only"
    if len(fwd_positions) > 0 and len(rev_positions) > 0:
        return "mixed_forward_reverse"
    return "none"

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

target_ids = set()
with open(TARGETS_TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        target_ids.add(str(row["plasmid_id"]))

print(f"Target plasmids: {len(target_ids)}")

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
    "orientation_class",
    "forward_positions_0based",
    "reverse_positions_0based",
    "all_positions_0based",
    "linear_distance_between_two_hits",
    "circular_distance_between_two_hits",
    "linear_distance_bin",
    "circular_distance_bin",
    "hit_label",
    "hit_orientation",
    "hit_index",
    "hit_position_0based",
    "hit_position_1based",
    "context_left",
    "motif_hit",
    "context_right",
]

rows_written = 0
seq_records_with_two = 0
plasmids_seen_with_two = set()
pattern_counter = Counter()
linear_distance_counter = Counter()
circular_distance_counter = Counter()

with open(OUT_TSV, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_fields)
    writer.writeheader()

    with open(BULK_JSON, "rb") as f:
        for plasmid in ijson.items(f, "plasmids.item"):
            pid = str(plasmid.get("id", ""))
            if pid not in target_ids:
                continue

            pname = plasmid.get("name", "")
            purl = plasmid.get("url", "")
            seqs_by_group = plasmid.get("sequences") or {}

            for group in SEQ_GROUPS:
                for seqrec in seqs_by_group.get(group, []) or []:
                    seq = (seqrec.get("sequence") or "").upper()
                    if not seq:
                        continue

                    n = len(seq)
                    fwd_positions = find_all_circular(seq, DELTA11_FORWARD)
                    rev_positions = find_all_circular(seq, DELTA11_REVERSE)
                    total = len(fwd_positions) + len(rev_positions)

                    if total != 2:
                        continue

                    seq_records_with_two += 1
                    plasmids_seen_with_two.add(pid)

                    all_labeled = []
                    for pos in fwd_positions:
                        all_labeled.append(("forward", pos))
                    for pos in rev_positions:
                        all_labeled.append(("reverse", pos))
                    all_labeled.sort(key=lambda x: x[1])

                    all_positions = [pos for orientation, pos in all_labeled]

                    if len(all_positions) == 2:
                        lin_dist = linear_distance(all_positions[0], all_positions[1])
                        circ_dist = circular_distance(all_positions[0], all_positions[1], n)
                        lin_dist_s = str(lin_dist)
                        circ_dist_s = str(circ_dist)
                    else:
                        lin_dist_s = ""
                        circ_dist_s = ""

                    orientation_class = classify_orientation(fwd_positions, rev_positions)
                    pattern = (len(fwd_positions), len(rev_positions), total, orientation_class)
                    pattern_counter[pattern] += 1
                    linear_distance_counter[distance_bin(lin_dist_s)] += 1
                    circular_distance_counter[distance_bin(circ_dist_s)] += 1

                    for hit_i, (orientation, pos) in enumerate(all_labeled, start=1):
                        motif = DELTA11_FORWARD if orientation == "forward" else DELTA11_REVERSE
                        left, hit, right = circular_context(seq, pos, len(motif), CONTEXT)

                        writer.writerow({
                            "plasmid_id": pid,
                            "plasmid_name": pname,
                            "plasmid_url": purl,
                            "sequence_group": group,
                            "sequence_id": seqrec.get("sequence_id", ""),
                            "sequence_description": seqrec.get("sequence_description", ""),
                            "sequence_length": n,
                            "forward_count": len(fwd_positions),
                            "reverse_count": len(rev_positions),
                            "total_count": total,
                            "orientation_class": orientation_class,
                            "forward_positions_0based": ",".join(map(str, fwd_positions)),
                            "reverse_positions_0based": ",".join(map(str, rev_positions)),
                            "all_positions_0based": ",".join(map(str, all_positions)),
                            "linear_distance_between_two_hits": lin_dist_s,
                            "circular_distance_between_two_hits": circ_dist_s,
                            "linear_distance_bin": distance_bin(lin_dist_s),
                            "circular_distance_bin": distance_bin(circ_dist_s),
                            "hit_label": f"hit_{hit_i}",
                            "hit_orientation": orientation,
                            "hit_index": hit_i,
                            "hit_position_0based": pos,
                            "hit_position_1based": pos + 1,
                            "context_left": left,
                            "motif_hit": hit,
                            "context_right": right,
                        })
                        rows_written += 1

print()
print(f"Target plasmids: {len(target_ids)}")
print(f"Target plasmids with at least one total=2 sequence record: {len(plasmids_seen_with_two)}")
print(f"Sequence records with total=2: {seq_records_with_two}")
print(f"Context rows written: {rows_written}")
print(f"Wrote: {OUT_TSV}")

print()
print("Pattern counts:")
for pattern, n in sorted(pattern_counter.items(), key=lambda x: (-x[1], x[0])):
    fwd, rev, total, cls = pattern
    print(f"  forward={fwd}\treverse={rev}\ttotal={total}\tclass={cls}\tsequence_records={n}")

print()
print("Linear distance bins:")
for b, n in sorted(linear_distance_counter.items()):
    print(f"  {b}\t{n}")

print()
print("Circular distance bins:")
for b, n in sorted(circular_distance_counter.items()):
    print(f"  {b}\t{n}")
