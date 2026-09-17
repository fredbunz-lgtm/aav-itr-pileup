#!/usr/bin/env python3

import csv
from collections import defaultdict, Counter

INFILE = "results/addgene_candidate_itr_hits.tsv"
SEQ_OUT = "results/addgene_candidate_itr_sequence_summary.tsv"
PLASMID_OUT = "results/addgene_candidate_itr_plasmid_summary.tsv"

GROUP_RANK = {
    "public_addgene_full_sequences": 1,
    "public_user_full_sequences": 2,
    "public_addgene_partial_sequences": 3,
    "public_user_partial_sequences": 4,
}

seq_hits = defaultdict(list)

with open(INFILE, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        key = (row["plasmid_id"], row["sequence_group"], row["sequence_id"])
        seq_hits[key].append(row)

seq_rows = []

for key, hits in seq_hits.items():
    first = hits[0]

    classes = Counter(h["itr_class"] for h in hits)
    refs = Counter(h["best_ref"] for h in hits)

    starts = []
    spans = []
    scores = []

    for h in hits:
        try:
            starts.append(int(float(h["start_0based"])))
        except Exception:
            pass
        try:
            spans.append(int(float(h["aligned_plasmid_span"])))
        except Exception:
            pass
        try:
            scores.append(float(h["score"]))
        except Exception:
            pass

    has_delta11_count = sum(
        1 for h in hits
        if str(h.get("has_delta11_junction", "")).lower() in {"true", "1", "yes"}
    )

    seq_rows.append({
        "plasmid_id": first["plasmid_id"],
        "plasmid_name": first["plasmid_name"],
        "plasmid_description": first["plasmid_description"],
        "plasmid_url": first["plasmid_url"],
        "sequence_group": first["sequence_group"],
        "sequence_id": first["sequence_id"],
        "sequence_description": first["sequence_description"],
        "sequence_length_observed": first["sequence_length_observed"],
        "genbank_api_url": first["genbank_api_url"],
        "itr_hit_count": len(hits),
        "has_2_or_more_itrs": len(hits) >= 2,
        "has_delta11_hit_count": has_delta11_count,
        "major_itr_class": classes.most_common(1)[0][0],
        "itr_class_counts": ";".join(f"{k}:{v}" for k, v in classes.most_common()),
        "best_ref_counts": ";".join(f"{k}:{v}" for k, v in refs.most_common()),
        "min_span": min(spans) if spans else "",
        "max_span": max(spans) if spans else "",
        "max_score": max(scores) if scores else "",
        "itr_starts_0based": ";".join(str(x) for x in sorted(starts)),
        "sequence_group_rank": GROUP_RANK.get(first["sequence_group"], 99),
    })

seq_rows.sort(key=lambda r: (
    int(r["plasmid_id"]) if r["plasmid_id"].isdigit() else 999999999,
    int(r["sequence_group_rank"]),
    str(r["sequence_id"]),
))

seq_fields = [
    "plasmid_id",
    "plasmid_name",
    "plasmid_description",
    "plasmid_url",
    "sequence_group",
    "sequence_id",
    "sequence_description",
    "sequence_length_observed",
    "genbank_api_url",
    "itr_hit_count",
    "has_2_or_more_itrs",
    "has_delta11_hit_count",
    "major_itr_class",
    "itr_class_counts",
    "best_ref_counts",
    "min_span",
    "max_span",
    "max_score",
    "itr_starts_0based",
    "sequence_group_rank",
]

with open(SEQ_OUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=seq_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(seq_rows)

by_plasmid = defaultdict(list)
for row in seq_rows:
    by_plasmid[row["plasmid_id"]].append(row)

plasmid_rows = []

for plasmid_id, rows in by_plasmid.items():
    max_any = max(int(r["itr_hit_count"]) for r in rows)
    min_any = min(int(r["itr_hit_count"]) for r in rows)

    full_rows = [
        r for r in rows
        if r["sequence_group"] in {
            "public_addgene_full_sequences",
            "public_user_full_sequences",
        }
    ]

    if full_rows:
        max_full = max(int(r["itr_hit_count"]) for r in full_rows)
        inferred = max_full
        basis = "max_full_sequence_itr_count"
    else:
        max_full = 0
        inferred = max_any
        basis = "max_any_sequence_itr_count"

    best = sorted(
        rows,
        key=lambda r: (
            int(r["sequence_group_rank"]),
            -int(r["itr_hit_count"]),
            -int(float(r["sequence_length_observed"] or 0)),
        )
    )[0]

    group_counts = Counter(r["sequence_group"] for r in rows)

    class_counts = Counter()
    for r in rows:
        for part in r["itr_class_counts"].split(";"):
            if not part:
                continue
            k, v = part.rsplit(":", 1)
            class_counts[k] += int(v)

    plasmid_rows.append({
        "plasmid_id": plasmid_id,
        "plasmid_name": best["plasmid_name"],
        "plasmid_description": best["plasmid_description"],
        "plasmid_url": best["plasmid_url"],
        "sequence_records_with_itr_hits": len(rows),
        "max_itr_hits_in_one_sequence": max_any,
        "min_itr_hits_in_one_sequence": min_any,
        "max_full_sequence_itr_hits": max_full,
        "inferred_plasmid_itr_count": inferred,
        "inference_basis": basis,
        "has_2_or_more_itrs": inferred >= 2,
        "sequence_records_disagree": min_any != max_any,
        "best_sequence_group": best["sequence_group"],
        "best_sequence_id": best["sequence_id"],
        "sequence_group_counts": ";".join(f"{k}:{v}" for k, v in group_counts.most_common()),
        "itr_class_counts_across_hit_sequences": ";".join(f"{k}:{v}" for k, v in class_counts.most_common()),
    })

plasmid_rows.sort(key=lambda r: int(r["plasmid_id"]) if r["plasmid_id"].isdigit() else 999999999)

plasmid_fields = [
    "plasmid_id",
    "plasmid_name",
    "plasmid_description",
    "plasmid_url",
    "sequence_records_with_itr_hits",
    "max_itr_hits_in_one_sequence",
    "min_itr_hits_in_one_sequence",
    "max_full_sequence_itr_hits",
    "inferred_plasmid_itr_count",
    "inference_basis",
    "has_2_or_more_itrs",
    "sequence_records_disagree",
    "best_sequence_group",
    "best_sequence_id",
    "sequence_group_counts",
    "itr_class_counts_across_hit_sequences",
]

with open(PLASMID_OUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=plasmid_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(plasmid_rows)

print(f"Sequence records with ITR hits: {len(seq_rows):,}")
print(f"Plasmids with ITR hits: {len(plasmid_rows):,}")
print(f"Wrote: {SEQ_OUT}")
print(f"Wrote: {PLASMID_OUT}")

counts = Counter(r["inferred_plasmid_itr_count"] for r in plasmid_rows)
print("\nInferred plasmid ITR counts:")
for k, v in sorted(counts.items(), key=lambda x: int(x[0])):
    print(f"{k}\t{v:,}")

disagree = sum(1 for r in plasmid_rows if r["sequence_records_disagree"])
print(f"\nPlasmids with disagreement among hit-containing sequence records: {disagree:,}")
