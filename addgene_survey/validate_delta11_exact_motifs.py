#!/usr/bin/env python3

import csv
from collections import Counter
import ijson

BULK_JSON = "addgene_bulk/plasmids_with_sequences_download"
CANDIDATE_IDS = "results/refined_aav_candidate_ids.txt"
OUT = "results/delta11_exact_motif_validation.tsv"

SEQ_GROUPS = [
    "public_addgene_full_sequences",
    "public_user_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_partial_sequences",
]

MOTIFS = {
    "delta11_forward": "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG",
    "delta11_reverse": "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA",
}


def clean(seq):
    return "".join(c for c in str(seq).upper() if c in "ACGTN")


def load_ids(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    candidate_ids = load_ids(CANDIDATE_IDS)

    rows = []
    plasmids_seen = 0
    candidate_seen = 0
    seqs_scanned = 0

    with open(BULK_JSON, "rb") as f:
        for plasmid in ijson.items(f, "plasmids.item"):
            plasmids_seen += 1

            if plasmids_seen % 10000 == 0:
                print(
                    f"Processed {plasmids_seen:,} plasmids; "
                    f"candidate plasmids seen {candidate_seen:,}; "
                    f"sequences scanned {seqs_scanned:,}",
                    flush=True,
                )

            pid = str(plasmid.get("id", ""))

            if pid not in candidate_ids:
                continue

            candidate_seen += 1
            seqs = plasmid.get("sequences") or {}

            for group in SEQ_GROUPS:
                for seq_entry in seqs.get(group) or []:
                    raw = seq_entry.get("sequence") or ""
                    if not raw:
                        continue

                    seqs_scanned += 1
                    seq = clean(raw)

                    # Use seq + seq so motifs that cross the circular origin are counted.
                    seq2 = seq + seq[:33]

                    forward_count = seq2.count(MOTIFS["delta11_forward"])
                    reverse_count = seq2.count(MOTIFS["delta11_reverse"])
                    total_count = forward_count + reverse_count

                    rows.append({
                        "plasmid_id": pid,
                        "plasmid_name": plasmid.get("name", ""),
                        "plasmid_url": plasmid.get("url", ""),
                        "sequence_group": group,
                        "sequence_id": seq_entry.get("sequence_id", ""),
                        "sequence_description": seq_entry.get("sequence_description", ""),
                        "sequence_length": len(seq),
                        "delta11_forward_exact_count": forward_count,
                        "delta11_reverse_exact_count": reverse_count,
                        "delta11_exact_total_count": total_count,
                        "has_delta11_exact": total_count > 0,
                    })

    fields = [
        "plasmid_id",
        "plasmid_name",
        "plasmid_url",
        "sequence_group",
        "sequence_id",
        "sequence_description",
        "sequence_length",
        "delta11_forward_exact_count",
        "delta11_reverse_exact_count",
        "delta11_exact_total_count",
        "has_delta11_exact",
    ]

    with open(OUT, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    seq_with = sum(1 for r in rows if r["has_delta11_exact"])
    plasmids_with = {r["plasmid_id"] for r in rows if r["has_delta11_exact"]}

    print()
    print("Done.")
    print(f"Candidate plasmids seen: {candidate_seen:,}")
    print(f"Sequences scanned: {seqs_scanned:,}")
    print(f"Sequence records with exact delta11 motif: {seq_with:,}")
    print(f"Plasmids with exact delta11 motif: {len(plasmids_with):,}")
    print(f"Wrote: {OUT}")

    by_total = Counter(r["delta11_exact_total_count"] for r in rows)

    print()
    print("Exact delta11 motif count per sequence record:")
    for k, v in sorted(by_total.items()):
        print(f"{k}\t{v:,}")


if __name__ == "__main__":
    main()
