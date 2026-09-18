import ijson
import csv
from pathlib import Path

INFILE = Path("addgene_bulk/plasmids_with_sequences_download")
OUTFILE = Path("results/metadata_aav_itr_candidates.tsv")

OUTFILE.parent.mkdir(parents=True, exist_ok=True)

terms = [
    "aav",
    "adeno-associated",
    "adeno associated",
    "adenoassociated",
    "itr",
    "inverted terminal repeat",
    "terminal repeat",
    "paav",
    "raav",
]

sequence_groups = [
    "public_addgene_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_full_sequences",
    "public_user_partial_sequences",
]

def flatten_text(obj):
    """
    Recursively collect primitive values from nested dict/list structures.
    """
    parts = []

    if obj is None:
        return ""

    if isinstance(obj, (str, int, float, bool)):
        return str(obj)

    if isinstance(obj, dict):
        for k, v in obj.items():
            parts.append(str(k))
            parts.append(flatten_text(v))
        return " ".join(parts)

    if isinstance(obj, list):
        for x in obj:
            parts.append(flatten_text(x))
        return " ".join(parts)

    return str(obj)

rows = []
plasmid_count = 0
hit_count = 0

with open(INFILE, "rb") as f:
    for plasmid in ijson.items(f, "plasmids.item"):
        plasmid_count += 1

        text = flatten_text(plasmid).lower()
        matched = sorted({term for term in terms if term in text})

        if matched:
            hit_count += 1

            sequences = plasmid.get("sequences") or {}
            seq_counts = {
                group: len(sequences.get(group) or [])
                for group in sequence_groups
            }

            total_sequences = sum(seq_counts.values())

            rows.append({
                "plasmid_id": plasmid.get("id"),
                "name": plasmid.get("name"),
                "description": plasmid.get("description"),
                "url": plasmid.get("url"),
                "matched_terms": ",".join(matched),
                "total_sequences": total_sequences,
                **seq_counts,
            })

        if plasmid_count % 10000 == 0:
            print(f"Processed {plasmid_count:,} plasmids; hits so far: {hit_count:,}")

with OUTFILE.open("w", newline="", encoding="utf-8") as out:
    fieldnames = [
        "plasmid_id",
        "name",
        "description",
        "url",
        "matched_terms",
        "total_sequences",
        *sequence_groups,
    ]
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("\nDone.")
print(f"Processed plasmids: {plasmid_count:,}")
print(f"Metadata candidates: {hit_count:,}")
print(f"Wrote: {OUTFILE}")
