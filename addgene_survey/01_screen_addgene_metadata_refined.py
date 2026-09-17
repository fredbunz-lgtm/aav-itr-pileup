import ijson
import csv
import re
from pathlib import Path

INFILE = Path("addgene_bulk/plasmids_with_sequences_download")
OUTFILE = Path("results/metadata_aav_vector_candidates_refined.tsv")

OUTFILE.parent.mkdir(parents=True, exist_ok=True)

sequence_groups = [
    "public_addgene_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_full_sequences",
    "public_user_partial_sequences",
]

positive_patterns = {
    "pAAV": re.compile(r"\bpAAV\b|\bpAAV[-_.]", re.IGNORECASE),
    "rAAV": re.compile(r"\brAAV\b|\brAAV[-_.]", re.IGNORECASE),
    "AAV_word": re.compile(r"\bAAV\b", re.IGNORECASE),
    "adeno_associated": re.compile(r"adeno[- ]associated", re.IGNORECASE),
    "inverted_terminal_repeat": re.compile(r"inverted terminal repeat", re.IGNORECASE),
    "AAV_ITR": re.compile(r"\bAAV.{0,80}\bITR\b|\bITR.{0,80}\bAAV\b", re.IGNORECASE),
}

exclude_patterns = {
    "AAVS1": re.compile(r"\bAAVS1\b", re.IGNORECASE),
    "in_vitro": re.compile(r"\bin vitro\b", re.IGNORECASE),
}

def flatten_metadata(obj):
    """
    Recursively collect metadata text, but skip raw DNA sequence strings.
    """
    parts = []

    if obj is None:
        return ""

    if isinstance(obj, (str, int, float, bool)):
        return str(obj)

    if isinstance(obj, dict):
        for k, v in obj.items():
            # Critical speed fix: do not search raw sequence bases.
            if k == "sequence":
                continue
            parts.append(str(k))
            parts.append(flatten_metadata(v))
        return " ".join(parts)

    if isinstance(obj, list):
        for x in obj:
            parts.append(flatten_metadata(x))
        return " ".join(parts)

    return str(obj)

rows = []
plasmid_count = 0
hit_count = 0

with open(INFILE, "rb") as f:
    for plasmid in ijson.items(f, "plasmids.item"):
        plasmid_count += 1

        text = flatten_metadata(plasmid)

        positive_hits = [
            name for name, pattern in positive_patterns.items()
            if pattern.search(text)
        ]

        if not positive_hits:
            if plasmid_count % 10000 == 0:
                print(f"Processed {plasmid_count:,} plasmids; refined hits so far: {hit_count:,}")
            continue

        exclude_hits = [
            name for name, pattern in exclude_patterns.items()
            if pattern.search(text)
        ]

        # Exclude likely AAVS1-only records.
        # These are genome-targeting plasmids, not necessarily AAV vectors.
        if "AAVS1" in exclude_hits and not any(
            h in positive_hits for h in ["pAAV", "rAAV", "adeno_associated", "inverted_terminal_repeat", "AAV_ITR"]
        ):
            if plasmid_count % 10000 == 0:
                print(f"Processed {plasmid_count:,} plasmids; refined hits so far: {hit_count:,}")
            continue

        sequences = plasmid.get("sequences") or {}
        seq_counts = {
            group: len(sequences.get(group) or [])
            for group in sequence_groups
        }

        total_sequences = sum(seq_counts.values())

        hit_count += 1

        rows.append({
            "plasmid_id": plasmid.get("id"),
            "name": plasmid.get("name"),
            "description": plasmid.get("description"),
            "url": plasmid.get("url"),
            "positive_hits": ",".join(positive_hits),
            "exclude_hits": ",".join(exclude_hits),
            "total_sequences": total_sequences,
            **seq_counts,
        })

        if plasmid_count % 10000 == 0:
            print(f"Processed {plasmid_count:,} plasmids; refined hits so far: {hit_count:,}")

with OUTFILE.open("w", newline="", encoding="utf-8") as out:
    fieldnames = [
        "plasmid_id",
        "name",
        "description",
        "url",
        "positive_hits",
        "exclude_hits",
        "total_sequences",
        *sequence_groups,
    ]
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("\nDone.")
print(f"Processed plasmids: {plasmid_count:,}")
print(f"Refined metadata candidates: {hit_count:,}")
print(f"Wrote: {OUTFILE}")
