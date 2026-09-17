import ijson
from collections import Counter

path = "addgene_bulk/plasmids_with_sequences_download"

sequence_groups = [
    "public_addgene_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_full_sequences",
    "public_user_partial_sequences",
]

plasmid_count = 0
sequence_count = 0
group_counts = Counter()
lengths = []

with open(path, "rb") as f:
    for plasmid in ijson.items(f, "plasmids.item"):
        plasmid_count += 1

        sequences = plasmid.get("sequences") or {}

        for group in sequence_groups:
            seq_list = sequences.get(group) or []
            group_counts[group] += len(seq_list)
            sequence_count += len(seq_list)

            for seq_entry in seq_list:
                length = seq_entry.get("length")
                if length is not None:
                    lengths.append(int(length))

        if plasmid_count % 10000 == 0:
            print(f"Processed {plasmid_count:,} plasmids...")

print("\nDone.")
print(f"Plasmids: {plasmid_count:,}")
print(f"Sequences: {sequence_count:,}")

print("\nSequences by group:")
for group, count in group_counts.items():
    print(f"  {group}: {count:,}")

if lengths:
    print("\nSequence lengths:")
    print(f"  min: {min(lengths):,}")
    print(f"  max: {max(lengths):,}")
    print(f"  total bp: {sum(lengths):,}")
