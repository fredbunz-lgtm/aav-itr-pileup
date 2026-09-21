#!/usr/bin/env python3
"""
scan_ori_motifs_batch.py

RECONSTRUCTION NOTICE
----------------------
This is a reconstructed batch wrapper for the original ori-motif pre-screen
step used in the AAV ITR Addgene survey pipeline (step 13). The original
batch script that looped this scan over thousands of plasmids was not
preserved/committed and could not be recovered. Only the underlying
single-FASTA-file tool, `scan_ori_motifs.py`, survived.

This script reuses the exact motif dictionary (`MOTIFS`) and the exact
circular motif-matching function (`find_all_circular`) from
`scan_ori_motifs.py` unchanged, and adds new batch I/O code (written during
this repository reorganization) to:

  1. Read the list of ITR-positive plasmids and their "best" chosen sequence
     (group + sequence_id) from `results/itr_positive_best_full_sequences.tsv`.
  2. Stream the Addgene bulk JSON download and, for each target plasmid,
     retrieve the specific sequence record that was already selected as
     "best" by `make_itr_positive_best_full_sequences.py`.
  3. Run the ori-motif scan (both strands, circular wraparound) on each
     retrieved sequence.
  4. Write one row per motif hit to `results/ori_motif_scan_results.tsv`.

This script was written after the fact, during documentation/reorganization
work, and was NOT used to generate the results reported in the bioRxiv
submission. It is provided so the pipeline is runnable end-to-end going
forward.
"""

import csv

import ijson

from scan_ori_motifs import MOTIFS, find_all_circular
from Bio.Seq import Seq

BULK_JSON = "addgene_bulk/plasmids_with_sequences_download"
BEST_SEQ_TSV = "results/itr_positive_best_full_sequences.tsv"
OUT_TSV = "results/ori_motif_scan_results.tsv"


def load_targets(path):
    """
    Read results/itr_positive_best_full_sequences.tsv and build a lookup of
    plasmid_id -> (sequence_group, sequence_id, plasmid_name), skipping any
    rows that have no chosen sequence at all.
    """
    targets = {}
    skipped = 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            group = row.get("sequence_group", "")
            seq_id = row.get("sequence_id", "")
            if not group or not seq_id:
                skipped += 1
                continue
            targets[str(row["plasmid_id"])] = (
                group,
                str(seq_id),
                row.get("plasmid_name", ""),
            )
    print(f"Targets loaded: {len(targets)} (skipped {skipped} with no chosen sequence)")
    return targets


def find_sequence(plasmid, group, sequence_id):
    """
    Given a plasmid record from the bulk JSON and a specific (group,
    sequence_id) pair, return the sequence string, or None if not found.
    """
    seqs_by_group = plasmid.get("sequences") or {}
    seqs = seqs_by_group.get(group) or []
    for s in seqs:
        if str(s.get("sequence_id", "")) == sequence_id:
            return s.get("sequence")
    return None


def scan_sequence(seq):
    """
    Run the ori-motif scan on a single sequence string, both strands,
    circular wraparound. Returns a list of hit dicts.
    """
    seq = seq.upper()
    n = len(seq)
    hits = []
    for label, motif in MOTIFS.items():
        rc = str(Seq(motif).reverse_complement())
        for strand, query in [("+", motif), ("-", rc)]:
            for pos in find_all_circular(seq, query):
                end = (pos + len(query)) % n
                hits.append({
                    "motif_label": label,
                    "strand": strand,
                    "start_0based": pos,
                    "end_0based_exclusive": end,
                    "length": len(query),
                })
    return hits


def main():
    targets = load_targets(BEST_SEQ_TSV)

    out_rows = []
    scanned = 0
    missing_sequence = 0
    plasmids_with_hits = 0
    total_hits = 0

    with open(BULK_JSON, "rb") as f:
        for plasmid in ijson.items(f, "plasmids.item"):
            pid = str(plasmid.get("id", ""))
            if pid not in targets:
                continue

            group, sequence_id, plasmid_name = targets[pid]
            sequence = find_sequence(plasmid, group, sequence_id)

            if not sequence:
                missing_sequence += 1
                continue

            scanned += 1
            hits = scan_sequence(sequence)

            if hits:
                plasmids_with_hits += 1
                total_hits += len(hits)
                for h in hits:
                    out_rows.append({
                        "plasmid_id": pid,
                        "plasmid_name": plasmid_name,
                        "motif_label": h["motif_label"],
                        "strand": h["strand"],
                        "start_0based": h["start_0based"],
                        "end_0based_exclusive": h["end_0based_exclusive"],
                        "length": h["length"],
                    })

    with open(OUT_TSV, "w", newline="") as f:
        fieldnames = [
            "plasmid_id",
            "plasmid_name",
            "motif_label",
            "strand",
            "start_0based",
            "end_0based_exclusive",
            "length",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Plasmids scanned: {scanned}")
    print(f"Plasmids missing sequence: {missing_sequence}")
    print(f"Plasmids with at least one hit: {plasmids_with_hits}")
    print(f"Total hits written: {total_hits}")
    print(f"Output written to: {OUT_TSV}")


if __name__ == "__main__":
    main()
