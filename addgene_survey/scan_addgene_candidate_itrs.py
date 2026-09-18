#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

import ijson

from test_itr_detection import (
    ITR_REFS,
    clean,
    make_aligner,
    aln_stats,
    classify_length,
    find_delta11_near,
    collapse_hits,
)


SEQ_GROUPS = [
    "public_addgene_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_full_sequences",
    "public_user_partial_sequences",
]


def detect_itrs_in_sequence(seq, record_id, min_score=180):
    """
    Adapted from test_itr_detection.py analyze_fasta(), but works directly
    on one raw sequence string instead of a FASTA file.
    """
    seq = clean(seq)
    n = len(seq)

    if n == 0:
        return []

    seq2 = seq + seq
    aligner = make_aligner()
    raw_hits = []

    for ref_name, ref_seq in ITR_REFS.items():
        for aln in aligner.align(ref_seq, seq2):
            if aln.score < min_score:
                continue

            stats = aln_stats(aln)

            # test_itr_detection.py versions may return either:
            #   dict with target_start/target_end/query_span
            # or:
            #   tuple/list containing query_start, query_end, target_start, target_end
            if isinstance(stats, dict):
                t_start = stats["target_start"]
                t_end = stats["target_end"]
                query_span = stats["query_span"]
            else:
                # test_itr_detection.py returns:
                #   q_aligned, t_start, t_end, t_span
                query_span, t_start, t_end, t_span = stats

            # Skip alignments that start in duplicated second copy.
            if t_start >= n:
                continue

            start = t_start % n
            end = t_end % n
            span = t_end - t_start

            raw_hits.append(
                {
                    "record_id": record_id,
                    "plasmid_length": n,
                    "start_0based": start,
                    "end_0based_exclusive": end,
                    "best_ref": ref_name,
                    "score": aln.score,
                    "aligned_ref_bases": query_span,
                    "aligned_plasmid_span": span,
                }
            )

    hits = collapse_hits(raw_hits, n, min_separation=80)

    for i, hit in enumerate(hits, start=1):
        start = int(hit["start_0based"])
        end = int(hit["end_0based_exclusive"])
        span = int(hit["aligned_plasmid_span"])

        delta = find_delta11_near(seq, start, end, margin=80)

        # test_itr_detection.py versions may return:
        #   None/empty
        #   a dict like {"name": ..., "position": ...}
        #   a list of dicts/tuples/strings
        has_delta = bool(delta)
        delta_name = ""
        delta_pos = ""

        if delta:
            first_delta = delta[0] if isinstance(delta, list) else delta

            if isinstance(first_delta, dict):
                delta_name = first_delta.get("name", "")
                delta_pos = first_delta.get("position", "")
            elif isinstance(first_delta, tuple):
                # Best-effort handling for tuple/list-style delta11 results.
                delta_name = str(first_delta[0]) if len(first_delta) > 0 else ""
                delta_pos = first_delta[1] if len(first_delta) > 1 else ""
            else:
                delta_name = str(first_delta)

        hit["itr_index"] = i
        hit["has_delta11_junction"] = has_delta
        hit["delta11_junction"] = delta_name
        hit["delta11_position_0based"] = delta_pos
        hit["itr_class"] = classify_length(span, has_delta)

    return hits


def load_candidate_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                ids.add(str(s))
    return ids


def iter_sequence_entries(plasmid):
    seqs = plasmid.get("sequences") or {}

    for group in SEQ_GROUPS:
        entries = seqs.get(group) or []
        for entry in entries:
            yield group, entry


def main():
    parser = argparse.ArgumentParser(
        description="Stream Addgene plasmids_with_sequences JSON and scan candidate plasmids for AAV ITRs."
    )
    parser.add_argument(
        "--json",
        default="addgene_bulk/plasmids_with_sequences_download",
        help="Addgene bulk plasmids_with_sequences JSON file",
    )
    parser.add_argument(
        "--candidate-ids",
        default="results/refined_aav_candidate_ids.txt",
        help="One Addgene plasmid ID per line",
    )
    parser.add_argument(
        "--out",
        default="results/addgene_candidate_itr_hits.tsv",
        help="Output TSV of ITR hits",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=180,
        help="Minimum local alignment score for ITR detection",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Progress message interval by plasmids processed from JSON",
    )

    args = parser.parse_args()

    candidate_ids = load_candidate_ids(args.candidate_ids)
    print(f"Loaded candidate plasmid IDs: {len(candidate_ids):,}", file=sys.stderr)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "plasmid_id",
        "plasmid_name",
        "plasmid_description",
        "plasmid_url",
        "sequence_group",
        "sequence_id",
        "sequence_description",
        "sequence_length_reported",
        "sequence_length_observed",
        "genbank_api_url",
        "genbank_url",
        "record_id",
        "plasmid_length",
        "start_0based",
        "end_0based_exclusive",
        "best_ref",
        "score",
        "aligned_ref_bases",
        "aligned_plasmid_span",
        "itr_index",
        "has_delta11_junction",
        "delta11_junction",
        "delta11_position_0based",
        "itr_class",
    ]

    plasmids_seen = 0
    candidate_plasmids_seen = 0
    sequences_scanned = 0
    plasmids_with_hits = set()
    total_hits = 0

    with open(args.json, "rb") as f, open(out_path, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for plasmid in ijson.items(f, "plasmids.item"):
            plasmids_seen += 1

            if plasmids_seen % args.progress_every == 0:
                print(
                    f"Processed {plasmids_seen:,} plasmids from JSON; "
                    f"candidate plasmids seen {candidate_plasmids_seen:,}; "
                    f"sequences scanned {sequences_scanned:,}; "
                    f"ITR hits {total_hits:,}",
                    file=sys.stderr,
                    flush=True,
                )

            plasmid_id = str(plasmid.get("id", ""))

            if plasmid_id not in candidate_ids:
                continue

            candidate_plasmids_seen += 1

            plasmid_name = plasmid.get("name", "")
            plasmid_description = plasmid.get("description", "")
            plasmid_url = plasmid.get("url", "")

            for seq_group, seq_entry in iter_sequence_entries(plasmid):
                raw_seq = seq_entry.get("sequence") or ""
                if not raw_seq:
                    continue

                sequences_scanned += 1

                sequence_id = seq_entry.get("sequence_id", "")
                sequence_description = seq_entry.get("sequence_description", "")
                sequence_length_reported = seq_entry.get("length", "")
                genbank_api_url = seq_entry.get("genbank_api_url", "")
                genbank_url = seq_entry.get("genbank_url", "")

                record_id = f"{plasmid_id}|{seq_group}|{sequence_id}"

                hits = detect_itrs_in_sequence(
                    raw_seq,
                    record_id=record_id,
                    min_score=args.min_score,
                )

                if hits:
                    plasmids_with_hits.add(plasmid_id)

                for hit in hits:
                    total_hits += 1

                    row = {
                        "plasmid_id": plasmid_id,
                        "plasmid_name": plasmid_name,
                        "plasmid_description": plasmid_description,
                        "plasmid_url": plasmid_url,
                        "sequence_group": seq_group,
                        "sequence_id": sequence_id,
                        "sequence_description": sequence_description,
                        "sequence_length_reported": sequence_length_reported,
                        "sequence_length_observed": len(clean(raw_seq)),
                        "genbank_api_url": genbank_api_url,
                        "genbank_url": genbank_url,
                    }
                    row.update(hit)
                    writer.writerow(row)

    print("Done.", file=sys.stderr)
    print(f"Plasmids processed from JSON: {plasmids_seen:,}", file=sys.stderr)
    print(f"Candidate plasmids seen: {candidate_plasmids_seen:,}", file=sys.stderr)
    print(f"Sequences scanned: {sequences_scanned:,}", file=sys.stderr)
    print(f"Plasmids with >=1 ITR hit: {len(plasmids_with_hits):,}", file=sys.stderr)
    print(f"Total ITR hits: {total_hits:,}", file=sys.stderr)
    print(f"Wrote: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
