#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import re
import argparse

ITR_REFS = {
    "145_ITR_Flip": "TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    "145_ITR_Flop": "TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT",
    "145RC_ITR_Flip": "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
    "145RC_ITR_Flop": "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAA",
}

DELTA11_JUNCTIONS = {
    "delta11_forward": "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG",
    "delta11_reverse": "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA",
}


def clean(seq):
    return re.sub("[^ACGTNacgtn]", "", str(seq)).upper()


def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1
    return aligner


def aln_stats(aln):
    q_blocks = aln.aligned[0]
    t_blocks = aln.aligned[1]

    q_aligned = sum(int(e - s) for s, e in q_blocks)
    t_start = int(t_blocks[0][0])
    t_end = int(t_blocks[-1][1])
    t_span = t_end - t_start

    return q_aligned, t_start, t_end, t_span


def circ_dist(a, b, n):
    d = abs(a - b)
    return min(d, n - d)


def classify_length(span, has_delta11):
    if has_delta11:
        return "119_delta11"
    elif 140 <= span <= 150:
        return "145_like"
    elif 126 <= span <= 134:
        return "130_like"
    elif 115 <= span <= 123:
        return "119_like_no_delta11_motif"
    else:
        return f"other_{span}bp"


def find_delta11_near(seq, start, end, margin=80):
    n = len(seq)

    if end >= start:
        itr_len = end - start
    else:
        itr_len = n - start + end

    window_start = start - margin
    window_len = itr_len + 2 * margin
    window = "".join(seq[(window_start + i) % n] for i in range(window_len))

    hits = []
    for name, motif in DELTA11_JUNCTIONS.items():
        p = window.find(motif)
        if p >= 0:
            hits.append((name, (window_start + p) % n))

    return hits


def collapse_hits(raw_hits, n, min_separation=80):
    raw_hits = sorted(raw_hits, key=lambda x: x["score"], reverse=True)
    kept = []

    for h in raw_hits:
        if all(circ_dist(h["start_0based"], k["start_0based"], n) >= min_separation for k in kept):
            kept.append(h)

    return sorted(kept, key=lambda x: x["start_0based"])


def analyze_fasta(path, min_score=180):
    aligner = make_aligner()
    rows = []

    for rec in SeqIO.parse(path, "fasta"):
        seq = clean(str(rec.seq))
        n = len(seq)
        seq2 = seq + seq

        raw_hits = []

        for ref_name, ref in ITR_REFS.items():
            ref = clean(ref)

            aln = aligner.align(ref, seq2)[0]
            q_aligned, t_start, t_end, t_span = aln_stats(aln)

            if aln.score < min_score:
                continue

            # Avoid duplicate alignments starting in second copy of circularized plasmid
            if t_start >= n:
                continue

            raw_hits.append({
                "record_id": rec.id,
                "plasmid_length": n,
                "start_0based": t_start % n,
                "end_0based_exclusive": t_end % n,
                "best_ref": ref_name,
                "score": float(aln.score),
                "aligned_ref_bases": q_aligned,
                "aligned_plasmid_span": t_span,
            })

        hits = collapse_hits(raw_hits, n)

        for idx, h in enumerate(hits, start=1):
            delta_hits = find_delta11_near(
                seq,
                h["start_0based"],
                h["end_0based_exclusive"]
            )

            has_delta = len(delta_hits) > 0

            h["itr_index"] = idx
            h["has_delta11_junction"] = has_delta
            h["delta11_junction"] = ";".join(x[0] for x in delta_hits) if delta_hits else ""
            h["delta11_position_0based"] = ";".join(str(x[1]) for x in delta_hits) if delta_hits else ""
            h["itr_class"] = classify_length(h["aligned_plasmid_span"], has_delta)

            rows.append(h)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("--out", default="itr_test_output.csv")
    parser.add_argument("--min-score", type=float, default=180)
    args = parser.parse_args()

    df = analyze_fasta(args.fasta, min_score=args.min_score)

    if df.empty:
        print("No ITR hits detected.")
    else:
        print(df.to_string(index=False))

    df.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
