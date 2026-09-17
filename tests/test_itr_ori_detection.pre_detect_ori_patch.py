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

ORI_REFS = {
    "pUC_pMB1_core": "TTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTA",
    "RNAII_short": "CTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAA",
    "colE1_RNAII": "ATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAG",
    "ori_adjacent_pBR322": "TCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTC",
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


def circ_dist_points(a, b, n):
    d = abs(a - b)
    return min(d, n - d)


def circular_interval_distance(a_start, a_end, b_start, b_end, n):
    """
    Approximate minimum circular distance between two intervals.
    Assumes intervals do not wrap for this test case.
    """
    # If overlapping
    if a_start <= b_end and b_start <= a_end:
        return 0

    distances = []
    for x in [a_start, a_end]:
        for y in [b_start, b_end]:
            distances.append(circ_dist_points(x % n, y % n, n))

    return min(distances)


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
        if all(circ_dist_points(h["start_0based"], k["start_0based"], n) >= min_separation for k in kept):
            kept.append(h)

    return sorted(kept, key=lambda x: x["start_0based"])


def detect_itrs(seq, rec_id, min_score=180):
    aligner = make_aligner()
    n = len(seq)
    seq2 = seq + seq
    raw_hits = []

    for ref_name, ref in ITR_REFS.items():
        ref = clean(ref)
        aln = aligner.align(ref, seq2)[0]
        q_aligned, t_start, t_end, t_span = aln_stats(aln)

        if aln.score < min_score:
            continue

        if t_start >= n:
            continue

        raw_hits.append({
            "record_id": rec_id,
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
        delta_hits = find_delta11_near(seq, h["start_0based"], h["end_0based_exclusive"])
        has_delta = len(delta_hits) > 0

        h["itr_index"] = idx
        h["has_delta11_junction"] = has_delta
        h["delta11_junction"] = ";".join(x[0] for x in delta_hits) if delta_hits else ""
        h["delta11_position_0based"] = ";".join(str(x[1]) for x in delta_hits) if delta_hits else ""
        h["itr_class"] = classify_length(h["aligned_plasmid_span"], has_delta)

    return hits


def detect_ori(seq, min_score=150):
    aligner = make_aligner()
    n = len(seq)
    seq2 = seq + seq

    best = None

    for ref_name, ref in ORI_REFS.items():
        ref = clean(ref)
        aln = aligner.align(ref, seq2)[0]
        q_aligned, t_start, t_end, t_span = aln_stats(aln)

        if aln.score < min_score:
            continue

        if t_start >= n:
            continue

        cand = {
            "ori_label": ref_name,
            "ori_start_0based": t_start % n,
            "ori_end_0based_exclusive": t_end % n,
            "ori_score": float(aln.score),
            "ori_aligned_ref_bases": q_aligned,
            "ori_aligned_plasmid_span": t_span,
        }

        if best is None or cand["ori_score"] > best["ori_score"]:
            best = cand

    return best


def analyze_fasta(path):
    rows = []
    summary_rows = []

    for rec in SeqIO.parse(path, "fasta"):
        seq = clean(str(rec.seq))
        n = len(seq)

        itrs = detect_itrs(seq, rec.id)
        ori = detect_ori(seq)

        if ori:
            for h in itrs:
                h["distance_to_ori_bp"] = circular_interval_distance(
                    h["start_0based"],
                    h["end_0based_exclusive"],
                    ori["ori_start_0based"],
                    ori["ori_end_0based_exclusive"],
                    n
                )
                h.update(ori)
        else:
            for h in itrs:
                h["distance_to_ori_bp"] = None
                h["ori_label"] = None
                h["ori_start_0based"] = None
                h["ori_end_0based_exclusive"] = None
                h["ori_score"] = None

        if itrs:
            prox = min(
                [h for h in itrs if h["distance_to_ori_bp"] is not None],
                key=lambda h: h["distance_to_ori_bp"],
                default=None
            )
        else:
            prox = None

        length_classes = [h["itr_class"] for h in itrs]

        summary_rows.append({
            "record_id": rec.id,
            "plasmid_length": n,
            "num_itr_hits": len(itrs),
            "itr_classes": ";".join(length_classes),
            "has_130_119_pair": ("130_like" in length_classes and "119_delta11" in length_classes),
            "ori_label": ori["ori_label"] if ori else None,
            "ori_start_0based": ori["ori_start_0based"] if ori else None,
            "ori_end_0based_exclusive": ori["ori_end_0based_exclusive"] if ori else None,
            "origin_proximal_itr_index": prox["itr_index"] if prox else None,
            "origin_proximal_itr_class": prox["itr_class"] if prox else None,
            "is_119_delta11_origin_proximal": (prox["itr_class"] == "119_delta11") if prox else None,
        })

        rows.extend(itrs)

    return pd.DataFrame(summary_rows), pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("--summary-out", default="itr_ori_summary.csv")
    parser.add_argument("--details-out", default="itr_ori_details.csv")
    args = parser.parse_args()

    summary, details = analyze_fasta(args.fasta)

    print("\nSUMMARY")
    print(summary.to_string(index=False))

    print("\nDETAILS")
    print(details.to_string(index=False))

    summary.to_csv(args.summary_out, index=False)
    details.to_csv(args.details_out, index=False)

    print(f"\nWrote {args.summary_out}")
    print(f"Wrote {args.details_out}")


if __name__ == "__main__":
    main()
