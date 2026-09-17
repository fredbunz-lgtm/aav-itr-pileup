#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

BEST_FULL_TSV = Path("results/itr_positive_best_full_sequences.tsv")
ITR_HITS_TSV = Path("results/addgene_candidate_itr_hits.tsv")
DELTA11_PLASMID_TSV = Path("results/delta11_exact_plasmid_summary.tsv")

OUT_DETAIL_TSV = Path("results/itr_length_isomer_detail.tsv")
OUT_PLASMID_TSV = Path("results/itr_length_isomer_plasmid_summary.tsv")
OUT_COUNTS_TSV = Path("results/itr_length_isomer_counts.tsv")


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def int_or_none(x):
    try:
        if x is None or x == "":
            return None
        return int(float(x))
    except Exception:
        return None


def load_best_full_sequences(path):
    rows = read_tsv(path)
    by_pid = {}
    for r in rows:
        pid = str(r.get("plasmid_id", ""))
        sid = str(r.get("sequence_id", ""))
        if pid and sid:
            by_pid[pid] = r
    return by_pid


def load_delta11_summary(path):
    rows = read_tsv(path)
    by_pid = {}
    for r in rows:
        pid = str(r.get("plasmid_id", ""))
        if pid:
            by_pid[pid] = r
    return by_pid


def infer_isomer(best_ref):
    if "Flip" in best_ref:
        return "flip"
    if "Flop" in best_ref:
        return "flop"
    return "unknown"


def infer_ref_orientation(best_ref):
    """
    Orientation of the best matching reference sequence used in the scan.
    This is distinct from biological flip/flop isomer state.
    """
    if best_ref.startswith("145RC_"):
        return "reverse_complement_reference"
    if best_ref.startswith("145_"):
        return "forward_reference"
    return "unknown"


def bool_from_text(x):
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def classify_itr_length(row):
    """
    ITR length/state classification from existing alignment fields.

    Main categories:
      145_full
      130_truncated
      119_delta11
      119_like_no_delta11_motif
      other_truncated_or_deleted
      ambiguous

    The exact Δ11 call takes precedence over length, because a 119-bp Δ11 ITR is
    mechanistically distinct from a generic 119-bp partial alignment.
    """

    ref_bases = int_or_none(row.get("aligned_ref_bases"))
    span = int_or_none(row.get("aligned_plasmid_span"))
    has_delta = bool_from_text(row.get("has_delta11_junction", ""))
    itr_class = row.get("itr_class", "")

    if ref_bases is None or span is None:
        return "ambiguous"

    if has_delta or itr_class == "119_delta11":
        return "119_delta11"

    # Strict 145-bp full-length ITR. This correctly captures pAAV2ST/Addgene 239400.
    if ref_bases >= 140 and span >= 140:
        return "145_full"

    # Common psub201-like 130-bp ITR. Use existing class or tight numerical range.
    if itr_class == "130_like" or (125 <= ref_bases <= 134 and 125 <= span <= 134):
        return "130_truncated"

    # 119-like but lacking the exact Δ11 junction motif.
    if itr_class == "119_like_no_delta11_motif" or (114 <= ref_bases <= 124 and 114 <= span <= 124):
        return "119_like_no_delta11_motif"

    return "other_truncated_or_deleted"


def sorted_pattern(values):
    values = [v for v in values if v]
    if not values:
        return "none"
    return "/".join(sorted(values))


def ordered_pattern(values):
    values = [v for v in values if v]
    if not values:
        return "none"
    return "/".join(values)


def get_delta_count(delta_row):
    if not delta_row:
        return ""
    for c in [
        "max_exact_delta11_count_per_sequence",
        "exact_delta11_count",
        "max_delta11_count",
        "delta11_count",
    ]:
        if c in delta_row and delta_row[c] != "":
            return delta_row[c]
    return ""

def main():
    best_by_pid = load_best_full_sequences(BEST_FULL_TSV)
    delta_by_pid = load_delta11_summary(DELTA11_PLASMID_TSV)

    # Choose the representative sequence for this structural analysis directly
    # from the ITR-hit table. The earlier best-full table is retained for
    # plasmid metadata, but a small number of plasmids have their ITR-containing
    # sequence in a different Addgene sequence record than the previously chosen
    # best/full record.
    all_hits_by_pid_sid = defaultdict(list)

    with open(ITR_HITS_TSV, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = str(row.get("plasmid_id", ""))
            sid = str(row.get("sequence_id", ""))
            if not pid or not sid:
                continue
            if pid not in best_by_pid:
                continue
            all_hits_by_pid_sid[(pid, sid)].append(row)

    selected_sid_by_pid = {}
    for pid in best_by_pid:
        candidates = [
            (sid, rows)
            for (p, sid), rows in all_hits_by_pid_sid.items()
            if p == pid
        ]
        if not candidates:
            selected_sid_by_pid[pid] = str(best_by_pid[pid].get("sequence_id", ""))
            continue

        # Prefer the sequence with the most ITR calls. Break ties by sequence
        # length when available, then by numeric sequence_id for determinism.
        def rank(item):
            sid, rows = item
            lengths = [int_or_none(r.get("sequence_length_observed")) for r in rows]
            max_len = max([x for x in lengths if x is not None], default=-1)
            sid_int = int_or_none(sid)
            return (len(rows), max_len, sid_int if sid_int is not None else -1)

        selected_sid_by_pid[pid] = max(candidates, key=rank)[0]

    detail_rows = []
    hits_by_pid = defaultdict(list)

    for (pid, sid), rows in all_hits_by_pid_sid.items():
        if sid != selected_sid_by_pid.get(pid):
            continue

        for row in rows:
            best_ref = row.get("best_ref", "")
            length_class = classify_itr_length(row)
            isomer = infer_isomer(best_ref)
            ref_orientation = infer_ref_orientation(best_ref)

            out = dict(row)
            out["itr_length_class"] = length_class
            out["itr_isomer"] = isomer
            out["best_ref_orientation"] = ref_orientation
            out["selected_representative_sequence_id"] = sid

            detail_rows.append(out)
            hits_by_pid[pid].append(out)

    plasmid_rows = []

    for pid, best_row in sorted(best_by_pid.items(), key=lambda x: int_or_none(x[0]) or 0):
        hits = sorted(
            hits_by_pid.get(pid, []),
            key=lambda r: int_or_none(r.get("start_0based")) if int_or_none(r.get("start_0based")) is not None else -1
        )

        length_classes_ordered = [h["itr_length_class"] for h in hits]
        length_classes_sorted = sorted(length_classes_ordered)

        isomers_ordered = [h["itr_isomer"] for h in hits]
        full145_isomers = [
            h["itr_isomer"] for h in hits
            if h["itr_length_class"] == "145_full"
        ]
        class_counts = Counter(length_classes_ordered)
        isomer_counts = Counter(isomers_ordered)

        delta_row = delta_by_pid.get(pid, {})
        exact_delta_count = get_delta_count(delta_row)

        plasmid_rows.append({
            "plasmid_id": pid,
            "sequence_id": best_row.get("sequence_id", ""),
            "plasmid_name": best_row.get("plasmid_name", ""),
            "sequence_group": best_row.get("sequence_group", ""),
            "sequence_length_observed": best_row.get("sequence_length_observed", ""),
            "best_sequence_itr_hit_count": len(hits),

            "itr_length_pattern_sorted": sorted_pattern(length_classes_sorted),
            "itr_length_pattern_ordered_by_coordinate": ordered_pattern(length_classes_ordered),
            "itr_isomer_pattern_sorted": sorted_pattern(isomers_ordered),
            "itr_isomer_pattern_ordered_by_coordinate": ordered_pattern(isomers_ordered),
            "full145_isomer_pattern_sorted": sorted_pattern(full145_isomers),

            "145_full_count": class_counts.get("145_full", 0),
            "130_truncated_count": class_counts.get("130_truncated", 0),
            "119_delta11_count": class_counts.get("119_delta11", 0),
            "119_like_no_delta11_motif_count": class_counts.get("119_like_no_delta11_motif", 0),
            "other_truncated_or_deleted_count": class_counts.get("other_truncated_or_deleted", 0),
            "ambiguous_count": class_counts.get("ambiguous", 0),

            "flip_count": isomer_counts.get("flip", 0),
            "flop_count": isomer_counts.get("flop", 0),
            "unknown_isomer_count": isomer_counts.get("unknown", 0),

            "exact_delta11_plasmid_count": exact_delta_count,
            "itr_length_class_counts": ";".join(f"{k}:{v}" for k, v in class_counts.most_common()),
            "itr_isomer_counts": ";".join(f"{k}:{v}" for k, v in isomer_counts.most_common()),
            "itr_starts_0based": ";".join(h.get("start_0based", "") for h in hits),
            "itr_best_refs": ";".join(h.get("best_ref", "") for h in hits),
            "itr_aligned_ref_bases": ";".join(h.get("aligned_ref_bases", "") for h in hits),
            "itr_aligned_plasmid_spans": ";".join(h.get("aligned_plasmid_span", "") for h in hits),
        })

    OUT_DETAIL_TSV.parent.mkdir(parents=True, exist_ok=True)

    detail_fields = list(detail_rows[0].keys()) if detail_rows else []
    if detail_fields:
        with open(OUT_DETAIL_TSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=detail_fields, delimiter="\t")
            w.writeheader()
            w.writerows(detail_rows)

    plasmid_fields = [
        "plasmid_id",
        "sequence_id",
        "plasmid_name",
        "sequence_group",
        "sequence_length_observed",
        "best_sequence_itr_hit_count",
        "itr_length_pattern_sorted",
        "itr_length_pattern_ordered_by_coordinate",
        "itr_isomer_pattern_sorted",
        "itr_isomer_pattern_ordered_by_coordinate",
        "full145_isomer_pattern_sorted",
        "145_full_count",
        "130_truncated_count",
        "119_delta11_count",
        "119_like_no_delta11_motif_count",
        "other_truncated_or_deleted_count",
        "ambiguous_count",
        "flip_count",
        "flop_count",
        "unknown_isomer_count",
        "exact_delta11_plasmid_count",
        "itr_length_class_counts",
        "itr_isomer_counts",
        "itr_starts_0based",
        "itr_best_refs",
        "itr_aligned_ref_bases",
        "itr_aligned_plasmid_spans",
    ]

    with open(OUT_PLASMID_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=plasmid_fields, delimiter="\t")
        w.writeheader()
        w.writerows(plasmid_rows)

    # Counts
    counts = []

    def add(section, label, value):
        counts.append({
            "section": section,
            "label": label,
            "count": value,
        })

    add("overall", "best_full_plasmids", len(best_by_pid))
    add("overall", "plasmids_with_itr_hits_on_best_full_sequence", len(hits_by_pid))
    add("overall", "itr_hits_on_best_full_sequences", len(detail_rows))

    per_itr_length = Counter(r["itr_length_class"] for r in detail_rows)
    per_itr_isomer = Counter(r["itr_isomer"] for r in detail_rows)
    per_itr_best_ref = Counter(r["best_ref"] for r in detail_rows)

    for k, v in per_itr_length.most_common():
        add("per_itr_length_class", k, v)

    for k, v in per_itr_isomer.most_common():
        add("per_itr_isomer_all_hits", k, v)

    for k, v in per_itr_best_ref.most_common():
        add("per_itr_best_ref", k, v)

    per_plasmid_itr_count = Counter(int(r["best_sequence_itr_hit_count"]) for r in plasmid_rows)
    per_plasmid_pattern = Counter(r["itr_length_pattern_sorted"] for r in plasmid_rows)
    per_plasmid_ordered_pattern = Counter(r["itr_length_pattern_ordered_by_coordinate"] for r in plasmid_rows)
    per_plasmid_isomer_pattern = Counter(r["itr_isomer_pattern_sorted"] for r in plasmid_rows)

    for k in sorted(per_plasmid_itr_count):
        add("per_plasmid_itr_count", str(k), per_plasmid_itr_count[k])

    for k, v in per_plasmid_pattern.most_common():
        add("per_plasmid_length_pattern_sorted_all", k, v)

    for k, v in per_plasmid_ordered_pattern.most_common():
        add("per_plasmid_length_pattern_ordered_all", k, v)

    for k, v in per_plasmid_isomer_pattern.most_common():
        add("per_plasmid_isomer_pattern_sorted_all", k, v)

    # Key two-ITR subset
    two_itr = [r for r in plasmid_rows if int(r["best_sequence_itr_hit_count"]) == 2]
    two_itr_pattern = Counter(r["itr_length_pattern_sorted"] for r in two_itr)
    two_itr_isomer_pattern = Counter(r["itr_isomer_pattern_sorted"] for r in two_itr)

    add("two_itr_subset", "two_itr_plasmids", len(two_itr))

    for k, v in two_itr_pattern.most_common():
        add("two_itr_length_pattern_sorted", k, v)

    for k, v in two_itr_isomer_pattern.most_common():
        add("two_itr_isomer_pattern_sorted", k, v)

    # Key biologically expected categories
    def n_pattern(pattern):
        return two_itr_pattern.get(pattern, 0)

    add("key_two_itr_patterns", "119_delta11/130_truncated", n_pattern("119_delta11/130_truncated"))
    add("key_two_itr_patterns", "130_truncated/130_truncated", n_pattern("130_truncated/130_truncated"))
    add("key_two_itr_patterns", "145_full/145_full", n_pattern("145_full/145_full"))
    add("key_two_itr_patterns", "119_delta11/119_delta11", n_pattern("119_delta11/119_delta11"))
    add("key_two_itr_patterns", "119_delta11/145_full", n_pattern("119_delta11/145_full"))
    add("key_two_itr_patterns", "130_truncated/145_full", n_pattern("130_truncated/145_full"))

    with open(OUT_COUNTS_TSV, "w", newline="") as f:
        fieldnames = ["section", "label", "count"]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(counts)

    print(f"Best full plasmids loaded: {len(best_by_pid):,}")
    print(f"Plasmids with ITR hits on best full sequence: {len(hits_by_pid):,}")
    print(f"ITR hits on best full sequences: {len(detail_rows):,}")
    print()

    print("Per-ITR length classes:")
    for k, v in per_itr_length.most_common():
        print(f"{k}\t{v:,}")
    print()

    print("Per-ITR isomer calls, all hits:")
    for k, v in per_itr_isomer.most_common():
        print(f"{k}\t{v:,}")
    print()

    print("Per-plasmid ITR count:")
    for k in sorted(per_plasmid_itr_count):
        print(f"{k}\t{per_plasmid_itr_count[k]:,}")
    print()

    print("Two-ITR length patterns:")
    for k, v in two_itr_pattern.most_common(30):
        print(f"{k}\t{v:,}")
    print()

    print("Two-ITR isomer patterns:")
    for k, v in two_itr_isomer_pattern.most_common(20):
        print(f"{k}\t{v:,}")
    print()

    print("Internal control Addgene 239400 / pAAV2ST:")
    control = [r for r in plasmid_rows if r["plasmid_id"] == "239400"]
    if control:
        r = control[0]
        print(f"pattern: {r['itr_length_pattern_sorted']}")
        print(f"ordered pattern: {r['itr_length_pattern_ordered_by_coordinate']}")
        print(f"isomer pattern: {r['itr_isomer_pattern_sorted']}")
        print(f"best refs: {r['itr_best_refs']}")
        print(f"aligned ref bases: {r['itr_aligned_ref_bases']}")
        print(f"aligned plasmid spans: {r['itr_aligned_plasmid_spans']}")
        print(f"exact delta11 count: {r['exact_delta11_plasmid_count']}")
    else:
        print("239400 not found")
    print()

    print(f"Wrote: {OUT_DETAIL_TSV}")
    print(f"Wrote: {OUT_PLASMID_TSV}")
    print(f"Wrote: {OUT_COUNTS_TSV}")


if __name__ == "__main__":
    main()
