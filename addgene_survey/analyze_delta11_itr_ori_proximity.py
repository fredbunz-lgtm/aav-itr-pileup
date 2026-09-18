#!/usr/bin/env python3

import csv
import os
import re
from pathlib import Path

import requests


BEST_FULL_TSV = Path("results/itr_positive_best_full_sequences.tsv")
ITR_HITS_TSV = Path("results/addgene_candidate_itr_hits.tsv")
DELTA11_PLASMID_TSV = Path("results/delta11_exact_plasmid_summary.tsv")
GENBANK_DIR = Path("genbank_files")

OUT_DETAIL_TSV = Path("results/delta11_itr_ori_proximity_detail.tsv")
OUT_SUMMARY_TSV = Path("results/delta11_itr_ori_proximity_summary.tsv")

API_BASE = "https://api.developers.addgene.org"
TOKEN_ENV = "ADDGENE_API_TOKEN"

DELTA11_FORWARD = "TGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCG"
DELTA11_REVERSE = "CGACCAAAGGTCGCCCGACGCCCGGGCGGCCTCA"


def clean_sequence(seq):
    if seq is None:
        return ""
    return re.sub(r"[^A-Za-z]", "", str(seq)).upper()


def load_delta11_positive_plasmids(path):
    ids = set()

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

        if "plasmid_id" not in fieldnames:
            raise ValueError(f"{path} missing plasmid_id")

        if "has_exact_delta11" in fieldnames:
            for row in reader:
                val = str(row["has_exact_delta11"]).strip().lower()
                if val in {"1", "true", "yes", "y"}:
                    ids.add(str(row["plasmid_id"]))
            return ids

        count_col = None
        for c in [
            "max_exact_delta11_count_per_sequence",
            "max_delta11_exact_total_count",
            "delta11_exact_total_count",
            "max_exact_delta11_count",
            "delta11_exact_count",
        ]:
            if c in fieldnames:
                count_col = c
                break

        if count_col is None:
            raise ValueError(f"No delta11 indicator/count column found in {path}")

        for row in reader:
            try:
                if int(row[count_col]) >= 1:
                    ids.add(str(row["plasmid_id"]))
            except Exception:
                pass

    return ids


def load_best_full_sequences(path):
    by_plasmid = {}

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = str(row["plasmid_id"])
            by_plasmid[pid] = row

    return by_plasmid


def load_itr_hits_for_best_sequences(path, best_by_plasmid):
    best_sequence_by_plasmid = {
        str(pid): str(row["sequence_id"])
        for pid, row in best_by_plasmid.items()
    }

    hits = {}

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = str(row["plasmid_id"])
            sid = str(row["sequence_id"])

            if pid not in best_sequence_by_plasmid:
                continue
            if sid != best_sequence_by_plasmid[pid]:
                continue

            try:
                start = int(row["start_0based"])
                end = int(row["end_0based_exclusive"])
                length = int(row["plasmid_length"])
            except Exception:
                continue

            row["_start"] = start
            row["_end"] = end
            row["_length"] = length
            row["_mid"] = circular_interval_midpoint(start, end, length)

            hits.setdefault(pid, []).append(row)

    for pid in hits:
        hits[pid].sort(key=lambda r: (r["_start"], r["_end"]))

    return hits


def circular_interval_midpoint(start, end, length):
    """
    Return midpoint coordinate for interval [start, end) on circular molecule.
    Current ITR hits should usually not wrap; handle wrap just in case.
    """
    if length <= 0:
        return None

    start %= length
    end %= length

    if start < end:
        return (start + end) / 2.0

    # wrapped interval
    span = (length - start) + end
    return (start + span / 2.0) % length


def circular_distance(a, b, length):
    if a is None or b is None or length <= 0:
        return None

    d = abs(float(a) - float(b))
    return min(d, length - d)


def genbank_path_for_sequence_id(sequence_id):
    return GENBANK_DIR / f"sequence_{sequence_id}.gb"


def download_genbank_if_needed(sequence_id):
    GENBANK_DIR.mkdir(exist_ok=True)

    path = genbank_path_for_sequence_id(sequence_id)
    if path.exists() and path.stat().st_size > 0:
        return path, "cached"

    token = os.environ.get(TOKEN_ENV)
    if not token:
        raise RuntimeError(f"{TOKEN_ENV} is not set and {path} is missing")

    url = f"{API_BASE}/download/genbank/{sequence_id}/"
    headers = {"Authorization": f"Token {token}"}

    r = requests.get(url, headers=headers, timeout=120)
    if r.status_code != 200:
        raise RuntimeError(
            f"Failed GenBank download sequence_id={sequence_id}: "
            f"HTTP {r.status_code}; {r.text[:300]}"
        )

    path.write_text(r.text)
    return path, "downloaded"


def parse_genbank_origin_sequence(gb_path):
    in_origin = False
    chunks = []

    with gb_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                chunks.append(line)

    return clean_sequence("".join(chunks))


def find_exact_motif_positions_circular(seq, motif):
    seq = clean_sequence(seq)
    motif = clean_sequence(motif)

    if not seq or not motif:
        return []

    n = len(seq)
    extended = seq + seq[: len(motif) - 1]
    out = []

    pos = 0
    while True:
        idx = extended.find(motif, pos)
        if idx == -1:
            break
        if idx < n:
            out.append(idx)
        pos = idx + 1

    return sorted(set(out))


def find_delta11_hits(seq):
    hits = []

    for p in find_exact_motif_positions_circular(seq, DELTA11_FORWARD):
        hits.append({"orientation": "forward", "position_0based": p})

    for p in find_exact_motif_positions_circular(seq, DELTA11_REVERSE):
        hits.append({"orientation": "reverse", "position_0based": p})

    hits.sort(key=lambda x: (x["position_0based"], x["orientation"]))
    return hits


def parse_location_numbers(location):
    """
    Extract coordinate numbers from GenBank feature location string.
    Converts to 0-based closed-ish points later by subtracting 1 from start.
    Handles simple cases like:
      6813..7401
      complement(6813..7401)
      join(1..20,7000..7401)
    """
    nums = [int(x) for x in re.findall(r"\d+", location)]
    return nums


def parse_rep_origin_features(gb_path):
    """
    Minimal GenBank FEATURES parser for rep_origin features.

    Returns list of dicts with:
      start_0based, end_0based_exclusive, midpoint_0based, label, note, raw_location
    """
    features = []
    current = None

    with gb_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("ORIGIN"):
                break

            # Feature line: five spaces, key, then location
            m = re.match(r"^     (\S+)\s+(.+?)\s*$", line)
            if m and not line[5:21].strip().startswith("/"):
                if current is not None:
                    features.append(current)

                key = m.group(1)
                loc = m.group(2).strip()

                if key == "rep_origin":
                    current = {
                        "key": key,
                        "raw_location": loc,
                        "label": "",
                        "note": "",
                    }
                else:
                    current = None
                continue

            if current is not None:
                q = line[21:].strip() if len(line) > 21 else line.strip()

                if q.startswith("/label="):
                    current["label"] = q.replace("/label=", "", 1).strip().strip('"')
                elif q.startswith("/note="):
                    current["note"] = q.replace("/note=", "", 1).strip().strip('"')
                elif current.get("note") and not q.startswith("/"):
                    current["note"] += " " + q.strip().strip('"')
                elif current.get("label") and not q.startswith("/"):
                    current["label"] += " " + q.strip().strip('"')

    if current is not None:
        features.append(current)

    parsed = []

    for feat in features:
        nums = parse_location_numbers(feat["raw_location"])
        if len(nums) < 2:
            continue

        start_1based = min(nums)
        end_1based = max(nums)

        start0 = start_1based - 1
        end0 = end_1based

        feat["start_0based"] = start0
        feat["end_0based_exclusive"] = end0
        parsed.append(feat)

    return parsed


def is_bacterial_plasmid_ori(feature):
    text = (feature.get("label", "") + " " + feature.get("note", "")).lower()

    # Exclude f1/phage origins when there is also a normal plasmid ori.
    if "f1" in text or "bacteriophage" in text or "phage" in text:
        return False

    bacterial_terms = [
        "cole1",
        "col e1",
        "pmb1",
        "pbr322",
        "puc",
        "rk2",
        "r6k",
        "ori",
        "origin of replication",
        "replication origin",
    ]

    return any(term in text for term in bacterial_terms)


def choose_bacterial_ori(rep_origins, plasmid_length):
    bacterial = [r for r in rep_origins if is_bacterial_plasmid_ori(r)]

    if not bacterial:
        return None

    # Prefer non-f1 bacterial ori. If multiple, pick longest feature.
    for r in bacterial:
        r["_span"] = abs(int(r["end_0based_exclusive"]) - int(r["start_0based"]))

    chosen = sorted(bacterial, key=lambda r: r["_span"], reverse=True)[0]
    chosen["midpoint_0based"] = circular_interval_midpoint(
        int(chosen["start_0based"]),
        int(chosen["end_0based_exclusive"]),
        plasmid_length,
    )
    return chosen


def position_within_or_near_itr(pos, itr, plasmid_length, padding=75):
    """
    Assign Δ11 motif to ITR if motif start is inside or near the ITR hit interval.
    """
    start = int(itr["_start"])
    end = int(itr["_end"])

    # Non-wrapping expected.
    if start <= end:
        if start - padding <= pos <= end + padding:
            return True
    else:
        # wrapped interval
        if pos >= start - padding or pos <= end + padding:
            return True

    return False


def assign_delta11_to_itr(delta_pos, itrs, plasmid_length):
    containing = [
        i for i, itr in enumerate(itrs)
        if position_within_or_near_itr(delta_pos, itr, plasmid_length, padding=75)
    ]

    if len(containing) == 1:
        return containing[0], "within_or_near"

    # Fallback: nearest ITR midpoint.
    distances = []
    for i, itr in enumerate(itrs):
        d = circular_distance(delta_pos, itr["_mid"], plasmid_length)
        distances.append((d, i))

    distances.sort()
    if len(distances) >= 2 and distances[0][0] == distances[1][0]:
        return None, "ambiguous_equal_distance"

    return distances[0][1], "nearest_midpoint"


def main():
    delta11_positive = load_delta11_positive_plasmids(DELTA11_PLASMID_TSV)
    best_by_plasmid = load_best_full_sequences(BEST_FULL_TSV)
    itr_hits_by_plasmid = load_itr_hits_for_best_sequences(ITR_HITS_TSV, best_by_plasmid)

    print(f"Delta11-positive plasmids: {len(delta11_positive)}")
    print(f"Best full sequence plasmids: {len(best_by_plasmid)}")
    print(f"Plasmids with ITR hits on best full sequence: {len(itr_hits_by_plasmid)}")

    rows = []
    counts = {}

    def inc(k):
        counts[k] = counts.get(k, 0) + 1

    analyzable_ids = sorted(delta11_positive & set(best_by_plasmid.keys()))

    for n, pid in enumerate(analyzable_ids, start=1):
        if n % 250 == 0:
            print(f"Processed {n} / {len(analyzable_ids)}")

        best = best_by_plasmid[pid]
        sid = str(best["sequence_id"])
        pname = best.get("plasmid_name", "")

        status = "started"
        reason = ""

        itrs = itr_hits_by_plasmid.get(pid, [])
        if len(itrs) != 2:
            status = "excluded"
            reason = f"best_sequence_itr_count_{len(itrs)}"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason))
            continue

        try:
            gb_path, gb_status = download_genbank_if_needed(sid)
            seq = parse_genbank_origin_sequence(gb_path)
        except Exception as e:
            status = "excluded"
            reason = "genbank_error"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason, error=repr(e)))
            continue

        plasmid_length = len(seq)
        if plasmid_length <= 0:
            status = "excluded"
            reason = "empty_genbank_sequence"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason))
            continue

        delta_hits = find_delta11_hits(seq)
        if len(delta_hits) != 1:
            status = "excluded"
            reason = f"delta11_exact_count_{len(delta_hits)}"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason))
            continue

        rep_origins = parse_rep_origin_features(gb_path)
        ori = choose_bacterial_ori(rep_origins, plasmid_length)

        if ori is None:
            status = "excluded"
            reason = "no_bacterial_ori"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason))
            continue

        delta_pos = int(delta_hits[0]["position_0based"])
        delta_orientation = delta_hits[0]["orientation"]

        delta_itr_index, assignment_method = assign_delta11_to_itr(
            delta_pos,
            itrs,
            plasmid_length,
        )

        if delta_itr_index is None:
            status = "excluded"
            reason = "ambiguous_delta11_itr_assignment"
            inc(reason)
            rows.append(make_row(pid, sid, pname, status, reason))
            continue

        ori_mid = ori["midpoint_0based"]

        itr0_dist = circular_distance(ori_mid, itrs[0]["_mid"], plasmid_length)
        itr1_dist = circular_distance(ori_mid, itrs[1]["_mid"], plasmid_length)

        if itr0_dist == itr1_dist:
            call = "tie"
        elif delta_itr_index == 0 and itr0_dist < itr1_dist:
            call = "delta11_itr_closest_to_ori"
        elif delta_itr_index == 1 and itr1_dist < itr0_dist:
            call = "delta11_itr_closest_to_ori"
        else:
            call = "delta11_itr_not_closest_to_ori"

        inc(call)

        rows.append(
            make_row(
                pid,
                sid,
                pname,
                "analyzed",
                call,
                plasmid_length=plasmid_length,
                delta11_position_0based=delta_pos,
                delta11_position_1based=delta_pos + 1,
                delta11_orientation=delta_orientation,
                delta11_itr_number=delta_itr_index + 1,
                delta11_itr_assignment_method=assignment_method,
                ori_start_0based=ori["start_0based"],
                ori_end_0based_exclusive=ori["end_0based_exclusive"],
                ori_midpoint_0based=ori_mid,
                ori_label=ori.get("label", ""),
                ori_note=ori.get("note", ""),
                itr1_start_0based=itrs[0]["_start"],
                itr1_end_0based_exclusive=itrs[0]["_end"],
                itr1_midpoint_0based=itrs[0]["_mid"],
                itr1_ori_distance=itr0_dist,
                itr2_start_0based=itrs[1]["_start"],
                itr2_end_0based_exclusive=itrs[1]["_end"],
                itr2_midpoint_0based=itrs[1]["_mid"],
                itr2_ori_distance=itr1_dist,
            )
        )

    fieldnames = [
        "plasmid_id",
        "sequence_id",
        "plasmid_name",
        "status",
        "reason_or_call",
        "error",
        "plasmid_length",
        "delta11_position_0based",
        "delta11_position_1based",
        "delta11_orientation",
        "delta11_itr_number",
        "delta11_itr_assignment_method",
        "ori_start_0based",
        "ori_end_0based_exclusive",
        "ori_midpoint_0based",
        "ori_label",
        "ori_note",
        "itr1_start_0based",
        "itr1_end_0based_exclusive",
        "itr1_midpoint_0based",
        "itr1_ori_distance",
        "itr2_start_0based",
        "itr2_end_0based_exclusive",
        "itr2_midpoint_0based",
        "itr2_ori_distance",
    ]

    OUT_DETAIL_TSV.parent.mkdir(exist_ok=True)

    with OUT_DETAIL_TSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with OUT_SUMMARY_TSV.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "count"])
        for k in sorted(counts):
            writer.writerow([k, counts[k]])

    print(f"Wrote: {OUT_DETAIL_TSV}")
    print(f"Wrote: {OUT_SUMMARY_TSV}")
    print()
    print("Summary:")
    for k in sorted(counts):
        print(f"{k}\t{counts[k]}")

    analyzed = (
        counts.get("delta11_itr_closest_to_ori", 0)
        + counts.get("delta11_itr_not_closest_to_ori", 0)
        + counts.get("tie", 0)
    )
    closest = counts.get("delta11_itr_closest_to_ori", 0)
    not_closest = counts.get("delta11_itr_not_closest_to_ori", 0)

    if analyzed:
        print()
        print(f"Analyzed clean cases: {analyzed}")
        print(f"Delta11 ITR closest to ori: {closest} / {analyzed} = {closest / analyzed:.1%}")
        print(f"Delta11 ITR not closest to ori: {not_closest} / {analyzed} = {not_closest / analyzed:.1%}")


def make_row(
    plasmid_id,
    sequence_id,
    plasmid_name,
    status,
    reason_or_call,
    error="",
    plasmid_length="",
    delta11_position_0based="",
    delta11_position_1based="",
    delta11_orientation="",
    delta11_itr_number="",
    delta11_itr_assignment_method="",
    ori_start_0based="",
    ori_end_0based_exclusive="",
    ori_midpoint_0based="",
    ori_label="",
    ori_note="",
    itr1_start_0based="",
    itr1_end_0based_exclusive="",
    itr1_midpoint_0based="",
    itr1_ori_distance="",
    itr2_start_0based="",
    itr2_end_0based_exclusive="",
    itr2_midpoint_0based="",
    itr2_ori_distance="",
):
    return {
        "plasmid_id": plasmid_id,
        "sequence_id": sequence_id,
        "plasmid_name": plasmid_name,
        "status": status,
        "reason_or_call": reason_or_call,
        "error": error,
        "plasmid_length": plasmid_length,
        "delta11_position_0based": delta11_position_0based,
        "delta11_position_1based": delta11_position_1based,
        "delta11_orientation": delta11_orientation,
        "delta11_itr_number": delta11_itr_number,
        "delta11_itr_assignment_method": delta11_itr_assignment_method,
        "ori_start_0based": ori_start_0based,
        "ori_end_0based_exclusive": ori_end_0based_exclusive,
        "ori_midpoint_0based": ori_midpoint_0based,
        "ori_label": ori_label,
        "ori_note": ori_note,
        "itr1_start_0based": itr1_start_0based,
        "itr1_end_0based_exclusive": itr1_end_0based_exclusive,
        "itr1_midpoint_0based": itr1_midpoint_0based,
        "itr1_ori_distance": itr1_ori_distance,
        "itr2_start_0based": itr2_start_0based,
        "itr2_end_0based_exclusive": itr2_end_0based_exclusive,
        "itr2_midpoint_0based": itr2_midpoint_0based,
        "itr2_ori_distance": itr2_ori_distance,
    }


if __name__ == "__main__":
    main()
