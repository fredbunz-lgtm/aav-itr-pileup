#!/usr/bin/env python3

import csv
import os
import re
import subprocess
import time
from pathlib import Path

IN_TSV = "results/itr_positive_best_full_sequences.tsv"
OUT_TSV = "results/test_genbank_ori_annotation_coverage_100.tsv"
GENBANK_DIR = Path("genbank_files")
N_TEST = 100

API_BASE = "https://api.developers.addgene.org/download/genbank"

ORIGIN_FEATURE_RE = re.compile(r"^\s{5}rep_origin\s+(.+)")
QUALIFIER_RE = re.compile(r'^\s{21}/([^=]+)=(?:"(.*)"|(.*))')
CONTINUATION_RE = re.compile(r'^\s{21}(.+)')

ORI_TEXT_RE = re.compile(
    r"\b(ori|origin|replication|pmb1|puc|pbr322|cole1|col e1)\b",
    re.IGNORECASE,
)

BACTERIAL_ORI_RE = re.compile(
    r"\b(pmb1|puc|pbr322|cole1|col e1|high-copy|copy-number|bacterial|e\.?\s*coli)\b",
    re.IGNORECASE,
)

EXCLUDE_ORI_RE = re.compile(
    r"\b(f1|sv40|aav|adeno|adenovirus|bacteriophage)\b",
    re.IGNORECASE,
)

def load_first_full_records(path, n):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("is_full_sequence_group") != "1":
                continue
            if not row.get("sequence_id"):
                continue
            rows.append(row)
            if len(rows) >= n:
                break
    return rows

def download_genbank(sequence_id):
    GENBANK_DIR.mkdir(exist_ok=True)
    out = GENBANK_DIR / f"sequence_{sequence_id}.gb"

    if out.exists() and out.stat().st_size > 0:
        return out, "cached"

    token = os.environ.get("ADDGENE_API_TOKEN")
    if not token:
        raise SystemExit("ADDGENE_API_TOKEN is not set. Run: source .env")

    url = f"{API_BASE}/{sequence_id}/"

    cmd = [
        "curl",
        "-L",
        "-sS",
        "-H",
        f"Authorization: Token {token}",
        "-o",
        str(out),
        url,
    ]

    subprocess.run(cmd, check=True)

    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"Downloaded empty/missing file for sequence_id {sequence_id}")

    return out, "downloaded"

def parse_rep_origin_features(gb_path):
    features = []
    current = None
    current_qual = None

    with open(gb_path, errors="replace") as f:
        for line in f:
            if line.startswith("ORIGIN"):
                break

            m = ORIGIN_FEATURE_RE.match(line)
            if m:
                if current is not None:
                    features.append(current)
                current = {
                    "location": m.group(1).strip(),
                    "qualifiers": {},
                    "raw_lines": [line.rstrip("\n")],
                }
                current_qual = None
                continue

            if current is None:
                continue

            # New feature starts; close current rep_origin.
            if re.match(r"^\s{5}\S", line) and not line.startswith("                     "):
                features.append(current)
                current = None
                current_qual = None
                continue

            current["raw_lines"].append(line.rstrip("\n"))

            qm = QUALIFIER_RE.match(line)
            if qm:
                key = qm.group(1)
                val = qm.group(2) if qm.group(2) is not None else qm.group(3)
                val = val or ""
                current["qualifiers"].setdefault(key, [])
                current["qualifiers"][key].append(val)
                current_qual = key
                continue

            cm = CONTINUATION_RE.match(line)
            if cm and current_qual:
                extra = cm.group(1).strip()
                if current["qualifiers"].get(current_qual):
                    current["qualifiers"][current_qual][-1] += " " + extra.strip('"')

    if current is not None:
        features.append(current)

    return features

def feature_text(feature):
    parts = [feature.get("location", "")]
    for key, vals in feature.get("qualifiers", {}).items():
        for val in vals:
            parts.append(f"{key}={val}")
    return " ".join(parts)

def classify_ori_feature(feature):
    txt = feature_text(feature)

    has_ori_text = bool(ORI_TEXT_RE.search(txt))
    has_bacterial = bool(BACTERIAL_ORI_RE.search(txt))
    has_exclude = bool(EXCLUDE_ORI_RE.search(txt))

    if has_bacterial:
        return "bacterial_plasmid_ori"
    if has_ori_text and not has_exclude:
        return "possible_plasmid_ori"
    if has_ori_text:
        return "non_bacterial_or_uncertain_ori"
    return "rep_origin_no_ori_text"

def main():
    records = load_first_full_records(IN_TSV, N_TEST)
    print(f"Testing first {len(records)} full sequence records")

    out_fields = [
        "plasmid_id",
        "plasmid_name",
        "plasmid_url",
        "sequence_id",
        "sequence_length",
        "download_status",
        "genbank_path",
        "rep_origin_count",
        "bacterial_plasmid_ori_count",
        "possible_plasmid_ori_count",
        "non_bacterial_or_uncertain_ori_count",
        "ori_feature_summaries",
    ]

    n_with_rep_origin = 0
    n_with_bacterial = 0
    n_with_possible_or_bacterial = 0

    with open(OUT_TSV, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=out_fields)
        writer.writeheader()

        for i, row in enumerate(records, start=1):
            sid = row["sequence_id"]
            gb_path, status = download_genbank(sid)

            features = parse_rep_origin_features(gb_path)
            classes = [classify_ori_feature(feat) for feat in features]

            bacterial_count = sum(1 for c in classes if c == "bacterial_plasmid_ori")
            possible_count = sum(1 for c in classes if c == "possible_plasmid_ori")
            uncertain_count = sum(1 for c in classes if c == "non_bacterial_or_uncertain_ori")

            if features:
                n_with_rep_origin += 1
            if bacterial_count:
                n_with_bacterial += 1
            if bacterial_count or possible_count:
                n_with_possible_or_bacterial += 1

            summaries = []
            for feat, cls in zip(features, classes):
                txt = feature_text(feat)
                txt = " ".join(txt.split())
                if len(txt) > 300:
                    txt = txt[:300] + "..."
                summaries.append(f"{cls}: {txt}")

            writer.writerow({
                "plasmid_id": row["plasmid_id"],
                "plasmid_name": row["plasmid_name"],
                "plasmid_url": row["plasmid_url"],
                "sequence_id": sid,
                "sequence_length": row["sequence_length"],
                "download_status": status,
                "genbank_path": str(gb_path),
                "rep_origin_count": len(features),
                "bacterial_plasmid_ori_count": bacterial_count,
                "possible_plasmid_ori_count": possible_count,
                "non_bacterial_or_uncertain_ori_count": uncertain_count,
                "ori_feature_summaries": " || ".join(summaries),
            })

            if i % 10 == 0:
                print(f"Processed {i}/{len(records)}")

            if status == "downloaded":
                time.sleep(0.1)

    def pct(n, d):
        return "NA" if d == 0 else f"{100.0 * n / d:.1f}%"

    print()
    print(f"Records tested: {len(records)}")
    print(f"With any rep_origin feature: {n_with_rep_origin} / {len(records)} = {pct(n_with_rep_origin, len(records))}")
    print(f"With bacterial plasmid ori feature: {n_with_bacterial} / {len(records)} = {pct(n_with_bacterial, len(records))}")
    print(f"With bacterial or possible plasmid ori feature: {n_with_possible_or_bacterial} / {len(records)} = {pct(n_with_possible_or_bacterial, len(records))}")
    print(f"Wrote: {OUT_TSV}")

if __name__ == "__main__":
    main()
