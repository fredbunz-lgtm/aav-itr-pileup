#!/usr/bin/env python3
"""
AAV Genome QC Pipeline v3
==========================
Analyzes Oxford Nanopore sequencing data from encapsidated AAV viral genomes.

Key features vs wf-aav-qc:
  - Automatic ITR detection from the cis plasmid FASTA (no BED file required)
  - Full cis plasmid used as reference (detects backbone/reverse packaging)
  - Contamination analysis against helper, rep-cap, and host references
  - Strand-specific coverage across the transgene insert
  - All outputs in CSV format suitable for publication and downstream analysis
  - Publication-quality figures

Usage:
    python aav_genome_qc.py \\
        --fastq-dir /path/to/barcode03/ \\
        --cis-plasmid pAAV2ST-Pcsk9sg1_full_plasmid.fasta \\
        --rep-cap rep_cap.fasta \\
        --name pAAV2ST-Pcsk9sg1 \\
        --output-dir ./results/

    # With custom helper/host references:
    python aav_genome_qc.py \\
        --fastq-dir /path/to/barcode03/ \\
        --cis-plasmid pAAV2ST-Pcsk9sg1_full_plasmid.fasta \\
        --rep-cap rep_cap.fasta \\
        --helper my_helper.fasta \\
        --host my_host.fasta \\
        --name pAAV2ST-Pcsk9sg1 \\
        --output-dir ./results/

Arguments:
    --fastq-dir     Directory containing fastq.gz files for one sample
    --cis-plasmid   Full cis plasmid FASTA (ITRs auto-detected)
    --rep-cap       Rep-Cap plasmid FASTA (required, varies by serotype)
    --helper        Helper plasmid FASTA (default: bundled pHelper)
    --host          Host genome FASTA (default: bundled E. coli K-12)
    --serotype      AAV serotype for ITR detection (default: AAV2)
    --name          Sample name for output files
    --output-dir    Output directory (created if absent)
    --itr-tol       ITR boundary tolerance in bp (default: 100)
    --min-mapq      Minimum mapping quality (default: 20)
    --no-contamination  Skip contamination analysis
    --validate-tsv  Optional: wf-aav-qc per-read TSV for validation

Outputs:
    {name}_detected_itrs.bed        Auto-detected ITR coordinates
    {name}_contamination.csv        Read mapping to each reference
    {name}_genome_types.csv         Structural classification counts
    {name}_pileup.csv               Per-position coverage (with strand)
    {name}_itr_coverage.csv         ITR coverage statistics
    {name}_truncations.csv          Per-read alignment coordinates
    {name}_summary.csv              Single-row QC summary
    {name}_structures.png           Genome structure figure
    {name}_coverage.png             Coverage and base composition
    {name}_strand_coverage.png      Strand-specific coverage
    {name}_contamination.png        Contamination breakdown
    {name}_truncations.png          Truncation hotspot
    {name}_itr_detail.png           ITR coverage detail

NOTE: Flip/flop isomer quantification is not performed on packaged viral
genomes. The ITR hairpin structure causes insufficient and unreliable
coverage at inner loop positions in single-stranded AAV genomes sequenced
by ONT. Flip/flop analysis should be performed on the plasmid stock using
aav_plasmid_workflow.py.
"""

import argparse
import collections
import gzip
import os
import sys
import time
from pathlib import Path

import mappy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────
BASES    = 'ATGC'
ITR_TOL  = 100
MIN_MAPQ = 20

# Bundled reference paths (relative to script location)
SCRIPT_DIR   = Path(__file__).parent
REF_DIR      = SCRIPT_DIR / 'references'
DEFAULT_HELPER = REF_DIR / 'ref_helper_pHelper.fasta'
DEFAULT_HOST   = REF_DIR / 'ref_host_ecoli_K12.fasta'

# AAV2 ITR detection signatures
ITR_SIGNATURES_AAV2 = [
    'CCTCTCTGCGCGCTCGCTCGCTCACTGAGGCC',  # outer arm core
    'AGGAACCCCTAGTGATGGAG',               # D-sequence
    'GGCCACTCCCTCTCTGCGCG',               # arm start
]

D_SEQ_PATTERNS = [
    'AGGAACCCCTAGTGATGGAG',
    'CTCCATCACTAGGGGTTCCT',
    'GGAACCCCTAGTGATGGAGT',
    'ACTCCATCACTAGGGGTTCC',
]

# ── Utilities ──────────────────────────────────────────────────────────────────
def rc(s):
    return s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def load_fasta(fasta_path):
    """Load single-contig FASTA, replace IUPAC ambiguity codes."""
    lines = Path(fasta_path).read_text(errors='replace').strip().splitlines()
    seq = ''.join(l for l in lines if not l.startswith('>')).upper()
    for bad in 'MRWSYKVHDBN':
        if bad in seq:
            print(f'  WARNING: replacing ambiguous base {bad} with A in {fasta_path}')
            seq = seq.replace(bad, 'A')
    return seq


def read_fastq_gz_dir(fastq_dir):
    """Stream reads from all fastq/fastq.gz files in a directory."""
    fastq_dir = Path(fastq_dir)
    files = sorted(fastq_dir.glob('*.fastq.gz')) + \
            sorted(fastq_dir.glob('*.fastq'))
    if not files:
        raise FileNotFoundError(f'No fastq/fastq.gz files in {fastq_dir}')
    for fpath in files:
        opener = gzip.open if str(fpath).endswith('.gz') else open
        with opener(fpath, 'rt') as f:
            while True:
                h = f.readline().strip()
                if not h:
                    break
                seq  = f.readline().strip().upper()
                f.readline()  # + strand indicator
                qual = f.readline().strip()
                if seq:
                    yield h[1:].split()[0], seq, qual


# ── ITR auto-detection ─────────────────────────────────────────────────────────
def find_best_palindrome(s, min_arm=30, max_arm=45, max_mm=3):
    best = (0, -1, -1, 999)
    for armlen in range(max_arm, min_arm - 1, -1):
        for sa in range(0, len(s) - 2 * armlen):
            arm_a    = s[sa:sa + armlen]
            arm_a_rc = rc(arm_a)
            for sb in range(sa + armlen + 1, len(s) - armlen + 1):
                arm_b = s[sb:sb + armlen]
                if len(arm_b) < armlen:
                    continue
                mm = sum(x != y for x, y in zip(arm_a_rc, arm_b))
                if mm <= max_mm and armlen > best[0]:
                    best = (armlen, sa, sb, mm)
        if best[0] >= armlen:
            break
    return best


def extend_to_d_sequence(seq, palindrome_start, palindrome_end, window=60):
    """Extend palindrome boundaries to include adjacent D-sequence."""
    R = len(seq)
    true_start = palindrome_start
    true_end   = palindrome_end
    # Search before palindrome
    region_before = seq[max(0, palindrome_start - window):palindrome_start + 20]
    for pat in D_SEQ_PATTERNS:
        if pat in region_before:
            idx     = region_before.find(pat)
            abs_pos = max(0, palindrome_start - window) + idx
            if abs_pos < palindrome_start:
                true_start = min(true_start, abs_pos)
    # Search after palindrome
    region_after = seq[palindrome_end - 20:min(R, palindrome_end + window)]
    for pat in D_SEQ_PATTERNS:
        if pat in region_after:
            idx     = region_after.find(pat)
            abs_pos = (palindrome_end - 20) + idx + len(pat)
            if abs_pos > palindrome_end:
                true_end = max(true_end, abs_pos)
    return true_start, true_end


def auto_detect_itrs(seq, serotype='AAV2', verbose=True):
    """
    Automatically detect ITR positions in a plasmid FASTA.
    Returns list of dicts with ITR coordinates and structural information.
    Raises ValueError if fewer than 2 ITRs are detected.
    """
    R    = len(seq)
    sigs = ITR_SIGNATURES_AAV2  # extend for other serotypes as needed

    # Find all signature hits
    all_hits = []
    for sig in sigs:
        for strand, target in [('+', seq), ('-', rc(seq))]:
            pos = 0
            while True:
                idx = target.find(sig, pos)
                if idx == -1:
                    break
                fwd_pos = idx if strand == '+' else R - idx - len(sig)
                all_hits.append(fwd_pos)
                pos = idx + 1

    if not all_hits:
        raise ValueError(
            f'No ITR signatures found. Is this an AAV2-derived vector? '
            f'Serotype specified: {serotype}')

    # Cluster hits into ITR groups
    all_hits = sorted(set(all_hits))
    groups   = []
    current  = [all_hits[0]]
    for h in all_hits[1:]:
        if h - current[-1] < 200:
            current.append(h)
        else:
            groups.append(current)
            current = [h]
    groups.append(current)

    if verbose:
        print(f'  Found {len(groups)} ITR candidate region(s)')

    # Refine each group using palindrome search + D-sequence extension
    itr_regions = []
    for group in groups:
        window_start = max(0, min(group) - 30)
        window_end   = min(R, max(group) + 200)
        window_seq   = seq[window_start:window_end]

        armlen, sa, sb, mm = find_best_palindrome(window_seq)
        if armlen == 0:
            continue

        pal_start_abs = window_start + sa
        pal_end_abs   = window_start + sb + armlen
        inner_seq     = window_seq[sa + armlen:sb]
        true_start, true_end = extend_to_d_sequence(seq, pal_start_abs, pal_end_abs)

        itr_regions.append({
            'start_0based':     true_start,
            'end_0based':       true_end,
            'start_1based':     true_start + 1,
            'end_1based':       true_end,
            'length':           true_end - true_start,
            'palindrome_start': pal_start_abs + 1,
            'palindrome_end':   pal_end_abs,
            'arm_len':          armlen,
            'arm_mm':           mm,
            'inner_loop_flip':  inner_seq,
            'inner_loop_flop':  rc(inner_seq),
            'inner_len':        len(inner_seq),
            'n_hits':           len(group),
        })

    itr_regions.sort(key=lambda x: x['start_0based'])

    if len(itr_regions) < 2:
        raise ValueError(
            f'Expected 2 ITRs but detected {len(itr_regions)}. '
            f'Check that the cis plasmid contains both ITRs.')

    return itr_regions


# ── pAAV-MCS deletion warning ──────────────────────────────────────────────────
# The plasmid pAAV-MCS and all derivatives contain an 11-bp deletion
# (AAAGCCCGGGC) within the inner loop of one ITR, resulting in a 24-bp
# inner loop instead of the canonical 31 bp.
DELETED_INNER_LOOP    = 'CGGGCGTCGGGCGACCTTTGGTCG'   # 24 bp deleted variant
DELETED_INNER_LOOP_RC = 'CGACCAAAGGTCGCCCGACGCCCG'   # RC of above
DELETED_SEQUENCE      = 'AAAGCCCGGGC'                 # the 11 bp absent


def check_deleted_itr(itr_region, itr_number):
    """
    Check if a detected ITR matches the pAAV-MCS 11-bp deletion variant.
    Issues a warning if detected. Analysis continues normally regardless.
    """
    inner    = itr_region.get('inner_loop_flip', '')
    inner_rc = itr_region.get('inner_loop_flop', '')

    is_deleted = (
        (DELETED_INNER_LOOP in inner or DELETED_INNER_LOOP_RC in inner or
         DELETED_INNER_LOOP in inner_rc or DELETED_INNER_LOOP_RC in inner_rc)
        and DELETED_SEQUENCE not in inner
        and DELETED_SEQUENCE not in inner_rc
    )

    if is_deleted:
        print(f'\n  *** WARNING: ITR{itr_number} appears to carry the pAAV-MCS '
              f'11-bp deletion ***')
        print(f'      Detected inner loop ({len(inner)} bp): {inner}')
        print(f'      Canonical inner loop (31 bp) should contain: {DELETED_SEQUENCE}')
        print(f'      pAAV-MCS and its derivatives (pAAV-CMV-MCS, pAAV-EF1a-MCS,')
        print(f'      etc.) carry this deletion in one ITR. It may reduce packaging')
        print(f'      efficiency and ITR functionality. Analysis continues normally.\n')
        return True
    return False


def write_bed(itr_regions, plasmid_name, out_path):
    """Write auto-detected ITR coordinates as BED file."""
    lines = []
    for i, r in enumerate(itr_regions[:2]):
        lines.append(
            f"{plasmid_name}\t{r['start_0based']}\t{r['end_0based']}\tITR{i+1}")
    Path(out_path).write_text('\n'.join(lines) + '\n')
    print(f'  Saved: {out_path}')


# ── Contamination analysis ─────────────────────────────────────────────────────
def build_contamination_index(cis_seq, helper_path, host_path, rep_cap_path,
                               min_mapq=MIN_MAPQ, verbose=True):
    """
    Build mappy aligners for all reference sequences.
    Returns dict of {label: aligner}.
    """
    aligners = {}

    # Cis plasmid (transgene)
    if verbose:
        print(f'  Building index: cis plasmid ({len(cis_seq):,} bp)')
    aligners['transgene'] = mappy.Aligner(seq=cis_seq, preset='map-ont', best_n=1)

    # Helper plasmid
    if helper_path and Path(helper_path).exists():
        helper_seq = load_fasta(helper_path)
        if verbose:
            print(f'  Building index: helper plasmid ({len(helper_seq):,} bp)')
        aligners['helper'] = mappy.Aligner(seq=helper_seq, preset='map-ont', best_n=1)
    else:
        if verbose:
            print(f'  WARNING: helper reference not found at {helper_path}')
            print(f'           Skipping helper contamination analysis')

    # Host genome
    if host_path and Path(host_path).exists():
        if verbose:
            host_size = Path(host_path).stat().st_size / 1e6
            print(f'  Building index: host genome ({host_size:.1f} MB) — may take 30s...')
        aligners['host'] = mappy.Aligner(str(host_path), preset='map-ont', best_n=1)
    else:
        if verbose:
            print(f'  WARNING: host reference not found at {host_path}')
            print(f'           Skipping host contamination analysis')

    # Rep-Cap plasmid
    if rep_cap_path and Path(rep_cap_path).exists():
        rep_cap_seq = load_fasta(rep_cap_path)
        if verbose:
            print(f'  Building index: rep-cap plasmid ({len(rep_cap_seq):,} bp)')
        aligners['rep_cap'] = mappy.Aligner(seq=rep_cap_seq, preset='map-ont', best_n=1)
    else:
        if verbose:
            print(f'  WARNING: rep-cap reference not found at {rep_cap_path}')
            print(f'           Skipping rep-cap contamination analysis')

    return aligners


def assign_read_to_reference(seq, aligners, min_mapq=MIN_MAPQ):
    """
    Align a read against all references and assign to best match.
    Returns (label, score) where label is the reference category.
    """
    best_label = 'unmapped'
    best_score = -1

    for label, aligner in aligners.items():
        if aligner is None:
            continue
        for hit in aligner.map(seq):
            if hit.mapq >= min_mapq and hit.blen > best_score:
                best_score = hit.blen
                best_label = label

    return best_label, best_score


# ── Genome structure classification ───────────────────────────────────────────
def classify_read(hit, R, itr1_start, itr1_end, itr2_start, itr2_end,
                  tol=ITR_TOL):
    if hit is None:
        return 'Unmapped', 'Unmapped', -1, -1

    r_st  = hit.r_st
    r_en  = hit.r_en
    q_len = hit.blen + hit.q_st

    covers_itr1    = r_st <= itr1_start + tol
    covers_itr2    = r_en >= itr2_end   - tol
    starts_in_itr1 = r_st <= itr1_end   + tol
    ends_in_itr2   = r_en >= itr2_start - tol
    insert_len     = itr2_end - itr1_start

    is_sc  = q_len > insert_len * 1.7
    is_sbg = insert_len * 0.45 < q_len < insert_len * 0.85

    if is_sc and covers_itr1 and covers_itr2:
        return 'Full scAAV', 'Full scAAV', r_st, r_en
    if covers_itr1 and covers_itr2:
        return 'Full ssAAV', 'Full ssAAV', r_st, r_en
    if starts_in_itr1 and not ends_in_itr2:
        if is_sbg:
            sym = 'symmetric' if abs(r_st - itr1_start) < tol * 2 else 'asymmetric'
            return 'Partial scAAV', f"SBG 5` {sym}", r_st, r_en
        return 'Partial ssAAV', "5` ICG", r_st, r_en
    if not starts_in_itr1 and ends_in_itr2:
        if is_sbg:
            sym = 'symmetric' if abs(r_en - itr2_end) < tol * 2 else 'asymmetric'
            return 'Partial scAAV', f"SBG 3` {sym}", r_st, r_en
        return 'Partial ssAAV', "3` ICG", r_st, r_en
    if not starts_in_itr1 and not ends_in_itr2:
        if r_en - r_st < 100:
            return 'ITR region only', 'ITR single strand', r_st, r_en
        return 'Partial ssAAV', 'Partial ICG - no ITRs', r_st, r_en

    return 'Transgene unclassified', 'Transgene unclassified', r_st, r_en


# ── Main analysis loop ─────────────────────────────────────────────────────────
def run_analysis(fastq_dir, cis_seq, R, itr1, itr2, aligners,
                 tol=ITR_TOL, min_mapq=MIN_MAPQ,
                 min_base_qual=20, verbose=True):
    """
    Single pass through all reads:
    1. Assign each read to best reference (contamination)
    2. For transgene reads: classify genome structure + accumulate pileup
    """
    itr1_start, itr1_end = itr1['start_0based'], itr1['end_0based']
    itr2_start, itr2_end = itr2['start_0based'], itr2['end_0based']

    counts     = {b: np.zeros(R, dtype=np.int32) for b in BASES}
    counts_fwd = {b: np.zeros(R, dtype=np.int32) for b in BASES}
    counts_rev = {b: np.zeros(R, dtype=np.int32) for b in BASES}
    coverage   = np.zeros(R, dtype=np.int32)
    cov_fwd    = np.zeros(R, dtype=np.int32)
    cov_rev    = np.zeros(R, dtype=np.int32)

    contam_counts = collections.Counter()
    records       = []
    gtype_counts  = collections.Counter()

    total = 0
    t0    = time.time()

    for read_id, seq, qual in read_fastq_gz_dir(fastq_dir):
        total += 1
        if verbose and total % 10000 == 0:
            print(f'  Processed {total:,} reads ({time.time()-t0:.1f}s)', flush=True)

        # Step 1: assign to best reference
        ref_label, ref_score = assign_read_to_reference(seq, aligners, min_mapq)
        contam_counts[ref_label] += 1

        if ref_label != 'transgene':
            records.append({
                'read_id':        read_id,
                'read_len':       len(seq),
                'reference':      ref_label,
                'genome_type':    ref_label,
                'genome_subtype': ref_label,
                'r_start':        -1,
                'r_end':          -1,
                'mapq':           0,
                'strand':         0,
            })
            continue

        # Step 2: structural classification (transgene reads only)
        best_hit   = None
        best_score = -1
        for hit in aligners['transgene'].map(seq):
            if hit.mapq >= min_mapq and hit.blen > best_score:
                best_score = hit.blen
                best_hit   = hit

        if best_hit is None:
            gtype_counts['Unmapped'] += 1
            records.append({
                'read_id': read_id, 'read_len': len(seq),
                'reference': 'transgene', 'genome_type': 'Unmapped',
                'genome_subtype': 'Unmapped', 'r_start': -1,
                'r_end': -1, 'mapq': 0, 'strand': 0,
            })
            continue

        gtype, gsubtype, r_st, r_en = classify_read(
            best_hit, R, itr1_start, itr1_end, itr2_start, itr2_end, tol)
        gtype_counts[gtype] += 1

        strand = best_hit.strand
        records.append({
            'read_id': read_id, 'read_len': len(seq),
            'reference': 'transgene', 'genome_type': gtype,
            'genome_subtype': gsubtype, 'r_start': r_st,
            'r_end': r_en, 'mapq': best_hit.mapq, 'strand': strand,
        })

        # Step 3: pileup accumulation
        # Only accumulate bases within the ITR-to-ITR insert region.
        # Reads mapping entirely outside this window are counted as
        # reverse_packaging in the contamination stats but excluded from pileup.
        insert_start = itr1_start
        insert_end   = itr2_end

        # Check if this read's alignment overlaps the insert region
        if best_hit.r_en <= insert_start or best_hit.r_st >= insert_end:
            # Read maps entirely outside the insert — reverse packaging
            contam_counts['reverse_packaging'] += 1
            # Update the record to flag it
            records[-1]['genome_type']    = 'Reverse packaging'
            records[-1]['genome_subtype'] = 'Reverse packaging'
            gtype_counts['Reverse packaging'] += 1
            gtype_counts[gtype] -= 1  # undo the earlier count
            continue

        # For reverse strand hits, mappy reports q_st/q_en from the 3' end
        # of the original read. We must use rc(seq) and start at
        # len(seq) - hit.q_en, not hit.q_st.
        if strand == 1:
            query = seq
            q_pos = best_hit.q_st
        else:
            query = rc(seq)
            q_pos = len(seq) - best_hit.q_en
        r_pos = best_hit.r_st

        for length, op in best_hit.cigar:   # NOTE: mappy returns (length, op)
            if op in (0, 7, 8):             # M / = / X
                for i in range(length):
                    rp = r_pos + i
                    qp = q_pos + i
                    # Only record bases within insert region
                    if insert_start <= rp < insert_end and qp < len(query):
                        # Apply per-base quality filter
                        bq = ord(qual[qp]) - 33 if qp < len(qual) else 0
                        if bq < min_base_qual:
                            continue
                        base = query[qp]
                        if base in counts:
                            counts[base][rp] += 1
                            coverage[rp]     += 1
                            if strand == 1:
                                cov_fwd[rp]        += 1
                                counts_fwd[base][rp] += 1
                            else:
                                cov_rev[rp]        += 1
                                counts_rev[base][rp] += 1
                r_pos += length
                q_pos += length
            elif op == 1:                   # I
                q_pos += length
            elif op in (2, 3):              # D / N
                r_pos += length
            elif op in (4, 5):              # S / H
                q_pos += length

    if verbose:
        print(f'\n  Total reads: {total:,}')
        print(f'\n  Reference assignment:')
        for ref, count in sorted(contam_counts.items(), key=lambda x: -x[1]):
            print(f'    {ref}: {count:,} ({100*count/max(total,1):.1f}%)')
        transgene_total = contam_counts.get('transgene', 0)
        print(f'\n  Transgene genome type breakdown:')
        for gtype, count in sorted(gtype_counts.items(), key=lambda x: -x[1]):
            print(f'    {gtype}: {count:,} ({100*count/max(transgene_total,1):.1f}%)')
        print(f'\n  Coverage: mean={coverage.mean():.0f}x  min={coverage.min()}x')

    return (coverage, counts, counts_fwd, counts_rev, cov_fwd, cov_rev,
            pd.DataFrame(records), contam_counts, total)


# ── ITR coverage analysis ──────────────────────────────────────────────────────
def analyze_itr_coverage(coverage, counts_fwd, counts_rev, cov_fwd, cov_rev, cis_seq, itr1, itr2, R, min_cov=100):
    print(f'  NOTE: Flip/flop analysis not performed on packaged genomes.')
    print(f'        ITR hairpin structure causes unreliable coverage at inner')
    print(f'        loop positions. Use aav_plasmid_workflow.py for flip/flop.')
    rows = []
    for label, itr in [('ITR1', itr1), ('ITR2', itr2)]:
        s, e    = itr['start_0based'], itr['end_0based']
        cov     = coverage[s:e]
        fwd     = cov_fwd[s:e]
        rev     = cov_rev[s:e]
        n_ok    = (cov >= min_cov).sum()
        fwd_m   = float(fwd.mean())
        rev_m   = float(rev.mean())
        rows.append({
            'itr':                  label,
            'start_1based':         s + 1,
            'end_1based':           e,
            'length_bp':            e - s,
            'inner_loop_seq':       itr['inner_loop_flip'],
            'inner_loop_len':       itr['inner_len'],
            'mean_coverage':        round(float(cov.mean()), 1),
            'min_coverage':         int(cov.min()),
            'max_coverage':         int(cov.max()),
            'pct_above_mincov':     round(100 * n_ok / max(len(cov), 1), 1),
            'mean_fwd_coverage':    round(fwd_m, 1),
            'mean_rev_coverage':    round(rev_m, 1),
            'strand_ratio_fwd_rev': round(fwd_m / max(rev_m, 0.1), 2),
        })
        print(f'  {label}: mean={cov.mean():.0f}x  min={cov.min()}x  '
              f'fwd/rev={rows[-1]["strand_ratio_fwd_rev"]:.2f}  '
              f'>={min_cov}x: {n_ok}/{len(cov)} ({rows[-1]["pct_above_mincov"]:.0f}%)')
    return pd.DataFrame(rows)


# ── Figures ────────────────────────────────────────────────────────────────────
BG        = 'white'; PANEL = '#f8f8f8'; GRID = '#dddddd'; TEXT = 'black'
BCOL      = {'A': '#00aa00', 'T': '#dd0000', 'G': '#111111', 'C': '#0066cc'}
DEV_COLOR = '#4488cc'; FLAG_COLOR = '#cc2200'
GTYPE_COLORS = {
    'Full ssAAV':             '#2a9d3f',
    'Full scAAV':             '#1a6fcc',
    'Partial ssAAV':          '#cc4422',
    'Partial scAAV':          '#dd8800',
    'Complex':                '#8855bb',
    'ITR region only':        '#888888',
    'Transgene unclassified': '#aaaaaa',
    'Unmapped':               '#cccccc',
}
CONTAM_COLORS = {
    'transgene': '#2a9d3f',
    'helper':    '#dd8800',
    'rep_cap':   '#8855bb',
    'host':      '#cc4422',
    'unmapped':  '#aaaaaa',
}


def _style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TEXT, labelsize=8)
    for lbl in (ax.yaxis.label, ax.xaxis.label, ax.title):
        lbl.set_color(TEXT)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.grid(axis='y', color=GRID, lw=0.5, alpha=0.6)


def _add_itr_spans(ax, itr1, itr2):
    ax.axvspan(itr1['start_1based'], itr1['end_1based'],
               color='#dd8800', alpha=0.25, label='ITR1')
    ax.axvspan(itr2['start_1based'], itr2['end_1based'],
               color='#2a9d3f', alpha=0.25, label='ITR2')


def plot_contamination(contam_counts, total, name, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — Read Assignment ({total:,} reads)',
                 color=TEXT, fontsize=13, fontweight='bold')

    labels = [k for k in ['transgene', 'helper', 'rep_cap', 'host', 'unmapped']
              if k in contam_counts]
    counts = [contam_counts[k] for k in labels]
    colors = [CONTAM_COLORS.get(k, '#8b949e') for k in labels]
    pcts   = [100 * c / max(total, 1) for c in counts]

    # Bar chart
    ax = axes[0]
    bars = ax.barh(range(len(labels)), counts, color=colors, alpha=0.9, height=0.6)
    for i, (count, pct) in enumerate(zip(counts, pcts)):
        ax.text(count + total * 0.005, i,
                f'{count:,} ({pct:.1f}%)', va='center', color=TEXT, fontsize=9)
    nice_labels = {'transgene': 'Transgene (cis plasmid)',
                   'helper': 'Helper plasmid', 'rep_cap': 'Rep-Cap plasmid',
                   'host': 'Host (E. coli)', 'unmapped': 'Unmapped'}
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([nice_labels.get(l, l) for l in labels],
                        color=TEXT, fontsize=9)
    ax.set_xlabel('Read count', color=TEXT)
    ax.set_title('Read assignment', color=TEXT, fontweight='bold')
    ax.set_xlim(0, total * 1.25)
    _style(ax)

    # Pie chart
    ax = axes[1]
    ax.set_facecolor(PANEL)
    wedges, texts = ax.pie(
        counts, labels=None, colors=colors,
        startangle=90,
        wedgeprops=dict(linewidth=1, edgecolor='white'))
    # Add percentage labels outside the pie with leader lines
    for i, (wedge, count) in enumerate(zip(wedges, counts)):
        pct = 100 * count / max(total, 1)
        if pct < 0.05:
            continue
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = 1.25 * np.cos(np.radians(angle))
        y = 1.25 * np.sin(np.radians(angle))
        ha = 'left' if x > 0 else 'right'
        ax.annotate(f'{pct:.1f}%',
                    xy=(np.cos(np.radians(angle)) * 0.9,
                        np.sin(np.radians(angle)) * 0.9),
                    xytext=(x, y),
                    ha=ha, va='center', fontsize=8, color=TEXT,
                    arrowprops=dict(arrowstyle='-', color='grey', lw=0.8))
    ax.legend([nice_labels.get(l, l) for l in labels],
              loc='lower center', bbox_to_anchor=(0.5, -0.18),
              fontsize=8, facecolor=PANEL, labelcolor=TEXT, ncol=2)
    ax.set_title('Read distribution', color=TEXT, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_structures(reads_df, name, out_path):
    transgene_df = reads_df[reads_df['reference'] == 'transgene']
    if len(transgene_df) == 0:
        print('  No transgene reads for structure plot')
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — AAV Genome Structure ({len(transgene_df):,} transgene reads)',
                 color=TEXT, fontsize=13, fontweight='bold')

    for ax, col, title, excl in [
        (axes[0], 'genome_type', 'Genome Types', ['Unmapped']),
        (axes[1], 'genome_subtype', 'Genome Subtypes (top 15)',
         ['Unmapped', 'Transgene unclassified']),
    ]:
        vc = transgene_df[~transgene_df[col].isin(excl)][col].value_counts()
        if col == 'genome_subtype':
            vc = vc.head(15)
        total = vc.sum()
        if col == 'genome_type':
            colors = [GTYPE_COLORS.get(t, '#8b949e') for t in vc.index]
        else:
            colors = []
            for sub in vc.index:
                parent = transgene_df[transgene_df['genome_subtype'] == sub
                                      ]['genome_type'].iloc[0]
                colors.append(GTYPE_COLORS.get(parent, '#8b949e'))
        ax.barh(range(len(vc)), vc.values, color=colors, alpha=0.9, height=0.6)
        for i, count in enumerate(vc.values):
            ax.text(count + total * 0.01, i,
                    f'{count:,} ({100*count/max(total,1):.1f}%)',
                    va='center', color=TEXT, fontsize=9)
        ax.set_yticks(range(len(vc)))
        ax.set_yticklabels(vc.index, color=TEXT, fontsize=9)
        ax.set_xlabel('Read count', color=TEXT)
        ax.set_title(title, color=TEXT, fontweight='bold')
        ax.set_xlim(0, total * 1.3)
        _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_coverage(coverage, counts, counts_fwd, counts_rev,
                  cov_fwd, cov_rev, cis_seq, R, itr1, itr2, name, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    pos   = np.arange(1, R + 1)

    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True)
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — Packaged Genome Coverage ({R:,} bp cis plasmid)',
                 color=TEXT, fontsize=13, fontweight='bold')

    # Mark transgene insert region
    insert_start = itr1['start_1based']
    insert_end   = itr2['end_1based']

    ax = axes[0]
    ax.fill_between(pos, coverage, color=DEV_COLOR, alpha=0.6)
    ax.axvspan(insert_start, insert_end, color='#1f6feb', alpha=0.1,
               label='Transgene insert')
    _add_itr_spans(ax, itr1, itr2)
    ax.set_ylabel('Coverage (x)', color=TEXT)
    ax.set_title('Sequencing depth', color=TEXT)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    ax  = axes[1]
    bot = np.zeros(R)
    for base, color in BCOL.items():
        ax.bar(pos, counts[base], bottom=bot, color=color, label=base, width=1, alpha=0.9)
        bot += counts[base]
    _add_itr_spans(ax, itr1, itr2)
    ax.set_ylabel('Read count', color=TEXT)
    ax.set_title('Base composition', color=TEXT)
    ax.legend(loc='lower right', ncol=4, fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    ax = axes[2]
    ref_cnt   = np.array([counts[cis_seq[i]][i] if cis_seq[i] in counts else 0
                          for i in range(R)], dtype=np.int32)
    dev_pct   = (1 - ref_cnt / cov_s) * 100
    # Mask positions with no coverage
    covered   = coverage >= 1
    flag_mask = (dev_pct >= 10) & (coverage >= 50)
    clean     = covered & ~flag_mask
    ax.bar(pos[clean],      dev_pct[clean],      width=1, color=DEV_COLOR, alpha=0.4)
    ax.bar(pos[flag_mask],  dev_pct[flag_mask],  width=1, color=FLAG_COLOR,
           alpha=0.9, label='>=10% deviation')
    _add_itr_spans(ax, itr1, itr2)
    ax.axhline(10.0, color=FLAG_COLOR, ls='--', lw=1, alpha=0.7)
    ax.set_ylabel('Non-ref bases (%)', color=TEXT)
    ax.set_xlabel('Position', color=TEXT)
    ax.set_title('Deviation from reference', color=TEXT)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_strand_coverage(coverage, cov_fwd, cov_rev, itr1, itr2, name, out_path):
    pos = np.arange(1, len(coverage) + 1)
    fig, axes = plt.subplots(2, 1, figsize=(18, 9), sharex=True)
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — Strand-Specific Coverage',
                 color=TEXT, fontsize=13, fontweight='bold')

    # Panel 1: forward and reverse coverage as lines on same axes
    ax = axes[0]
    ax.plot(pos, cov_fwd, color='#1a6fcc', lw=0.8, alpha=0.9, label='Forward')
    ax.plot(pos, cov_rev, color='#cc4422', lw=0.8, alpha=0.9, label='Reverse')
    _add_itr_spans(ax, itr1, itr2)
    ax.set_ylabel('Coverage (x)', color=TEXT)
    ax.set_title('Forward and reverse strand coverage', color=TEXT)
    ax.legend(fontsize=9, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    # Panel 2: strand balance ratio
    ax    = axes[1]
    total = np.where(coverage > 0, coverage, 1)
    fwd_frac = cov_fwd / total * 100
    rev_frac = cov_rev / total * 100
    ax.plot(pos, fwd_frac, color='#1a6fcc', lw=0.8, alpha=0.9, label='Forward %')
    ax.plot(pos, rev_frac, color='#cc4422', lw=0.8, alpha=0.9, label='Reverse %')
    ax.axhline(50, color=GRID, lw=1, ls='--', alpha=0.8, label='50%')
    _add_itr_spans(ax, itr1, itr2)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Strand %', color=TEXT)
    ax.set_xlabel('Position', color=TEXT)
    ax.set_title('Strand balance (expected ~50% each for ssAAV population)', color=TEXT)
    ax.legend(fontsize=9, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_truncations(reads_df, R, itr1, itr2, name, out_path):
    """
    Single-panel overlaid histogram of read start and end positions.
    Sharp peaks at ITR boundaries indicate full-length genomes;
    internal peaks indicate truncation hotspots.
    """
    mapped = reads_df[(reads_df['reference'] == 'transgene') &
                      (reads_df['genome_type'] != 'Unmapped')]
    if len(mapped) == 0:
        return

    bins = np.arange(0, R + 10, 10)
    fig, ax = plt.subplots(1, 1, figsize=(18, 6))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — Read Start and End Positions ({len(mapped):,} transgene reads)',
                 color=TEXT, fontsize=13, fontweight='bold')

    ax.hist(mapped['r_start'], bins=bins, color='#f78166', alpha=0.8,
            label='Read start')
    ax.hist(mapped['r_end'],   bins=bins, color='#3fb950', alpha=0.8,
            label='Read end')
    _add_itr_spans(ax, itr1, itr2)
    ax.set_ylabel('Reads mapped', color=TEXT)
    ax.set_xlabel('Position', color=TEXT)
    ax.legend(fontsize=9, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_read_lengths(reads_df, itr1, itr2, name, out_path):
    """Read length distribution for transgene-mapped reads, colored by genome type."""
    transgene = reads_df[reads_df['reference'] == 'transgene'].copy()
    if len(transgene) == 0:
        return

    insert_len = itr2['end_0based'] - itr1['start_0based']

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(PANEL)

    # Plot histogram for each genome type
    type_order = ['Full ssAAV', 'Partial ssAAV', 'Partial scAAV',
                  'Full scAAV', 'Transgene unclassified']
    bins = np.arange(0, transgene['read_len'].max() + 200, 100)

    bottom = np.zeros(len(bins) - 1)
    for gtype in type_order:
        sub = transgene[transgene['genome_type'] == gtype]['read_len']
        if len(sub) == 0:
            continue
        counts, _ = np.histogram(sub, bins=bins)
        color = GTYPE_COLORS.get(gtype, '#888888')
        ax.bar(bins[:-1], counts, width=100, bottom=bottom,
               color=color, alpha=0.9, label=f'{gtype} (n={len(sub):,})')
        bottom += counts

    # Mark key lengths
    ax.axvline(insert_len, color='black', lw=1.5, ls='--',
               label=f'Full genome ({insert_len:,} bp)')
    ax.axvline(insert_len * 2, color='#888888', lw=1, ls=':',
               label=f'scAAV ({insert_len*2:,} bp)')

    ax.set_xlabel('Read length (bp)', color=TEXT)
    ax.set_ylabel('Read count', color=TEXT)
    ax.set_title(f'{name} — Read Length Distribution ({len(transgene):,} transgene reads)',
                 color=TEXT, fontweight='bold')
    ax.tick_params(colors=TEXT, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.grid(axis='y', color=GRID, lw=0.5, alpha=0.6)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


def write_read_lengths_csv(reads_df, name, out_path):
    """Write per-read length data with genome type for downstream analysis."""
    transgene = reads_df[reads_df['reference'] == 'transgene'][
        ['read_id', 'read_len', 'genome_type', 'genome_subtype']
    ].copy()
    transgene.to_csv(out_path, index=False)
    print(f'  Saved: {out_path}')


def plot_itr_detail(coverage, cov_fwd, cov_rev, itr1, itr2, name, out_path):
    pos = np.arange(1, len(coverage) + 1)
    R   = len(coverage)
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — ITR Coverage Detail',
                 color=TEXT, fontsize=13, fontweight='bold')

    for col_idx, (itr, label, color) in enumerate([
        (itr1, 'ITR1', '#ffa657'),
        (itr2, 'ITR2', '#3fb950'),
    ]):
        s, e = itr['start_0based'], itr['end_0based']
        lo   = max(0, s - 30)
        hi   = min(R, e + 30)
        mask = (pos >= lo + 1) & (pos <= hi)
        p_sub   = pos[mask]
        cov_sub = coverage[mask]
        fwd_sub = cov_fwd[mask]
        rev_sub = cov_rev[mask]
        itr_sub = (p_sub >= s + 1) & (p_sub <= e)

        ax = axes[0][col_idx]
        ax.fill_between(p_sub, cov_sub, color=DEV_COLOR, alpha=0.4, label='Total')
        ax.fill_between(p_sub, fwd_sub, color='#3fb950', alpha=0.7, label='Forward')
        ax.fill_between(p_sub, rev_sub, color='#f78166', alpha=0.7, label='Reverse')
        ax.axvline(s + 1, color=color, lw=2, ls='--', alpha=0.8,
                   label=f'{label} boundary')
        ax.axvline(e, color=color, lw=2, ls='--', alpha=0.8)
        mean_c = cov_sub[itr_sub].mean() if itr_sub.sum() > 0 else 0
        min_c  = cov_sub[itr_sub].min()  if itr_sub.sum() > 0 else 0
        fwd_m  = fwd_sub[itr_sub].mean() if itr_sub.sum() > 0 else 0
        rev_m  = rev_sub[itr_sub].mean() if itr_sub.sum() > 0 else 0
        ratio  = fwd_m / max(rev_m, 0.1)
        info   = f'Mean: {mean_c:.0f}x\nMin: {min_c}x\nFwd/Rev: {ratio:.1f}\nInner loop: {itr["inner_len"]} bp'
        ax.text(0.02, 0.96, info, transform=ax.transAxes, color=TEXT, fontsize=8,
                va='top', bbox=dict(boxstyle='round', facecolor=PANEL, alpha=0.7))
        ax.set_title(f'{label} — Coverage by Strand', color=TEXT, fontweight='bold')
        ax.set_ylabel('Coverage (x)', color=TEXT)
        ax.legend(fontsize=7, facecolor=PANEL, labelcolor=TEXT, ncol=4)
        _style(ax)

        ax    = axes[1][col_idx]
        total = np.where(cov_sub > 0, cov_sub, 1)
        ax.fill_between(p_sub, fwd_sub / total,
                        color='#3fb950', alpha=0.8, label='Forward fraction')
        ax.fill_between(p_sub, -(rev_sub / total),
                        color='#f78166', alpha=0.8, label='Reverse fraction')
        ax.axhline(0,    color=TEXT, lw=0.5, alpha=0.5)
        ax.axhline(0.5,  color=GRID, lw=1, ls='--', alpha=0.7)
        ax.axhline(-0.5, color=GRID, lw=1, ls='--', alpha=0.7)
        ax.axvline(s + 1, color=color, lw=2, ls='--', alpha=0.8)
        ax.axvline(e,     color=color, lw=2, ls='--', alpha=0.8)
        ax.set_ylim(-1.1, 1.1)
        ax.set_title(f'{label} — Strand Balance', color=TEXT, fontweight='bold')
        ax.set_ylabel('Strand fraction', color=TEXT)
        ax.set_xlabel('Position', color=TEXT)
        ax.legend(fontsize=7, facecolor=PANEL, labelcolor=TEXT)
        _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {out_path}')


# ── CSV outputs ────────────────────────────────────────────────────────────────
def write_contamination_csv(contam_counts, total, name, out_path):
    nice = {'transgene': 'Transgene (cis plasmid)',
            'helper': 'Helper plasmid', 'rep_cap': 'Rep-Cap plasmid',
            'host': 'Host (E. coli)', 'unmapped': 'Unmapped',
            'reverse_packaging': 'Reverse packaging (backbone)'}
    rows = []
    for ref in ['transgene', 'helper', 'rep_cap', 'host',
                'reverse_packaging', 'unmapped']:
        count = contam_counts.get(ref, 0)
        rows.append({
            'sample':      name,
            'reference':   nice.get(ref, ref),
            'count':       count,
            'pct_of_total': round(100 * count / max(total, 1), 2),
        })
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f'  Saved: {out_path}')


def write_genome_types_csv(reads_df, name, out_path):
    transgene_df = reads_df[reads_df['reference'] == 'transgene']
    total = len(transgene_df)
    rows  = []
    for gtype in transgene_df['genome_type'].unique():
        sub = transgene_df[transgene_df['genome_type'] == gtype]
        for gsubtype in sub['genome_subtype'].unique():
            count = (sub['genome_subtype'] == gsubtype).sum()
            rows.append({
                'sample':         name,
                'genome_type':    gtype,
                'genome_subtype': gsubtype,
                'count':          count,
                'pct_of_transgene': round(100 * count / max(total, 1), 2),
                'pct_of_mapped':  round(100 * count / max(
                    (transgene_df['genome_type'] != 'Unmapped').sum(), 1), 2),
            })
    pd.DataFrame(rows).sort_values(['genome_type', 'count'],
                                    ascending=[True, False]).to_csv(out_path, index=False)
    print(f'  Saved: {out_path}')


def write_pileup_csv(coverage, counts, counts_fwd, counts_rev,
                     cov_fwd, cov_rev, cis_seq, R, itr1, itr2, out_path):
    RC_BASE = {'A':'T','T':'A','G':'C','C':'G'}
    cov_s      = np.where(coverage > 0, coverage, 1)
    cov_fwd_s  = np.where(cov_fwd > 0, cov_fwd, 1)
    cov_rev_s  = np.where(cov_rev > 0, cov_rev, 1)

    # Standard deviation (total non-ref / total)
    ref_cnt   = np.array([counts[cis_seq[i]][i] if cis_seq[i] in counts else 0
                          for i in range(R)], dtype=np.int32)
    dev_pct   = (1 - ref_cnt / cov_s) * 100

    # Strand-aware deviation:
    # fwd reads vs reference; rev reads vs RC(reference)
    fwd_ref_cnt = np.array([counts_fwd[cis_seq[i]][i] if cis_seq[i] in counts_fwd else 0
                            for i in range(R)], dtype=np.int32)
    rev_ref_cnt = np.array([counts_rev[RC_BASE.get(cis_seq[i], cis_seq[i])][i]
                            if RC_BASE.get(cis_seq[i], cis_seq[i]) in counts_rev else 0
                            for i in range(R)], dtype=np.int32)
    dev_fwd   = (1 - fwd_ref_cnt / cov_fwd_s) * 100
    dev_rev   = (1 - rev_ref_cnt / cov_rev_s) * 100
    # Combined strand-aware: non-ref on fwd + non-rc-ref on rev / total
    dev_strand_aware = ((cov_fwd - fwd_ref_cnt) + (cov_rev - rev_ref_cnt)) / cov_s * 100

    flag_mask        = (dev_pct >= 10) & (coverage >= 50)
    flag_strand_aware = (dev_strand_aware >= 10) & (coverage >= 50)

    s1, e1 = itr1['start_0based'], itr1['end_0based']
    s2, e2 = itr2['start_0based'], itr2['end_0based']

    def region(p):
        if s1 + 1 <= p <= e1: return 'ITR1'
        if s2 + 1 <= p <= e2: return 'ITR2'
        if e1 < p <= s2:      return 'Transgene_insert'
        return 'Backbone'

    pos_1 = np.arange(1, R + 1)
    df_full = pd.DataFrame({
        'Position':             pos_1,
        'Region':               [region(p) for p in pos_1],
        'Reference':            list(cis_seq),
        'Coverage':             coverage,
        'Coverage_fwd':         cov_fwd,
        'Coverage_rev':         cov_rev,
        'A':                    counts['A'],
        'T':                    counts['T'],
        'G':                    counts['G'],
        'C':                    counts['C'],
        'A_fwd':                counts_fwd['A'],
        'T_fwd':                counts_fwd['T'],
        'G_fwd':                counts_fwd['G'],
        'C_fwd':                counts_fwd['C'],
        'A_rev':                counts_rev['A'],
        'T_rev':                counts_rev['T'],
        'G_rev':                counts_rev['G'],
        'C_rev':                counts_rev['C'],
        'Deviation_pct':        np.round(dev_pct, 3),
        'Deviation_fwd':        np.round(dev_fwd, 3),
        'Deviation_rev':        np.round(dev_rev, 3),
        'Deviation_strand_aware': np.round(dev_strand_aware, 3),
        'Flag_10pct':           flag_mask,
        'Flag_strand_aware':    flag_strand_aware,
    })
    # Restrict output to ITR-to-ITR insert region only
    s1 = itr1['start_0based']
    e2 = itr2['end_0based']
    df_insert = df_full[(df_full['Position'] >= s1 + 1) &
                        (df_full['Position'] <= e2)].copy()
    df_insert.to_csv(out_path, index=False)
    print(f'  Saved: {out_path}')


def write_summary_csv(reads_df, coverage, cov_fwd, cov_rev,
                      contam_counts, itr_cov_df, name, total, out_path):
    tdf     = reads_df[reads_df['reference'] == 'transgene']
    mapped  = (tdf['genome_type'] != 'Unmapped').sum()
    full    = (tdf['genome_type'] == 'Full ssAAV').sum()
    partial = (tdf['genome_type'] == 'Partial ssAAV').sum()
    sc      = tdf['genome_type'].isin(['Full scAAV', 'Partial scAAV']).sum()
    fwd_t   = int(cov_fwd.sum())
    rev_t   = int(cov_rev.sum())

    row = {
        'sample':                name,
        'total_reads':           total,
        'pct_transgene':         round(100*contam_counts.get('transgene',0)/max(total,1),2),
        'pct_helper':            round(100*contam_counts.get('helper',0)/max(total,1),2),
        'pct_rep_cap':           round(100*contam_counts.get('rep_cap',0)/max(total,1),2),
        'pct_host':              round(100*contam_counts.get('host',0)/max(total,1),2),
        'pct_unmapped':          round(100*contam_counts.get('unmapped',0)/max(total,1),2),
        'pct_reverse_packaging': round(100*contam_counts.get('reverse_packaging',0)/max(total,1),2),
        'transgene_reads':       contam_counts.get('transgene', 0),
        'pct_full_ssAAV':        round(100*full/max(mapped,1),2),
        'pct_partial_ssAAV':     round(100*partial/max(mapped,1),2),
        'pct_scAAV':             round(100*sc/max(mapped,1),2),
        'mean_coverage':         round(float(coverage.mean()),1),
        'min_coverage':          int(coverage.min()),
        'pct_fwd_reads':         round(100*fwd_t/max(fwd_t+rev_t,1),1),
        'pct_rev_reads':         round(100*rev_t/max(fwd_t+rev_t,1),1),
    }
    for _, r in itr_cov_df.iterrows():
        k = r['itr'].lower()
        row[f'{k}_mean_cov']      = r['mean_coverage']
        row[f'{k}_min_cov']       = r['min_coverage']
        row[f'{k}_fwd_rev_ratio'] = r['strand_ratio_fwd_rev']

    pd.DataFrame([row]).to_csv(out_path, index=False)
    print(f'  Saved: {out_path}')


def validate_against_tsv(reads_df, tsv_path):
    ref_df = pd.read_csv(tsv_path, sep='\t')
    ref_df.columns = [c.lower().replace(' ', '_') for c in ref_df.columns]
    ref_df = ref_df.rename(columns={'read': 'read_id'})
    merged = reads_df.merge(ref_df, on='read_id', how='inner')
    if len(merged) == 0:
        print('  WARNING: no matching read IDs')
        return
    total = len(merged)
    match = (merged['genome_type'] == merged['assigned_genome_type']).sum()
    print(f'\n  Validation vs wf-aav-qc: {match:,}/{total:,} ({100*match/total:.1f}%) agreement')
    print(f'  Discrepancies:')
    for (ours, theirs), count in \
            merged.groupby(['genome_type', 'assigned_genome_type']).size().items():
        if ours != theirs:
            print(f'    {ours} -> {theirs}: {count}')


# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='AAV genome QC pipeline v3',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--fastq-dir',        required=True,
                        help='Directory containing fastq.gz files')
    parser.add_argument('--cis-plasmid',      required=True,
                        help='Full cis plasmid FASTA (ITRs auto-detected)')
    parser.add_argument('--rep-cap',          required=True,
                        help='Rep-Cap plasmid FASTA')
    parser.add_argument('--helper',           default=str(DEFAULT_HELPER),
                        help=f'Helper plasmid FASTA (default: bundled pHelper)')
    parser.add_argument('--host',             default=str(DEFAULT_HOST),
                        help=f'Host genome FASTA (default: bundled E. coli K-12)')
    parser.add_argument('--serotype',         default='AAV2',
                        help='AAV serotype for ITR detection (default: AAV2)')
    parser.add_argument('--name',             default='sample')
    parser.add_argument('--output-dir',       default='.')
    parser.add_argument('--itr-tol',          type=int, default=ITR_TOL)
    parser.add_argument('--min-mapq',         type=int, default=MIN_MAPQ)
    parser.add_argument('--min-base-qual',    type=int, default=20,
                        help='Minimum base quality score (Phred, default 10)')
    parser.add_argument('--no-contamination', action='store_true',
                        help='Skip contamination analysis (faster)')
    parser.add_argument('--validate-tsv',     default=None)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'\n[1/6] Loading cis plasmid and detecting ITRs...')
    cis_seq = load_fasta(args.cis_plasmid)
    R       = len(cis_seq)
    print(f'  Cis plasmid: {R:,} bp')

    itr_regions = auto_detect_itrs(cis_seq, serotype=args.serotype, verbose=True)
    itr1, itr2  = itr_regions[0], itr_regions[1]

    print(f'\n  Auto-detected ITRs:')
    for i, r in enumerate([itr1, itr2]):
        print(f'    ITR{i+1}: pos {r["start_1based"]}-{r["end_1based"]} '
              f'({r["length"]} bp, arm={r["arm_len"]}bp)')
        print(f'           Inner loop ({r["inner_len"]} bp): {r["inner_loop_flip"]}')
        check_deleted_itr(r, i+1)

    # Write auto-detected BED file
    plasmid_name = Path(args.cis_plasmid).stem
    bed_path = out_dir / f'{args.name}_detected_itrs.bed'
    write_bed(itr_regions[:2], plasmid_name, bed_path)

    print(f'\n[2/6] Building alignment indices...')
    if args.no_contamination:
        print('  Contamination analysis skipped (--no-contamination)')
        aligners = {
            'transgene': mappy.Aligner(seq=cis_seq, preset='map-ont', best_n=1)
        }
    else:
        aligners = build_contamination_index(
            cis_seq,
            helper_path  = args.helper,
            host_path    = args.host,
            rep_cap_path = args.rep_cap,
            min_mapq     = args.min_mapq,
        )

    print(f'\n[3/6] Aligning reads...')
    (coverage, counts, counts_fwd, counts_rev, cov_fwd, cov_rev,
     reads_df, contam_counts, total) = run_analysis(
        args.fastq_dir, cis_seq, R, itr1, itr2, aligners,
        tol=args.itr_tol, min_mapq=args.min_mapq,
        min_base_qual=args.min_base_qual
    )

    if args.validate_tsv:
        validate_against_tsv(reads_df, args.validate_tsv)

    print(f'\n[4/6] Analyzing ITR coverage...')
    itr_cov_df = analyze_itr_coverage(
        coverage, counts_fwd, counts_rev, cov_fwd, cov_rev,
        cis_seq, itr1, itr2, R)

    print(f'\n[5/6] Writing CSV outputs...')
    write_contamination_csv(contam_counts, total, args.name,
                            out_dir / f'{args.name}_contamination.csv')
    write_genome_types_csv(reads_df, args.name,
                           out_dir / f'{args.name}_genome_types.csv')
    write_pileup_csv(coverage, counts, counts_fwd, counts_rev,
                     cov_fwd, cov_rev, cis_seq, R, itr1, itr2,
                     out_dir / f'{args.name}_pileup.csv')
    itr_cov_df.to_csv(out_dir / f'{args.name}_itr_coverage.csv', index=False)
    print(f'  Saved: {out_dir}/{args.name}_itr_coverage.csv')
    reads_df[reads_df['reference'] == 'transgene'][
        ['read_id','read_len','genome_type','genome_subtype',
         'r_start','r_end','mapq','strand']
    ].to_csv(out_dir / f'{args.name}_truncations.csv', index=False)
    print(f'  Saved: {out_dir}/{args.name}_truncations.csv')
    write_read_lengths_csv(reads_df, args.name,
                           out_dir / f'{args.name}_read_lengths.csv')
    write_summary_csv(reads_df, coverage, cov_fwd, cov_rev,
                      contam_counts, itr_cov_df, args.name, total,
                      out_dir / f'{args.name}_summary.csv')

    print(f'\n[6/6] Generating figures...')
    if not args.no_contamination:
        plot_contamination(contam_counts, total, args.name,
                           out_dir / f'{args.name}_contamination.png')
    plot_read_lengths(reads_df, itr1, itr2, args.name,
                      out_dir / f'{args.name}_read_lengths.png')
    plot_structures(reads_df, args.name,
                    out_dir / f'{args.name}_structures.png')
    plot_coverage(coverage, counts, counts_fwd, counts_rev,
                  cov_fwd, cov_rev, cis_seq, R, itr1, itr2, args.name,
                  out_dir / f'{args.name}_coverage.png')
    plot_strand_coverage(coverage, cov_fwd, cov_rev, itr1, itr2, args.name,
                         out_dir / f'{args.name}_strand_coverage.png')
    plot_truncations(reads_df, R, itr1, itr2, args.name,
                     out_dir / f'{args.name}_truncations.png')
    plot_itr_detail(coverage, cov_fwd, cov_rev, itr1, itr2, args.name,
                    out_dir / f'{args.name}_itr_detail.png')

    print(f'\nDone.\nResults in: {out_dir}\n')


if __name__ == '__main__':
    main()
