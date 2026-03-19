#!/usr/bin/env python3
"""
AAV Plasmid ONT Pileup Analysis Workflow
=========================================
Standard pipeline for per-position basecall heterogeneity analysis
of AAV plasmid vectors from Oxford Nanopore raw FASTQ reads.

Usage:
    python aav_plasmid_workflow.py --fastq reads.fastq.gz --fasta reference.fasta \
        --right-itr 1838 1982 --left-itr 84 228 \
        --name pAAVcis --output-dir ./results/

    python aav_plasmid_workflow.py --fastq reads.fastq.gz --fasta reference.fasta \
        --right-itr 4343 4479 --left-itr 8869 9012 \
        --name pSM630 --output-dir ./results/

Arguments:
    --fastq         Raw ONT reads (gzip or plain FASTQ)
    --fasta         Reference plasmid sequence (single contig)
    --right-itr     1-based inclusive coords of right ITR (start end)
    --left-itr      1-based inclusive coords of left ITR (start end)
    --name          Plasmid name (used in output filenames and titles)
    --output-dir    Directory for outputs (created if absent)
    --flag-pct      Deviation threshold for flagging (default 10.0)
    --min-cov       Minimum coverage to include a position (default 50)
    --no-flipflop   Skip flip/flop inner loop analysis

Outputs (in output-dir/):
    {name}_pileup.csv       Per-position pileup table
    {name}_pileup.png       Full-plasmid 3-panel figure
    {name}_itr_detail.png   ITR detail 2×2 figure
    {name}_flipflop.txt     Flip/flop quantification report
    {name}_summary.txt      Run summary statistics

Algorithm:
    Circular-aware seed-chain alignment with banded Needleman-Wunsch.
    Handles reads spanning the linearization point (wrap-around reads).
    All query bases are recorded at their reference-aligned positions.
    No external aligner required (minimap2 / samtools not needed).

Notes on coordinate offset:
    Azenta's Plasmid-EZ pipeline may linearize the plasmid at a different
    position than the uploaded FASTA. Before running, verify that the first
    ~40 bp of the Azenta reference column matches your FASTA from position 1.
    If there is an offset N, rotate the reference:
        ref = ref_orig[N:] + ref_orig[:N]
    and adjust ITR coordinates accordingly (subtract N, modulo plasmid length).
    The script reports the first 40 bp of the reference so you can verify.
"""

import argparse
import collections
import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────
BASES = 'ATGC'
K = 12          # k-mer length for seeding
STRIDE = 3      # sampling stride for seed placement
MAX_KMER_OCC = 3   # ignore k-mers occurring more than this many times (repeats)
BAND = 12       # Needleman-Wunsch band half-width for inter-anchor gaps
END_CLIP = 80   # max bp to align at read ends
WRAP_TOL = 100  # tolerance for wrap-around diagonal detection

# ── Utilities ──────────────────────────────────────────────────────────────────
def rc(s):
    return s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def load_reference(fasta_path):
    lines = Path(fasta_path).read_text().strip().splitlines()
    seq = ''.join(l for l in lines if not l.startswith('>')).upper()
    return seq


def load_reads(fastq_path):
    opener = gzip.open if str(fastq_path).endswith('.gz') else open
    reads = []
    with opener(fastq_path, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip().upper()
            f.readline()
            f.readline()
            if seq:
                reads.append(seq)
    return reads


def build_kmer_index(ref, k=K, max_occ=MAX_KMER_OCC):
    idx = collections.defaultdict(list)
    for i in range(len(ref) - k + 1):
        idx[ref[i:i+k]].append(i)
    # Remove high-frequency k-mers (repeats / low complexity)
    idx = {kmer: pos for kmer, pos in idx.items() if len(pos) <= max_occ}
    return idx


# ── Banded Needleman-Wunsch ────────────────────────────────────────────────────
def banded_nw_pairs(query, ref_seg, band=BAND):
    """Return list of (query_idx, ref_idx) aligned pairs. ref_idx=-1 for insertions."""
    nq, nr = len(query), len(ref_seg)
    if nq == 0 or nr == 0:
        return []
    NEG = -10**6
    dp = {(0, 0): 0}
    tb = {}
    for i in range(1, nq+1):
        dp[(i, 0)] = -2*i
        tb[(i, 0)] = 'U'
    for j in range(1, nr+1):
        dp[(0, j)] = -2*j
        tb[(0, j)] = 'L'
    ratio = nr / max(nq, 1)
    for i in range(1, nq+1):
        jc = round(i * ratio)
        for j in range(max(1, jc-band), min(nr, jc+band)+1):
            d = dp.get((i-1, j-1), NEG) + (1 if query[i-1] == ref_seg[j-1] else -1)
            u = dp.get((i-1, j), NEG) - 2
            l = dp.get((i, j-1), NEG) - 2
            best = max(d, u, l)
            dp[(i, j)] = best
            tb[(i, j)] = 'D' if best == d else ('U' if best == u else 'L')
    pairs = []
    i, j = nq, nr
    while i > 0 or j > 0:
        mv = tb.get((i, j), 'D')
        if mv == 'D':
            i -= 1; j -= 1
            pairs.append((i, j))
        elif mv == 'U':
            i -= 1
            pairs.append((i, -1))
        else:
            j -= 1
    pairs.reverse()
    return pairs


def fill_gap(pileup, s, qi0, ri0, qi1, ri1, R):
    if qi1 <= qi0 or ri1 <= ri0:
        return
    q_seg = s[qi0:qi1]
    r_seg_len = ri1 - ri0
    # Build a flat reference segment (no wrapping needed for short gaps)
    for dqi, dri in banded_nw_pairs(q_seg, 'N' * r_seg_len):
        if dri < 0:
            continue
        rp = ri0 + dri
        qp = qi0 + dqi
        if 0 <= rp < R and qp < len(s) and rp not in pileup:
            pileup[rp] = s[qp]


# ── Read alignment ─────────────────────────────────────────────────────────────
def align_read(seq, ref, idx, R, k=K, band=30, end_clip=END_CLIP):
    """
    Align a single read to the circular reference.
    Returns dict {ref_position: base} or None if alignment fails.
    """
    best_strand, best_cnt, best_hits = '+', 0, []
    for strand, s in (('+', seq), ('-', rc(seq))):
        hits = [
            (qi, ri)
            for qi in range(0, len(s)-k+1, STRIDE)
            for ri in idx.get(s[qi:qi+k], [])
        ]
        if len(hits) > best_cnt:
            best_cnt = len(hits)
            best_strand = strand
            best_hits = hits

    if best_cnt < 5:
        return None

    s = seq if best_strand == '+' else rc(seq)

    # Find dominant diagonal(s); detect wrap-around reads
    dc = collections.Counter(ri - qi for qi, ri in best_hits)
    top = [d for d, _ in dc.most_common(6)]
    main_diag = top[0]
    wrap_diag = None
    for d in top[1:]:
        if abs(abs(main_diag - d) - R) < WRAP_TOL:
            wrap_diag = d
            break

    if wrap_diag is not None:
        d_lo, d_hi = sorted([main_diag, wrap_diag])
        hits_lo = [(qi, ri) for qi, ri in best_hits if abs(ri-qi-d_lo) <= band]
        hits_hi = [(qi, ri) for qi, ri in best_hits if abs(ri-qi-d_hi) <= band]
        if hits_lo and hits_hi:
            trans = (max(qi for qi, ri in hits_hi) + min(qi for qi, ri in hits_lo)) // 2
            segments = [(0, trans, d_hi), (trans, len(s), d_lo)]
        else:
            segments = [(0, len(s), main_diag)]
    else:
        segments = [(0, len(s), main_diag)]

    pileup = {}
    for seg_start, seg_end, diag in segments:
        seg_hits = sorted(
            [(qi, ri) for qi, ri in best_hits
             if seg_start <= qi < seg_end and abs(ri-qi-diag) <= band],
            key=lambda x: x[0]
        )
        if not seg_hits:
            continue

        # Record k-mer anchor bases
        for qi, ri in seg_hits:
            for kp in range(k):
                rp = ri + kp; qp = qi + kp
                if 0 <= rp < R and qp < len(s) and rp not in pileup:
                    pileup[rp] = s[qp]

        # Fill inter-anchor gaps with banded NW
        anchors = [(qi, ri, qi+k, ri+k) for qi, ri in seg_hits]
        for i in range(len(anchors)-1):
            _, _, qi0, ri0 = anchors[i]
            qi1, ri1, _, _ = anchors[i+1]
            fill_gap(pileup, s, qi0, ri0, qi1, ri1, R)

        # Align read ends
        fqi, fri = seg_hits[0]
        cq = min(fqi, end_clip); cr = min(fri, end_clip)
        if cq > 0 and cr > 0:
            fill_gap(pileup, s, fqi-cq, fri-cr, fqi, fri, R)

        lqi, lri = seg_hits[-1]; lqi += k; lri += k
        cq = min(len(s)-lqi, end_clip); cr = min(R-lri, end_clip)
        if cq > 0 and cr > 0:
            fill_gap(pileup, s, lqi, lri, lqi+cq, lri+cr, R)

    return pileup if pileup else None


# ── Pileup accumulation ────────────────────────────────────────────────────────
def run_pileup(reads, ref, idx, R, verbose=True):
    counts = {b: np.zeros(R, dtype=np.int32) for b in BASES}
    coverage = np.zeros(R, dtype=np.int32)
    aligned = 0
    t0 = time.time()

    for i, seq in enumerate(reads):
        if verbose and i % 500 == 0:
            print(f'  Aligning read {i}/{len(reads)} ({time.time()-t0:.1f}s)', flush=True)
        proj = align_read(seq, ref, idx, R)
        if proj:
            aligned += 1
            for rpos, base in proj.items():
                if 0 <= rpos < R:
                    coverage[rpos] += 1
                    if base in counts:
                        counts[base][rpos] += 1

    if verbose:
        print(f'  Aligned: {aligned}/{len(reads)} ({100*aligned/len(reads):.1f}%) '
              f'in {time.time()-t0:.1f}s')
        print(f'  Mean coverage: {coverage.mean():.0f}x, min: {coverage.min()}x')

    return coverage, counts, aligned


# ── Flip/flop analysis ─────────────────────────────────────────────────────────
def find_palindromic_arms(itr_seq, min_arm=30, max_arm=44):
    """Find the longest palindromic arm pair within an ITR sequence."""
    best_len, best_i, best_j = 0, -1, -1
    for armlen in range(max_arm, min_arm-1, -1):
        for i in range(len(itr_seq)-armlen):
            seg = itr_seq[i:i+armlen]
            seg_rc = rc(seg)
            for j in range(i+armlen+1, len(itr_seq)-armlen+1):
                if itr_seq[j:j+armlen] == seg_rc:
                    if armlen > best_len:
                        best_len = armlen
                        best_i = i; best_j = j
        if best_len >= armlen:
            break
    return best_i, best_len, best_j  # arm_a_start, arm_len, arm_b_start


def analyze_flipflop(coverage, counts, ref, itr_lo, itr_hi, itr_name, cov_s):
    """
    Quantify flip/flop signal in an ITR inner loop.
    itr_lo, itr_hi: 0-based half-open indices [lo, hi)
    Returns a dict of results.
    """
    itr_seq = ref[itr_lo:itr_hi]
    arm_a_start, arm_len, arm_b_start = find_palindromic_arms(itr_seq)

    if arm_a_start < 0:
        return {'found_arms': False, 'itr_name': itr_name}

    inner_lo_rel = arm_a_start + arm_len
    inner_hi_rel = arm_b_start
    inner_seq = itr_seq[inner_lo_rel:inner_hi_rel]
    flop_seq = rc(inner_seq)

    inner_lo_abs = itr_lo + inner_lo_rel
    inner_hi_abs = itr_lo + inner_hi_rel

    flop_counts = []
    total_covs = []
    informative = 0

    results_per_pos = []
    for rel_i, (flip_b, flop_b) in enumerate(zip(inner_seq, flop_seq)):
        if flip_b == flop_b:
            continue  # cannot distinguish at this position
        abs_i = inner_lo_abs + rel_i
        a = counts['A'][abs_i]; t = counts['T'][abs_i]
        g = counts['G'][abs_i]; c = counts['C'][abs_i]
        cov = coverage[abs_i]
        flop_n = {'A': a, 'T': t, 'G': g, 'C': c}[flop_b]
        flip_n = {'A': a, 'T': t, 'G': g, 'C': c}[flip_b]
        flop_pct = 100 * flop_n / cov if cov > 0 else 0.0
        flop_counts.append(flop_n)
        total_covs.append(cov)
        informative += 1
        results_per_pos.append({
            'pos_1based': abs_i + 1,
            'inner_idx': rel_i,
            'flip_base': flip_b,
            'flop_base': flop_b,
            'flip_count': flip_n,
            'flop_count': flop_n,
            'coverage': cov,
            'flop_pct': round(flop_pct, 2),
        })

    mean_flop = 100 * sum(flop_counts) / max(sum(total_covs), 1)

    return {
        'found_arms': True,
        'itr_name': itr_name,
        'arm_len': arm_len,
        'arm_a_abs': (itr_lo + arm_a_start + 1, itr_lo + arm_a_start + arm_len),
        'arm_b_abs': (itr_lo + arm_b_start + 1, itr_lo + arm_b_start + arm_len),
        'inner_abs': (inner_lo_abs + 1, inner_hi_abs),
        'inner_seq_flip': inner_seq,
        'inner_seq_flop': flop_seq,
        'n_informative': informative,
        'mean_flop_pct': round(mean_flop, 2),
        'per_position': results_per_pos,
    }


# ── Output: CSV ────────────────────────────────────────────────────────────────
def write_csv(coverage, counts, ref, R, right_itr, left_itr, flag_pct, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    ref_cnt = np.array([counts[ref[i]][i] for i in range(R)], dtype=np.int32)
    dev_pct = (1 - ref_cnt / cov_s) * 100
    pos_1 = np.arange(1, R+1)

    def region(pos):
        if right_itr and right_itr[0] <= pos <= right_itr[1]:
            return 'Right_ITR'
        if left_itr and left_itr[0] <= pos <= left_itr[1]:
            return 'Left_ITR'
        return 'Backbone'

    flag_mask = (dev_pct >= flag_pct) & (coverage >= 50)

    df = pd.DataFrame({
        'Position':      pos_1,
        'Region':        [region(p) for p in pos_1],
        'Reference':     list(ref),
        'Coverage':      coverage,
        'A':             counts['A'],
        'T':             counts['T'],
        'G':             counts['G'],
        'C':             counts['C'],
        'Deviation_pct': np.round(dev_pct, 3),
        'Flag_10pct':    flag_mask,
    })
    df.to_csv(out_path, index=False)
    return df, dev_pct, flag_mask


# ── Output: Figures ────────────────────────────────────────────────────────────
BG = '#0d1117'; PANEL = '#161b22'; GRID = '#21262d'; TEXT = 'white'
BCOL = {'A': '#3fb950', 'T': '#f78166', 'G': '#d2a8ff', 'C': '#ffa657'}
DEV_COLOR = '#58a6ff'; FLAG_COLOR = '#f85149'
RITR_COLOR = '#ffa657'; LITR_COLOR = '#3fb950'


def _style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TEXT, labelsize=8)
    ax.yaxis.label.set_color(TEXT)
    ax.xaxis.label.set_color(TEXT)
    ax.title.set_color(TEXT)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.grid(axis='y', color=GRID, lw=0.5, alpha=0.6)


def plot_full_plasmid(coverage, counts, ref, R, dev_pct, flag_mask,
                      right_itr, left_itr, name, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    pos = np.arange(1, R+1)

    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True)
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} ONT Pileup — Full Plasmid ({R:,} bp)',
                 color=TEXT, fontsize=13, fontweight='bold')

    # Panel 1: coverage
    ax = axes[0]
    ax.fill_between(pos, coverage, color=DEV_COLOR, alpha=0.6)
    if right_itr:
        ax.axvspan(*right_itr, color=RITR_COLOR, alpha=0.25, label='Right ITR')
    if left_itr:
        ax.axvspan(*left_itr, color=LITR_COLOR, alpha=0.25, label='Left ITR')
    ax.set_ylabel('Coverage (×)', color=TEXT)
    ax.set_title('Coverage', color=TEXT)
    if right_itr or left_itr:
        ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    # Panel 2: base composition
    ax = axes[1]
    bot = np.zeros(R)
    for base, color in BCOL.items():
        frac = counts[base] / cov_s
        ax.bar(pos, frac, bottom=bot, color=color, label=base, width=1, alpha=0.9)
        bot += frac
    if right_itr:
        ax.axvspan(*right_itr, color=RITR_COLOR, alpha=0.2)
    if left_itr:
        ax.axvspan(*left_itr, color=LITR_COLOR, alpha=0.2)
    ax.set_ylabel('Base fraction', color=TEXT)
    ax.set_title('Base composition', color=TEXT)
    ax.legend(loc='lower right', ncol=4, fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    # Panel 3: deviation
    ax = axes[2]
    ritr_mask = np.zeros(R, bool)
    litr_mask = np.zeros(R, bool)
    if right_itr:
        ritr_mask[right_itr[0]-1:right_itr[1]] = True
    if left_itr:
        litr_mask[left_itr[0]-1:left_itr[1]] = True

    clean = ~flag_mask
    ax.bar(pos[clean], dev_pct[clean], width=1, color=DEV_COLOR, alpha=0.4)
    if right_itr:
        m = flag_mask & ritr_mask
        ax.bar(pos[m], dev_pct[m], width=1, color=RITR_COLOR, alpha=0.9, label='Right ITR')
    if left_itr:
        m = flag_mask & litr_mask
        ax.bar(pos[m], dev_pct[m], width=1, color=LITR_COLOR, alpha=0.9, label='Left ITR')
    other = flag_mask & ~ritr_mask & ~litr_mask
    ax.bar(pos[other], dev_pct[other], width=1, color=FLAG_COLOR, alpha=0.9, label='Other')
    ax.axhline(10.0, color=FLAG_COLOR, ls='--', lw=1, alpha=0.7, label='10% threshold')
    ax.set_ylabel('Non-ref bases (%)', color=TEXT)
    ax.set_xlabel('Position', color=TEXT)
    ax.set_title('Deviation from reference', color=TEXT)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT, ncol=4)
    _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()


def plot_itr_detail(coverage, counts, ref, R, dev_pct, flag_mask,
                    right_itr, left_itr, flipflop_results, name, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    pos = np.arange(1, R+1)
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — ITR Detail', color=TEXT, fontsize=13, fontweight='bold')

    itrs = []
    if right_itr:
        itrs.append((right_itr, 'Right ITR', RITR_COLOR))
    if left_itr:
        itrs.append((left_itr, 'Left ITR', LITR_COLOR))

    for col_idx, (itr, label, color) in enumerate(itrs[:2]):
        lo = max(1, itr[0] - 30)
        hi = min(R, itr[1] + 10)
        mask = (pos >= lo) & (pos <= hi)
        p_sub = pos[mask]
        dev_sub = dev_pct[mask]
        cov_sub = cov_s[mask]
        flag_sub = flag_mask[mask]

        # Top row: base composition
        ax = axes[0][col_idx]
        bot = np.zeros(mask.sum())
        for base, bc in BCOL.items():
            frac = counts[base][mask] / cov_sub
            ax.bar(p_sub, frac, bottom=bot, color=bc, label=base, width=1, alpha=0.9)
            bot += frac
        ax.axvline(itr[0], color=color, lw=2, ls='--', alpha=0.8, label='ITR boundary')
        ax.axvline(itr[1], color=color, lw=2, ls='--', alpha=0.8)

        # Mark inner loop if flip/flop results available
        ff = next((r for r in flipflop_results if r.get('found_arms') and label in r['itr_name']), None)
        if ff:
            ax.axvspan(ff['inner_abs'][0], ff['inner_abs'][1],
                       color='#ffffff', alpha=0.1, label=f'Inner loop ({ff["mean_flop_pct"]:.1f}% flop)')
        ax.set_title(f'{label} — Base Composition', color=TEXT, fontweight='bold')
        ax.set_ylabel('Fraction', color=TEXT)
        ax.legend(ncol=5, fontsize=7, facecolor=PANEL, labelcolor=TEXT)
        _style(ax)

        # Bottom row: deviation
        ax = axes[1][col_idx]
        ax.bar(p_sub[~flag_sub], dev_sub[~flag_sub], width=1, color=DEV_COLOR, alpha=0.6)
        ax.bar(p_sub[flag_sub], dev_sub[flag_sub], width=1, color=color, alpha=0.9)
        ax.axhline(10.0, color=FLAG_COLOR, ls='--', lw=1, alpha=0.7)
        ax.axvline(itr[0], color=color, lw=2, ls='--', alpha=0.8)
        ax.axvline(itr[1], color=color, lw=2, ls='--', alpha=0.8)
        if ff:
            ax.axvspan(ff['inner_abs'][0], ff['inner_abs'][1], color='#ffffff', alpha=0.1)

        itr_mask_sub = (p_sub >= itr[0]) & (p_sub <= itr[1])
        n_flag = (flag_sub & itr_mask_sub).sum()
        mean_dev = dev_sub[itr_mask_sub].mean() if itr_mask_sub.sum() > 0 else 0
        info = f'Mean: {mean_dev:.2f}%\nFlagged: {n_flag}/{itr_mask_sub.sum()}'
        if ff:
            info += f'\nFlop: {ff["mean_flop_pct"]:.1f}%'
        ax.text(0.02, 0.96, info, transform=ax.transAxes, color=TEXT, fontsize=9,
                va='top', bbox=dict(boxstyle='round', facecolor=PANEL, alpha=0.7))
        ax.set_title(f'{label} — Deviation', color=TEXT, fontweight='bold')
        ax.set_ylabel('Non-ref (%)', color=TEXT)
        ax.set_xlabel('Position', color=TEXT)
        _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()


# ── Output: flip/flop report ────────────────────────────────────────────────────
def write_flipflop_report(flipflop_results, aligned, n_reads, out_path):
    lines = [
        'FLIP/FLOP ANALYSIS REPORT',
        '=' * 60,
        '',
    ]
    for ff in flipflop_results:
        lines.append(f"ITR: {ff['itr_name']}")
        if not ff.get('found_arms'):
            lines.append('  No palindromic arms detected.')
            lines.append('')
            continue
        lines += [
            f"  Palindromic arm length: {ff['arm_len']} bp",
            f"  Arm A: pos {ff['arm_a_abs'][0]}–{ff['arm_a_abs'][1]}",
            f"  Arm B: pos {ff['arm_b_abs'][0]}–{ff['arm_b_abs'][1]}",
            f"  Inner loop: pos {ff['inner_abs'][0]}–{ff['inner_abs'][1]}",
            f"  Inner loop (FLIP/ref): {ff['inner_seq_flip']}",
            f"  Inner loop (FLOP/RC):  {ff['inner_seq_flop']}",
            f"  Informative positions: {ff['n_informative']}",
            f"  Mean flop fraction:    {ff['mean_flop_pct']:.2f}%",
            '',
            '  Per-position detail:',
            f"  {'Pos':>7}  {'Idx':>4}  {'Flip':>5}  {'Flop':>5}  {'Cov':>6}  {'Flop%':>6}",
        ]
        for p in ff['per_position']:
            lines.append(
                f"  {p['pos_1based']:>7}  {p['inner_idx']:>4}  "
                f"{p['flip_base']}={p['flip_count']:>4}  "
                f"{p['flop_base']}={p['flop_count']:>4}  "
                f"{p['coverage']:>6}  {p['flop_pct']:>5.1f}%"
            )
        lines.append('')

    Path(out_path).write_text('\n'.join(lines))


# ── Output: summary ─────────────────────────────────────────────────────────────
def write_summary(coverage, dev_pct, flag_mask, right_itr, left_itr,
                  aligned, n_reads, R, name, out_path):
    lines = [
        f'PILEUP SUMMARY: {name}',
        '=' * 60,
        f'Plasmid length:  {R:,} bp',
        f'Total reads:     {n_reads:,}',
        f'Aligned reads:   {aligned:,} ({100*aligned/n_reads:.1f}%)',
        f'Mean coverage:   {coverage.mean():.0f}x',
        f'Min coverage:    {coverage.min()}x',
        f'Max coverage:    {coverage.max()}x',
        f'Total flagged (≥10%, cov≥50): {flag_mask.sum()} positions',
        '',
    ]
    for label, lo, hi in (
        [('Right ITR', *right_itr)] if right_itr else [] +
        [('Left ITR',  *left_itr)]  if left_itr  else []
    ):
        mask = np.zeros(R, bool)
        mask[lo-1:hi] = True
        n_f = (flag_mask & mask).sum()
        m_d = dev_pct[mask].mean()
        lines.append(f'{label} ({lo}–{hi}, {hi-lo+1} bp):')
        lines.append(f'  Mean deviation: {m_d:.2f}%')
        lines.append(f'  Flagged:        {n_f}/{hi-lo+1}')
        lines.append('')

    bb = ~(np.zeros(R, bool))
    if right_itr:
        bb[right_itr[0]-1:right_itr[1]] = False
    if left_itr:
        bb[left_itr[0]-1:left_itr[1]] = False
    lines.append(f'Backbone ({bb.sum()} bp):')
    lines.append(f'  Mean deviation: {dev_pct[bb].mean():.2f}%')
    lines.append(f'  Flagged:        {(flag_mask & bb).sum()}')

    Path(out_path).write_text('\n'.join(lines))


# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='AAV plasmid ONT pileup analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--fastq',      required=True)
    parser.add_argument('--fasta',      required=True)
    parser.add_argument('--right-itr',  nargs=2, type=int, metavar=('START', 'END'))
    parser.add_argument('--left-itr',   nargs=2, type=int, metavar=('START', 'END'))
    parser.add_argument('--name',       default='plasmid')
    parser.add_argument('--output-dir', default='.')
    parser.add_argument('--flag-pct',   type=float, default=10.0)
    parser.add_argument('--min-cov',    type=int,   default=50)
    parser.add_argument('--no-flipflop', action='store_true')
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    right_itr = tuple(args.right_itr) if args.right_itr else None
    left_itr  = tuple(args.left_itr)  if args.left_itr  else None

    # ── Load inputs ──
    print(f'\n[1/5] Loading reference and reads...')
    ref = load_reference(args.fasta)
    R = len(ref)
    print(f'  Reference: {R:,} bp')
    print(f'  First 40 bp: {ref[:40]}')
    print(f'  (Verify this matches Azenta pos 1-40; if offset, rotate reference manually)')

    reads = load_reads(args.fastq)
    import numpy as np
    lens = np.array([len(r) for r in reads])
    print(f'  Reads: {len(reads):,} | median {int(np.median(lens))} bp | '
          f'min {lens.min()} | max {lens.max()}')

    # ── Build index & run pileup ──
    print(f'\n[2/5] Building k-mer index (k={K})...')
    idx = build_kmer_index(ref)

    print(f'\n[3/5] Aligning reads...')
    coverage, counts, aligned = run_pileup(reads, ref, idx, R)

    # ── Compute deviation ──
    cov_s = np.where(coverage > 0, coverage, 1)
    ref_cnt = np.array([counts[ref[i]][i] for i in range(R)], dtype=np.int32)
    dev_pct = (1 - ref_cnt / cov_s) * 100
    flag_mask = (dev_pct >= args.flag_pct) & (coverage >= args.min_cov)

    # ── Flip/flop analysis ──
    flipflop_results = []
    if not args.no_flipflop:
        print(f'\n[4/5] Analyzing flip/flop signal...')
        for label, itr in [('Right ITR', right_itr), ('Left ITR', left_itr)]:
            if itr is None:
                continue
            ff = analyze_flipflop(coverage, counts, ref,
                                  itr[0]-1, itr[1], label, cov_s)
            flipflop_results.append(ff)
            if ff.get('found_arms'):
                print(f'  {label}: inner loop flop fraction = {ff["mean_flop_pct"]:.2f}% '
                      f'({ff["n_informative"]} informative positions)')
            else:
                print(f'  {label}: no palindromic arms found')
    else:
        print(f'\n[4/5] Flip/flop analysis skipped.')

    # ── Write outputs ──
    print(f'\n[5/5] Writing outputs to {out_dir}/ ...')

    csv_path = out_dir / f'{args.name}_pileup.csv'
    df, _, _ = write_csv(coverage, counts, ref, R, right_itr, left_itr,
                         args.flag_pct, csv_path)
    print(f'  CSV: {csv_path}')

    fig1_path = out_dir / f'{args.name}_pileup.png'
    plot_full_plasmid(coverage, counts, ref, R, dev_pct, flag_mask,
                      right_itr, left_itr, args.name, fig1_path)
    print(f'  Full-plasmid figure: {fig1_path}')

    if right_itr or left_itr:
        fig2_path = out_dir / f'{args.name}_itr_detail.png'
        plot_itr_detail(coverage, counts, ref, R, dev_pct, flag_mask,
                        right_itr, left_itr, flipflop_results, args.name, fig2_path)
        print(f'  ITR detail figure: {fig2_path}')

    if flipflop_results:
        ff_path = out_dir / f'{args.name}_flipflop.txt'
        write_flipflop_report(flipflop_results, aligned, len(reads), ff_path)
        print(f'  Flip/flop report: {ff_path}')

    sum_path = out_dir / f'{args.name}_summary.txt'
    write_summary(coverage, dev_pct, flag_mask, right_itr, left_itr,
                  aligned, len(reads), R, args.name, sum_path)
    print(f'  Summary: {sum_path}')

    print(f'\nDone.\n')


if __name__ == '__main__':
    main()

