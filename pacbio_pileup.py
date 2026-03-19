#!/usr/bin/env python3
"""
PacBio Subreads Pileup Analysis — Linear Fragment
===================================================
Aligns PacBio subreads BAM to a LINEAR reference using mappy,
generates per-position pileup, and performs flip/flop analysis.

Designed for gel-isolated restriction fragments where the reference
is linear (no circular wrapping needed).

Usage:
    python pacbio_pileup.py \
        --bam pSM620.subreads.bam \
        --fasta "pSM620.7 fragment.fasta" \
        --right-itr 42 185 \
        --left-itr 4583 4724 \
        --name pSM620_fragment \
        --output-dir ./pacbio_results/

Arguments:
    --bam           PacBio subreads BAM (unaligned)
    --fasta         Linear reference FASTA
    --right-itr     1-based inclusive coords of right ITR
    --left-itr      1-based inclusive coords of left ITR
    --name          Output file prefix
    --output-dir    Output directory
    --flag-pct      Deviation threshold for flagging (default 10.0)
    --min-len       Minimum subread length to align (default 500)
    --max-reads     Subsample to this many reads (default 200000)

Outputs:
    {name}_pileup.csv       Per-position pileup table
    {name}_pileup.png       Full-fragment figure
    {name}_itr_detail.png   ITR detail figure
    {name}_flipflop.txt     Flip/flop report
    {name}_summary.txt      Run summary
"""

import argparse
import collections
import random
import time
from pathlib import Path

import mappy
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────
BASES     = 'ATGC'
MIN_LEN   = 500
MAX_READS = 200000
MIN_MAPQ  = 0

# ── Utilities ──────────────────────────────────────────────────────────────────
def rc(s):
    return s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def load_reference(fasta_path):
    lines = Path(fasta_path).read_text().strip().splitlines()
    seq = ''.join(l for l in lines if not l.startswith('>')).upper()
    iupac = {c: 'A' for c in 'MRWSYKVHDBN'}
    for bad, good in iupac.items():
        if bad in seq:
            print(f'  WARNING: replacing ambiguous base {bad} with {good}')
            seq = seq.replace(bad, good)
    return seq


# ── Pileup ─────────────────────────────────────────────────────────────────────
def run_pileup(bam_path, ref, R, min_len=MIN_LEN, max_reads=MAX_READS, verbose=True):
    counts   = {b: np.zeros(R, dtype=np.int32) for b in BASES}
    coverage = np.zeros(R, dtype=np.int32)

    aligner = mappy.Aligner(seq=ref, preset='map-pb', best_n=1)
    if not aligner:
        raise RuntimeError('Failed to build mappy index')

    if verbose:
        print('  Counting reads...', flush=True)
    bam   = pysam.AlignmentFile(bam_path, 'rb', check_sq=False)
    total = bam.count(until_eof=True)
    bam.close()
    if verbose:
        print(f'  Total subreads in BAM: {total:,}')

    random.seed(42)
    keep_prob = min(1.0, max_reads / total) if total > max_reads else 1.0
    if keep_prob < 1.0 and verbose:
        print(f'  Subsampling: keeping {keep_prob:.3f} of reads (~{max_reads:,})')

    aligned = skipped_short = skipped_nomatch = processed = 0
    t0 = time.time()

    bam = pysam.AlignmentFile(bam_path, 'rb', check_sq=False)
    for read in bam.fetch(until_eof=True):
        if keep_prob < 1.0 and random.random() > keep_prob:
            continue
        seq = read.query_sequence
        if seq is None or len(seq) < min_len:
            skipped_short += 1
            continue

        processed += 1
        if verbose and processed % 50000 == 0:
            print(f'  Processed {processed:,} reads ({time.time()-t0:.1f}s)', flush=True)

        best_hit   = None
        best_score = -1
        for hit in aligner.map(seq):
            if hit.mapq >= MIN_MAPQ and hit.blen > best_score:
                best_score = hit.blen
                best_hit   = hit

        if best_hit is None:
            skipped_nomatch += 1
            continue

        aligned += 1
        query  = seq if best_hit.strand == 1 else rc(seq)
        r_pos  = best_hit.r_st
        q_pos  = best_hit.q_st

        for op, length in best_hit.cigar:
            if op in (0, 7, 8):
                for i in range(length):
                    rp = r_pos + i
                    qp = q_pos + i
                    if 0 <= rp < R and qp < len(query):
                        base = query[qp]
                        if base in counts:
                            counts[base][rp] += 1
                            coverage[rp]     += 1
                r_pos += length
                q_pos += length
            elif op == 1:
                q_pos += length
            elif op in (2, 3):
                r_pos += length
            elif op in (4, 5):
                q_pos += length

    bam.close()

    if verbose:
        print(f'  Aligned: {aligned:,}/{processed:,} '
              f'({100*aligned/max(processed,1):.1f}%) in {time.time()-t0:.1f}s')
        print(f'  Skipped (too short <{min_len}bp): {skipped_short:,}')
        print(f'  Skipped (no alignment): {skipped_nomatch:,}')
        print(f'  Mean coverage: {coverage.mean():.0f}x  '
              f'min: {coverage.min()}x  max: {coverage.max()}x')

    return coverage, counts, aligned, processed


# ── Flip/flop ──────────────────────────────────────────────────────────────────
def find_best_palindrome(s, min_arm=25, max_arm=45, max_mm=3):
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


def analyze_flipflop(coverage, counts, ref, itr_lo, itr_hi, itr_name):
    itr_seq          = ref[itr_lo:itr_hi]
    armlen, sa, sb, mm = find_best_palindrome(itr_seq)
    if armlen == 0:
        return {'found_arms': False, 'itr_name': itr_name}

    inner_seq    = itr_seq[sa + armlen:sb]
    flop_seq     = rc(inner_seq)
    inner_lo_abs = itr_lo + sa + armlen
    inner_hi_abs = itr_lo + sb

    flop_counts = []
    total_covs  = []
    results_per_pos = []

    for rel_i, (flip_b, flop_b) in enumerate(zip(inner_seq, flop_seq)):
        if flip_b == flop_b:
            continue
        abs_i  = inner_lo_abs + rel_i
        bases  = {b: counts[b][abs_i] for b in BASES}
        cov    = coverage[abs_i]
        flop_n = bases[flop_b]
        flip_n = bases[flip_b]
        flop_pct = 100 * flop_n / cov if cov > 0 else 0.0
        flop_counts.append(flop_n)
        total_covs.append(cov)
        results_per_pos.append({
            'pos_1based': abs_i + 1,
            'inner_idx':  rel_i,
            'flip_base':  flip_b,
            'flop_base':  flop_b,
            'flip_count': flip_n,
            'flop_count': flop_n,
            'coverage':   cov,
            'flop_pct':   round(flop_pct, 2),
        })

    mean_flop = 100 * sum(flop_counts) / max(sum(total_covs), 1)

    return {
        'found_arms':     True,
        'itr_name':       itr_name,
        'arm_len':        armlen,
        'arm_mm':         mm,
        'inner_abs':      (inner_lo_abs + 1, inner_hi_abs),
        'inner_seq_flip': inner_seq,
        'inner_seq_flop': flop_seq,
        'n_informative':  len(results_per_pos),
        'mean_flop_pct':  round(mean_flop, 2),
        'per_position':   results_per_pos,
    }


# ── Figures ────────────────────────────────────────────────────────────────────
BG        = '#0d1117'; PANEL = '#161b22'; GRID = '#21262d'; TEXT = 'white'
BCOL      = {'A': '#3fb950', 'T': '#f78166', 'G': '#d2a8ff', 'C': '#ffa657'}
DEV_COLOR = '#58a6ff'; FLAG_COLOR = '#f85149'
RITR_COL  = '#ffa657'; LITR_COL   = '#3fb950'


def _style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TEXT, labelsize=8)
    for lbl in (ax.yaxis.label, ax.xaxis.label, ax.title):
        lbl.set_color(TEXT)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.grid(axis='y', color=GRID, lw=0.5, alpha=0.6)


def plot_full(coverage, counts, ref, R, dev_pct, flag_mask,
              right_itr, left_itr, name, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    pos   = np.arange(1, R + 1)

    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True)
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} PacBio Pileup — Fragment ({R:,} bp)',
                 color=TEXT, fontsize=13, fontweight='bold')

    ax = axes[0]
    ax.fill_between(pos, coverage, color=DEV_COLOR, alpha=0.6)
    if right_itr:
        ax.axvspan(*right_itr, color=RITR_COL, alpha=0.25, label='Right ITR')
    if left_itr:
        ax.axvspan(*left_itr, color=LITR_COL, alpha=0.25, label='Left ITR')
    ax.set_ylabel('Coverage (x)', color=TEXT)
    ax.set_title('Coverage', color=TEXT)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    ax  = axes[1]
    bot = np.zeros(R)
    for base, color in BCOL.items():
        frac = counts[base] / cov_s
        ax.bar(pos, frac, bottom=bot, color=color, label=base, width=1, alpha=0.9)
        bot += frac
    if right_itr:
        ax.axvspan(*right_itr, color=RITR_COL, alpha=0.2)
    if left_itr:
        ax.axvspan(*left_itr, color=LITR_COL, alpha=0.2)
    ax.set_ylabel('Base fraction', color=TEXT)
    ax.set_title('Base composition', color=TEXT)
    ax.legend(loc='lower right', ncol=4, fontsize=8, facecolor=PANEL, labelcolor=TEXT)
    _style(ax)

    ax        = axes[2]
    ritr_mask = np.zeros(R, bool)
    litr_mask = np.zeros(R, bool)
    if right_itr:
        ritr_mask[right_itr[0] - 1:right_itr[1]] = True
    if left_itr:
        litr_mask[left_itr[0] - 1:left_itr[1]] = True

    ax.bar(pos[~flag_mask], dev_pct[~flag_mask], width=1, color=DEV_COLOR, alpha=0.4)
    for mask, color, label in [
        (flag_mask & ritr_mask,              RITR_COL,   'Right ITR'),
        (flag_mask & litr_mask,              LITR_COL,   'Left ITR'),
        (flag_mask & ~ritr_mask & ~litr_mask, FLAG_COLOR, 'Other'),
    ]:
        if mask.any():
            ax.bar(pos[mask], dev_pct[mask], width=1, color=color, alpha=0.9, label=label)
    ax.axhline(10.0, color=FLAG_COLOR, ls='--', lw=1, alpha=0.7)
    ax.set_ylabel('Non-ref bases (%)', color=TEXT)
    ax.set_xlabel('Position', color=TEXT)
    ax.set_title('Deviation from reference', color=TEXT)
    ax.legend(fontsize=8, facecolor=PANEL, labelcolor=TEXT, ncol=4)
    _style(ax)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved: {out_path}')


def plot_itr_detail(coverage, counts, ref, R, dev_pct, flag_mask,
                    right_itr, left_itr, flipflop_results, name, out_path):
    cov_s = np.where(coverage > 0, coverage, 1)
    pos   = np.arange(1, R + 1)

    itrs = []
    if right_itr:
        itrs.append((right_itr, 'Right ITR', RITR_COL))
    if left_itr:
        itrs.append((left_itr,  'Left ITR',  LITR_COL))

    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.patch.set_facecolor(BG)
    fig.suptitle(f'{name} — ITR Detail (PacBio)',
                 color=TEXT, fontsize=13, fontweight='bold')

    for col_idx, (itr, label, color) in enumerate(itrs[:2]):
        lo   = max(1, itr[0] - 30)
        hi   = min(R, itr[1] + 10)
        mask = (pos >= lo) & (pos <= hi)
        p_sub    = pos[mask]
        dev_sub  = dev_pct[mask]
        cov_sub  = cov_s[mask]
        flag_sub = flag_mask[mask]
        itr_sub  = (p_sub >= itr[0]) & (p_sub <= itr[1])

        ff = next((r for r in flipflop_results
                   if r.get('found_arms') and label in r['itr_name']), None)

        ax  = axes[0][col_idx]
        bot = np.zeros(mask.sum())
        for base, bc in BCOL.items():
            frac = counts[base][mask] / cov_sub
            ax.bar(p_sub, frac, bottom=bot, color=bc, label=base, width=1, alpha=0.9)
            bot += frac
        ax.axvline(itr[0], color=color, lw=2, ls='--', alpha=0.8, label='ITR boundary')
        ax.axvline(itr[1], color=color, lw=2, ls='--', alpha=0.8)
        if ff:
            ax.axvspan(ff['inner_abs'][0], ff['inner_abs'][1],
                       color='#ffffff', alpha=0.1,
                       label=f'Inner loop ({ff["mean_flop_pct"]:.1f}% flop)')
        ax.set_title(f'{label} — Base Composition', color=TEXT, fontweight='bold')
        ax.set_ylabel('Fraction', color=TEXT)
        ax.legend(ncol=5, fontsize=7, facecolor=PANEL, labelcolor=TEXT)
        _style(ax)

        ax = axes[1][col_idx]
        ax.bar(p_sub[~flag_sub], dev_sub[~flag_sub], width=1, color=DEV_COLOR, alpha=0.6)
        ax.bar(p_sub[flag_sub],  dev_sub[flag_sub],  width=1, color=color,     alpha=0.9)
        ax.axhline(10.0, color=FLAG_COLOR, ls='--', lw=1, alpha=0.7)
        ax.axvline(itr[0], color=color, lw=2, ls='--', alpha=0.8)
        ax.axvline(itr[1], color=color, lw=2, ls='--', alpha=0.8)
        if ff:
            ax.axvspan(ff['inner_abs'][0], ff['inner_abs'][1],
                       color='#ffffff', alpha=0.1)

        n_flag   = (flag_sub & itr_sub).sum()
        mean_dev = dev_sub[itr_sub].mean() if itr_sub.sum() > 0 else 0
        info = f'Mean: {mean_dev:.2f}%\nFlagged: {n_flag}/{itr_sub.sum()}'
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
    print(f'  Saved: {out_path}')


# ── Text outputs ───────────────────────────────────────────────────────────────
def write_flipflop_report(flipflop_results, out_path):
    lines = ['FLIP/FLOP ANALYSIS REPORT (PacBio)', '=' * 60, '']
    for ff in flipflop_results:
        lines.append(f"ITR: {ff['itr_name']}")
        if not ff.get('found_arms'):
            lines += ['  No palindromic arms detected.', '']
            continue
        lines += [
            f"  Arm length: {ff['arm_len']} bp ({ff['arm_mm']} mismatches)",
            f"  Inner loop: pos {ff['inner_abs'][0]}-{ff['inner_abs'][1]}",
            f"  Inner loop (FLIP/ref): {ff['inner_seq_flip']}",
            f"  Inner loop (FLOP/RC):  {ff['inner_seq_flop']}",
            f"  Informative positions: {ff['n_informative']}",
            f"  Mean flop fraction:    {ff['mean_flop_pct']:.2f}%",
            '',
            f"  {'Pos':>7}  {'Idx':>4}  {'Flip':>6}  {'Flop':>6}  {'Cov':>6}  {'Flop%':>6}",
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
    print(f'  Saved: {out_path}')


def write_csv(coverage, counts, ref, R, right_itr, left_itr, flag_pct, out_path):
    cov_s   = np.where(coverage > 0, coverage, 1)
    ref_cnt = np.array([counts[ref[i]][i] if ref[i] in counts else 0
                        for i in range(R)], dtype=np.int32)
    dev_pct   = (1 - ref_cnt / cov_s) * 100
    flag_mask = (dev_pct >= flag_pct) & (coverage >= 50)

    def region(p):
        if right_itr and right_itr[0] <= p <= right_itr[1]:
            return 'Right_ITR'
        if left_itr  and left_itr[0]  <= p <= left_itr[1]:
            return 'Left_ITR'
        return 'Fragment'

    pos_1 = np.arange(1, R + 1)
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
    print(f'  Saved: {out_path}')
    return df, dev_pct, flag_mask


def write_summary(coverage, dev_pct, flag_mask, right_itr, left_itr,
                  aligned, processed, R, name, out_path):
    lines = [
        f'PACBIO FRAGMENT PILEUP SUMMARY: {name}',
        '=' * 60,
        f'Fragment length: {R:,} bp',
        f'Processed reads: {processed:,}',
        f'Aligned reads:   {aligned:,} ({100*aligned/max(processed,1):.1f}%)',
        f'Mean coverage:   {coverage.mean():.0f}x',
        f'Min coverage:    {coverage.min()}x',
        f'Max coverage:    {coverage.max()}x',
        f'Flagged (>=10%, cov>=50): {flag_mask.sum()} positions',
        '',
    ]
    for label, itr in [('Right ITR', right_itr), ('Left ITR', left_itr)]:
        if itr is None:
            continue
        lo, hi = itr
        m = np.zeros(R, bool)
        m[lo - 1:hi] = True
        lines += [
            f'{label} ({lo}-{hi}, {hi-lo+1} bp):',
            f'  Mean deviation: {dev_pct[m].mean():.2f}%',
            f'  Flagged:        {(flag_mask & m).sum()}/{m.sum()}',
            '',
        ]
    Path(out_path).write_text('\n'.join(lines))
    print(f'  Saved: {out_path}')


# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='PacBio subreads pileup — linear fragment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--bam',        required=True)
    parser.add_argument('--fasta',      required=True)
    parser.add_argument('--right-itr',  nargs=2, type=int, metavar=('START', 'END'))
    parser.add_argument('--left-itr',   nargs=2, type=int, metavar=('START', 'END'))
    parser.add_argument('--name',       default='pacbio_fragment')
    parser.add_argument('--output-dir', default='.')
    parser.add_argument('--flag-pct',   type=float, default=10.0)
    parser.add_argument('--min-len',    type=int,   default=MIN_LEN)
    parser.add_argument('--max-reads',  type=int,   default=MAX_READS)
    args = parser.parse_args()

    out_dir   = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    right_itr = tuple(args.right_itr) if args.right_itr else None
    left_itr  = tuple(args.left_itr)  if args.left_itr  else None

    print(f'\n[1/5] Loading reference...')
    ref = load_reference(args.fasta)
    R   = len(ref)
    print(f'  Fragment length: {R:,} bp')
    print(f'  First 40 bp: {ref[:40]}')

    print(f'\n[2/5] Building mappy index (map-pb preset)...')

    print(f'\n[3/5] Aligning subreads...')
    coverage, counts, aligned, processed = run_pileup(
        args.bam, ref, R,
        min_len=args.min_len,
        max_reads=args.max_reads,
    )

    print(f'\n[4/5] Analyzing flip/flop signal...')
    flipflop_results = []
    for label, itr in [('Right ITR', right_itr), ('Left ITR', left_itr)]:
        if itr is None:
            continue
        ff = analyze_flipflop(coverage, counts, ref, itr[0] - 1, itr[1], label)
        flipflop_results.append(ff)
        if ff.get('found_arms'):
            print(f'  {label}: {ff["mean_flop_pct"]:.2f}% flop '
                  f'({ff["n_informative"]} informative positions)')
        else:
            print(f'  {label}: no palindromic arms found')

    print(f'\n[5/5] Writing outputs to {out_dir}/ ...')
    csv_path = out_dir / f'{args.name}_pileup.csv'
    df, dev_pct, flag_mask = write_csv(
        coverage, counts, ref, R, right_itr, left_itr, args.flag_pct, csv_path)

    plot_full(coverage, counts, ref, R, dev_pct, flag_mask,
              right_itr, left_itr, args.name,
              out_dir / f'{args.name}_pileup.png')

    if right_itr or left_itr:
        plot_itr_detail(coverage, counts, ref, R, dev_pct, flag_mask,
                        right_itr, left_itr, flipflop_results, args.name,
                        out_dir / f'{args.name}_itr_detail.png')

    write_flipflop_report(flipflop_results,
                          out_dir / f'{args.name}_flipflop.txt')
    write_summary(coverage, dev_pct, flag_mask, right_itr, left_itr,
                  aligned, processed, R, args.name,
                  out_dir / f'{args.name}_summary.txt')

    print(f'\nDone.\n')


if __name__ == '__main__':
    main()
