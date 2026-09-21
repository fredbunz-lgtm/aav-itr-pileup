# Addgene Survey Pipeline

This directory contains the scripts used to survey publicly available AAV
plasmids deposited at Addgene, in order to characterize the prevalence and
configuration of inverted terminal repeat (ITR) sequences — including
Delta-11 (Δ11) ITR variants — across the Addgene plasmid repository. Results
from this pipeline are reported in the associated bioRxiv submission.

**Note on step ordering:** The step order below was reconstructed from each
script's content, inputs/outputs, and file naming, since a record of the
exact original run sequence was not preserved. It reflects the best
available reconstruction of the logical pipeline flow, not a verified
execution log.

## Prerequisites

- Run all scripts from the repository root (e.g., `python addgene_survey/step_script.py`),
  so that relative paths to `results/` and `addgene_bulk/` resolve correctly.
- Requires the Addgene bulk plasmid data download, saved locally at
  `addgene_bulk/plasmids_with_sequences_download`. This file is large (~1.8GB)
  and is not included in this repository — it must be downloaded directly
  from Addgene.
- Requires Python packages: `biopython`, `ijson`.
- One script (`analyze_delta11_itr_ori_proximity.py`) fetches GenBank records
  live from Addgene's API and requires a personal API token set as the
  environment variable `ADDGENE_API_TOKEN`. Never commit this token or paste
  it into shell history, scripts, or version control.
- Output files are written to `results/`, which is excluded from version
  control (`.gitignore`) since these are large, regeneratable data files.
  The `results/*.tsv` files currently present locally are the authoritative
  record backing the published bioRxiv results.

## Pipeline Steps

1. **`inventory_addgene_bulk.py`**
   Reads the raw Addgene bulk JSON download and prints summary statistics
   (plasmid counts, sequence counts by group, sequence length ranges) to
   the console. No output file — purely informational, used to get
   oriented with the bulk dataset.

2. **`screen_addgene_metadata.py`**
   Performs an initial metadata-based screen of the bulk data to identify
   candidate AAV/ITR-relevant plasmids (e.g., by name/description keyword
   matching), without yet examining actual DNA sequences.

3. **`screen_addgene_metadata_refined.py`**
   A refined version of the metadata screen above, tightening the
   candidate-selection criteria based on lessons learned from the initial
   pass.

4. **`scan_addgene_candidate_itrs.py`** (uses **`itr_detection.py`** as a
   shared library)
   Runs actual ITR sequence detection against the candidate plasmids'
   sequences, using alignment-based detection logic defined in
   `itr_detection.py` (reference ITR sequences, aligner setup, hit
   classification, and collapsing overlapping hits). This is the core
   sequence-level ITR detection step.

5. **`summarize_addgene_itr_hits.py`**
   Summarizes the raw ITR detection hits from the previous step into a
   more condensed, human-readable report.

6. **`make_itr_positive_best_full_sequences.py`**
  For each ITR-positive plasmid, selects a single "best" sequence record
   to use in all downstream steps. Preference order: Addgene-verified full
   sequences, then user-submitted full sequences, then partial sequences
   (Addgene, then user), taking the longest sequence within the first
   non-empty preferred group. Writes
   `results/itr_positive_best_full_sequences.tsv`, which records which
   group/sequence_id was chosen for each plasmid — this becomes the
   reference lookup table for later steps that need to re-fetch specific
   sequences from the bulk JSON.
7. **`summarize_itr_length_isomer_configurations.py`**
   Summarizes ITR length and isomer configuration patterns across the
   ITR-positive plasmid set.

8. **`validate_delta11_exact_motifs.py`**
   Validates exact Delta-11 (Δ11) motif matches within the detected ITR
   sequences.

9. **`summarize_delta11_exact_plasmids.py`**
   Summarizes the set of plasmids confirmed to have exact Δ11 motif
   matches.

10. **`summarize_delta11_orientation_counts.py`**
    Summarizes orientation counts (e.g., forward vs. reverse-complement)
    of Δ11 motif matches.

11. **`review_addgene_delta11_total2_positions.py`**
    Reviews positional patterns of Δ11 "total2" configurations across
    plasmids, supporting closer inspection of specific configuration
    types.

12. **`scan_ori_motifs.py`** (single-FASTA tool) + **`scan_ori_motifs_batch.py`**
    (reconstructed batch wrapper)
    Performs a quick pre-screen for known plasmid origin-of-replication
    (ori) related sequence motifs (e.g., pUC/pMB1/colE1-associated
    elements), scanning both strands with circular wraparound handling.

    - `scan_ori_motifs.py` is the original single-FASTA-file command-line
      tool: `python scan_ori_motifs.py plasmid.fasta`, printing hits to
      stdout for one sequence at a time.
    - `scan_ori_motifs_batch.py` is a **reconstruction**, written during
      later repository reorganization. The original batch script that
      looped this scan over the full set of ITR-positive plasmids was not
      preserved and could not be recovered. This reconstruction reuses the
      original motif dictionary and matching logic unchanged (imported
      directly from `scan_ori_motifs.py`), and adds new batch I/O code to:
      read `results/itr_positive_best_full_sequences.tsv`, look up each
      plasmid's already-selected "best" sequence in the bulk JSON, run the
      motif scan, and write one row per hit to
      `results/ori_motif_scan_results.tsv`. It was verified to run
      successfully end-to-end (7,041/7,041 targets matched, 27,351 hits
      written, ~9 seconds), but it was **not** used to produce the results
      reported in the bioRxiv submission — it is provided so the full
      pipeline remains runnable going forward.

13. **`analyze_delta11_itr_ori_proximity.py`**
    Analyzes the proximity between Δ11 ITR motifs and ori-related motifs
    within each plasmid. This script fetches GenBank records live from
    Addgene's API and therefore requires a valid `ADDGENE_API_TOKEN`
    environment variable to run.

14. **`compare_itr_delta11_plasmids.py`**
    Performs a final comparison across the ITR/Δ11-characterized plasmid
    set, consolidating results from the earlier summary and analysis
    steps.
