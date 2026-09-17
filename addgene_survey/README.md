# Addgene AAV ITR Analysis Scripts

This folder contains the main scripts used for the Addgene AAV ITR analysis workflow.

## Pipeline overview

The numbered scripts are intended to be run in order:

- `00_inventory_addgene_bulk.py` — collect / inventory Addgene bulk metadata
- `01_screen_addgene_metadata_refined.py` — screen and refine candidate metadata
- `02_scan_addgene_candidate_itrs.py` — scan candidate sequences for ITR-like hits
- `03_summarize_addgene_itr_hits.py` — summarize ITR scan results
- `04_select_itr_positive_representative_sequences.py` — select representative positive sequences
- `05_validate_delta11_exact_motifs.py` — validate exact delta11 motifs
- `06_summarize_delta11_exact_plasmids.py` — summarize exact motif plasmids
- `07_summarize_delta11_orientation_counts.py` — summarize orientation counts
- `08_compare_itr_delta11_plasmids.py` — compare ITR and delta11 plasmid sets
- `09_analyze_delta11_itr_ori_proximity.py` — analyze proximity between ITR and ori motifs
- `10_summarize_itr_length_isomer_configurations.py` — summarize ITR length/isomer configurations

## Supporting code

- `itr_detection.py` contains reusable ITR/origin detection logic shared across scripts.
- `qc/` contains quality-control scripts, including:
  - `qc_inspect_two_delta11_motif_plasmids.py`

## Notes

- Raw Addgene GenBank downloads (`*.gb`) are not tracked in GitHub.
- Large generated outputs are stored outside the repo or excluded via `.gitignore`.
- The `archive/old_scripts/` folder contains older script versions kept only for reference.
