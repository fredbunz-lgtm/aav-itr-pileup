# AAV Plasmid ITR Heterogeneity Pipeline

Custom Python pipelines for per-position basecall heterogeneity analysis
of AAV plasmid vectors from Oxford Nanopore Technologies (ONT) raw FASTQ
reads, and PacBio subreads BAM files. Developed for the analysis described
in:

Jaskula-Ranga et al. "Stable propagation of inverted terminal repeats in
adeno-associated virus plasmid vectors." (2025)

SRA accession: PRJNA1439403

---

## Overview

These scripts align raw nanopore or PacBio reads to a circular or linear
plasmid reference without external aligners (no minimap2 or samtools
required for the ONT pipeline). At each reference position, per-base
counts and deviation from the reference are tabulated. Flip/flop isomer
fractions are estimated at ITR inner loop positions by identifying the
palindromic outer stem and comparing observed minor alleles to the
predicted flop sequence.

---

## Scripts

### aav_plasmid_workflow.py
ONT pileup pipeline for circular plasmid references.
- No external dependencies beyond Python scientific stack
- Circular-aware seed-chain alignment with banded Needleman-Wunsch
- Handles reads spanning the reference linearization point
- Automatic flip/flop inner loop quantification
- Outputs: per-position CSV, full-plasmid figure, ITR detail figure,
  flip/flop report, run summary

### pacbio_pileup.py
PacBio subreads pileup pipeline for linear fragment references.
- Uses mappy (minimap2 Python bindings) for alignment
- Designed for gel-isolated restriction fragments
- Same output format as aav_plasmid_workflow.py for direct comparison

---

## Dependencies

### ONT pipeline (aav_plasmid_workflow.py)
```
numpy
pandas
matplotlib
```

### PacBio pipeline (pacbio_pileup.py)
```
numpy
pandas
matplotlib
mappy
pysam
```

Install all dependencies with:
```bash
pip install numpy pandas matplotlib mappy pysam
```

---

## Usage

### ONT circular plasmid analysis
```bash
python aav_plasmid_workflow.py \
  --fastq reads.fastq.gz \
  --fasta reference.fasta \
  --left-itr <start> <end> \
  --right-itr <start> <end> \
  --name sample_name \
  --output-dir ./results/
```

### PacBio linear fragment analysis
```bash
python pacbio_pileup.py \
  --bam subreads.bam \
  --fasta fragment_reference.fasta \
  --left-itr <start> <end> \
  --right-itr <start> <end> \
  --name sample_name \
  --output-dir ./results/
```

ITR coordinates are **1-based inclusive**, matching standard genome
browser conventions. If using a BED file (0-based half-open), add 1 to
the start coordinate and leave the end coordinate unchanged.

---

## Commands used in the paper

### pAAVcis
```bash
python aav_plasmid_workflow.py \
  --fastq psAAVcis.fastq \
  --fasta pAAVcis.fasta \
  --left-itr 423 566 \
  --right-itr 2177 2320 \
  --name pAAVcis \
  --output-dir ./pAAVcis_results/
```

### pSM620 (barcode07)
```bash
python aav_plasmid_workflow.py \
  --fastq barcode07_pSM620.fastq \
  --fasta pSM620_ref_rotated.fasta \
  --left-itr 2081 2223 \
  --right-itr 6614 6747 \
  --name pSM620_barcode07 \
  --output-dir ./pSM620_barcode07_results/
```

### pSM620 (barcode41)
```bash
python aav_plasmid_workflow.py \
  --fastq barcode41_pSM620-1.fastq \
  --fasta pSM620_ref_rotated.fasta \
  --left-itr 2081 2223 \
  --right-itr 6614 6747 \
  --name pSM620_barcode41 \
  --output-dir ./pSM620_barcode41_results/
```

### pSM630 (72h serial culture)
```bash
python aav_plasmid_workflow.py \
  --fastq pSM630.fastq \
  --fasta pSM620_ref_rotated.fasta \
  --left-itr 2081 2223 \
  --right-itr 6614 6747 \
  --name pSM630_72h \
  --output-dir ./pSM630_72h_results/
```

### Clone4 (pSM620.4)
```bash
python aav_plasmid_workflow.py \
  --fastq Clone4.fastq \
  --fasta pSM620_ref_rotated.fasta \
  --left-itr 2081 2223 \
  --right-itr 6614 6747 \
  --name Clone4 \
  --output-dir ./Clone4_results/
```

### pBR-AAV2(-)
```bash
python aav_plasmid_workflow.py \
  --fastq barcode08_pBR-AAV2.fastq \
  --fasta pAAV2-AAV2-.fasta \
  --left-itr 1 145 \
  --right-itr 4535 4679 \
  --name pBR-AAV2 \
  --output-dir ./pBR-AAV2_results/
```

### pAAV2ST
```bash
python aav_plasmid_workflow.py \
  --fastq barcode13_Spacer-3.fastq \
  --fasta pAAV2ST.fasta \
  --left-itr 2216 2360 \
  --right-itr 3970 4114 \
  --name pAAV2ST_spacer3 \
  --output-dir ./pAAV2ST_spacer3_results/
```

### pAAV2ST-Pcsk9sg1
```bash
python aav_plasmid_workflow.py \
  --fastq pAAV2ST-Pcsk9sg1.fastq \
  --fasta pAAV2ST-Pcsk9sg1.fasta \
  --left-itr 1 145 \
  --right-itr 4607 4750 \
  --name pAAV2ST-Pcsk9sg1 \
  --output-dir ./pAAV2ST_results/
```

---

## Important notes

**Reference linearization:** For circular plasmid references where the
linearization point falls within or adjacent to an ITR, wrap-around
alignment artifacts can inflate apparent heterogeneity at ITR positions.
Rotate the reference to place the linearization point >2 kb from both
ITRs before running the pipeline. The pSM620_ref_rotated.fasta file
provided in this repository has been rotated accordingly (cut at original
position 7000).

**Coordinate offset:** Some reference FASTAs may be linearized at a
different position than the Azenta Plasmid EZ coordinate frame. Verify
that the first 40 bp printed at runtime matches the expected reference
sequence. If there is an offset N, rotate the reference:
```python
ref = ref_orig[N:] + ref_orig[:N]
```
and adjust ITR coordinates accordingly.

**IUPAC ambiguity codes:** The pipeline replaces non-ACGT bases with A.
Check the runtime output for warnings about ambiguous bases in the
reference.

**Dam methylation artifacts:** GATC sites produce systematic ONT
basecall errors (10-40% apparent deviation) due to N6-methyladenine.
Identify by GATC sequence context; do not interpret as genuine sequence
heterogeneity.

---

## Output files

| File | Description |
|------|-------------|
| `{name}_pileup.csv` | Per-position table (position, coverage, A/T/G/C counts, deviation %) |
| `{name}_pileup.png` | Full-plasmid 3-panel figure (coverage, base composition, deviation) |
| `{name}_itr_detail.png` | ITR detail 2×2 figure with flip/flop inner loop highlighted |
| `{name}_flipflop.txt` | Per-position flip/flop quantification report |
| `{name}_summary.txt` | Run summary statistics |

---

## License

MIT License. See LICENSE file for details.
