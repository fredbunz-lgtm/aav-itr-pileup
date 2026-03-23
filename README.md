# AAV ITR Pileup Pipeline

Custom Python pipelines for nanopore sequencing analysis of AAV plasmid vectors and encapsidated viral genomes. Developed to support:

> Jaskula-Ranga V, Li S, Pandya A, Bunz F. *Stable propagation of inverted terminal repeats in adeno-associated virus plasmid vectors.* bioRxiv (2026).

---

## Overview

This repository contains three analysis scripts:

| Script | Input | Purpose |
|--------|-------|---------|
| `aav_plasmid_workflow.py` | Plasmid ONT FASTQ | Per-position sequence heterogeneity and ITR integrity in plasmid stocks |
| `aav_genome_qc.py` | Viral genome ONT FASTQ | Packaged genome structure, coverage, and contamination |
| `pacbio_pileup.py` | PacBio CLR BAM | Per-position pileup from PacBio subreads |

All scripts are pure Python, require no external aligners, and produce publication-quality figures alongside CSV outputs suitable for downstream analysis.

---

## Installation

Requires Python 3.8+ and the following packages:

```bash
pip install mappy numpy pandas matplotlib pysam
```

Clone the repository:

```bash
git clone https://github.com/[username]/aav-itr-pileup.git
cd aav-itr-pileup
```

### Bundled reference sequences

The `references/` folder contains default reference sequences used by `aav_genome_qc.py` for contamination analysis:

| File | Description |
|------|-------------|
| `references/ref_helper_pHelper.fasta` | pHelper plasmid (Agilent/Stratagene, accession AF369965.1, 11,635 bp). Contains adenovirus 5 helper functions (E2A, E4, VA RNA) required for AAV packaging. Standard across virtually all helper-free AAV production systems. |
| `references/ref_host_ecoli_K12.fasta` | *E. coli* K-12 MG1655 genome (RefSeq NC_000913.3, 4.6 Mb). Suitable for detection of encapsidated host DNA from any common laboratory *E. coli* strain. |

These defaults can be overridden at runtime using `--helper` and `--host`. The rep-cap plasmid reference must always be provided by the user as it varies by AAV serotype and construct.

---

## 1. Plasmid ITR Analysis (`aav_plasmid_workflow.py`)

Analyzes ONT sequencing data from AAV plasmid DNA preparations. Uses a custom circular-aware seed-chain alignment algorithm (no external aligners required). Quantifies flip/flop isomer fractions at ITR inner loop positions with single-molecule resolution.

### Usage

```bash
python aav_plasmid_workflow.py \
    --fastq psAAVcis.fastq \
    --fasta pAAVcis.fasta \
    --left-itr 423 566 \
    --right-itr 2177 2320 \
    --name pAAVcis \
    --output-dir ./pAAVcis_results/
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--fastq` | ONT FASTQ file (can be gzipped) |
| `--fasta` | Plasmid reference FASTA (circular) |
| `--left-itr` | Left ITR start and end positions (1-based inclusive) |
| `--right-itr` | Right ITR start and end positions (1-based inclusive) |
| `--name` | Sample name for output files |
| `--output-dir` | Output directory (created if absent) |
| `--min-base-qual` | Minimum per-base Phred quality score (bases below threshold excluded from pileup, default 20) |

### Outputs

| File | Description |
|------|-------------|
| `{name}_pileup.csv` | Per-position base counts and deviation from reference |
| `{name}_flipflop.txt` | Flip/flop analysis report with per-position detail |
| `{name}_summary.txt` | Alignment and coverage summary |
| `{name}_pileup.png` | Coverage and deviation figure |
| `{name}_itr_detail.png` | ITR detail with flip/flop annotation |
| `{name}_flipflop.png` | Flip/flop summary figure |

### Run commands used in the paper

```bash
# pAAVcis
python aav_plasmid_workflow.py \
    --fastq "psAAVcis.fastq" --fasta "pAAVcis.fasta" \
    --left-itr 423 566 --right-itr 2177 2320 \
    --name pAAVcis --output-dir ./pAAVcis_results/

# pSM620 (barcode07)
python aav_plasmid_workflow.py \
    --fastq "barcode07_pSM620.fastq" --fasta "pSM620_ref_rotated.fasta" \
    --left-itr 2081 2223 --right-itr 6614 6747 \
    --name pSM620_barcode07 --output-dir ./pSM620_barcode07_results/

# pSM630 (72h serial culture)
python aav_plasmid_workflow.py \
    --fastq "pSM630.fastq" --fasta "pSM620_ref_rotated.fasta" \
    --left-itr 2081 2223 --right-itr 6614 6747 \
    --name pSM630_72h --output-dir ./pSM630_72h_results/

# pBR-AAV2(-)
python aav_plasmid_workflow.py \
    --fastq "barcode08_pBR-AAV2.fastq" --fasta "pAAV2-AAV2-.fasta" \
    --left-itr 1 145 --right-itr 4607 4750 \
    --name pBR-AAV2 --output-dir ./pBR-AAV2_results/

# pAAV2ST (Spacer-3, GFP entry vector)
python aav_plasmid_workflow.py \
    --fastq "barcode13_Spacer-3.fastq" --fasta "pAAV2ST.fasta" \
    --left-itr 2216 2360 --right-itr 3970 4114 \
    --name pAAV2ST_spacer3 --output-dir ./pAAV2ST_spacer3_results/

# pAAV2ST-Pcsk9sg1
python aav_plasmid_workflow.py \
    --fastq "pAAV2ST-Pcsk9sg1.fastq" --fasta "pAAV2ST-Pcsk9sg1.fasta" \
    --left-itr 1 145 --right-itr 4607 4750 \
    --name pAAV2ST-Pcsk9sg1 --output-dir ./pAAV2ST_results/
```

### Note on reference preparation

The pSM620 reference contains an IUPAC ambiguity code (M) at position 794 which was replaced with A prior to analysis. The reference was also rotated to avoid a linearization artifact at the ITR boundary. The rotated reference (`pSM620_ref_rotated.fasta`) is deposited in the SRA alongside the sequencing data.

---

## 2. Packaged Genome QC (`aav_genome_qc.py`)

Analyzes ONT sequencing data from DNase I-treated, purified AAV viral preparations. Automatically detects ITR positions from the cis plasmid FASTA — no BED file is required. Performs multi-reference contamination analysis to quantify encapsidated transgene, helper plasmid, rep-cap plasmid, host genomic DNA, and reverse-packaged backbone sequences.

### Key features

- **Automatic ITR detection** from any AAV2-derived cis plasmid FASTA
- **No BED file required** — ITR coordinates are determined algorithmically from conserved AAV2 palindromic arm sequences
- **Contamination analysis** against helper, rep-cap, and host references
- **Backbone exclusion** — pileup restricted to ITR-to-ITR insert; backbone-mapping reads reported separately as reverse packaging
- **pAAV-MCS deletion warning** — automatically detects the 11-bp deletion (AAAGCCCGGGC) present in pAAV-MCS and all derivatives
- **Bundled default references** — pHelper (Agilent V005569) and *E. coli* K-12 MG1655 included in `references/` folder
- **Strand-specific coverage** — forward and reverse read coverage tracked separately at each position
- **All outputs in CSV format** — suitable for publication tables and downstream analysis in R or Excel

### Usage

```bash
# Standard recombinant vector (with full contamination analysis)
python aav_genome_qc.py \
    --fastq-dir /path/to/barcode03/ \
    --cis-plasmid pAAV2ST-Pcsk9sg1.fasta \
    --rep-cap ref_rep_cap.fasta \
    --name pAAV2ST-Pcsk9sg1 \
    --output-dir ./barcode03_results/

# With custom helper and host references
python aav_genome_qc.py \
    --fastq-dir /path/to/barcode03/ \
    --cis-plasmid pAAV2ST-Pcsk9sg1.fasta \
    --rep-cap ref_rep_cap.fasta \
    --helper ref_helper.fasta \
    --host ref_host.fasta \
    --name pAAV2ST-Pcsk9sg1 \
    --output-dir ./barcode03_results/

# Skip contamination analysis
# Use when the packaged insert IS the rep-cap sequence (e.g. wild-type AAV rescue),
# or when contamination references are unavailable
python aav_genome_qc.py \
    --fastq-dir /path/to/barcode01/ \
    --cis-plasmid pSM620.fasta \
    --rep-cap ref_rep_cap.fasta \
    --no-contamination \
    --name pSM620_virus \
    --output-dir ./barcode01_results/
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--fastq-dir` | required | Directory containing fastq.gz files for one sample |
| `--cis-plasmid` | required | Full cis plasmid FASTA (ITRs auto-detected) |
| `--rep-cap` | required | Rep-Cap plasmid FASTA |
| `--helper` | bundled pHelper | Helper plasmid FASTA |
| `--host` | bundled E. coli K-12 | Host genome FASTA |
| `--serotype` | AAV2 | AAV serotype for ITR detection |
| `--name` | sample | Sample name for output files |
| `--output-dir` | . | Output directory (created if absent) |
| `--itr-tol` | 100 | ITR boundary tolerance in bp |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--no-contamination` | False | Skip contamination analysis |
| `--min-base-qual` | 10 | Minimum per-base Phred quality score (bases below threshold excluded from pileup, default 20) |
| `--validate-tsv` | None | wf-aav-qc per-read TSV for validation against existing results |

### Outputs

| File | Description |
|------|-------------|
| `{name}_detected_itrs.bed` | Auto-detected ITR coordinates (BED format, 0-based half-open) |
| `{name}_contamination.csv` | Read counts per reference (transgene/helper/rep-cap/host/backbone/unmapped) |
| `{name}_genome_types.csv` | Structural classification counts and percentages |
| `{name}_pileup.csv` | Per-position base counts, strand coverage, and deviation (ITR-to-ITR region only) |
| `{name}_itr_coverage.csv` | ITR coverage depth and strand balance statistics |
| `{name}_truncations.csv` | Per-read alignment coordinates and genome type assignments |
| `{name}_summary.csv` | Single-row QC summary for multi-sample comparison |
| `{name}_contamination.png` | Read assignment bar chart and pie chart |
| `{name}_structures.png` | Genome type and subtype breakdown |
| `{name}_coverage.png` | Coverage depth, base composition, and deviation from reference |
| `{name}_strand_coverage.png` | Forward and reverse strand coverage across the insert |
| `{name}_truncations.png` | Read start/end position distributions (truncation hotspot map) |
| `{name}_itr_detail.png` | ITR coverage detail with strand breakdown |

### Genome type classifications

| Type | Subtype | Description |
|------|---------|-------------|
| Full ssAAV | Full ssAAV | Complete ITR-to-ITR single-stranded genome |
| Partial ssAAV | 5' ICG | Incomplete genome, truncated at 3' end |
| Partial ssAAV | 3' ICG | Incomplete genome, truncated at 5' end |
| Partial ssAAV | Partial ICG - no ITRs | Internal fragment, neither ITR present |
| Partial scAAV | SBG 5'/3' symmetric/asymmetric | Snapback genome variants |
| Full scAAV | Full scAAV | Self-complementary genome |
| Reverse packaging | Reverse packaging | Maps to plasmid backbone only |
| Unmapped | Unmapped | No alignment to transgene reference |

### Notes

**Flip/flop analysis** is not performed on packaged viral genomes. The ITR hairpin structure at the ends of encapsidated genomes causes insufficient and unreliable coverage at inner loop positions with ONT sequencing. Flip/flop analysis should be performed on the plasmid stock using `aav_plasmid_workflow.py`.

**Contamination in AAV preparations** reflects encapsidated DNA rather than surface contamination. All DNA detected after DNase I treatment is physically packaged inside capsids. The helper, rep-cap, and host fractions represent low-frequency non-specific packaging events.

**Wild-type AAV preparations** where the packaged insert contains rep-cap sequence should be analyzed with `--no-contamination` to avoid misclassification of packaged genomes as rep-cap contamination.

**The pAAV-MCS deletion warning** fires when the pipeline detects an inner loop of 24 bp (vs the canonical 31 bp) caused by an 11-bp deletion (AAAGCCCGGGC) within the ITR palindrome. This deletion is present in pAAV-MCS (Agilent) and all derivative vectors. Analysis continues normally but the user is alerted.

### Run commands used in the paper

```bash
# pAAV2ST-Pcsk9sg1 packaged genomes (barcode03)
python aav_genome_qc.py \
    --fastq-dir /path/to/fastq_pass/barcode03/ \
    --cis-plasmid pAAV2ST-Pcsk9sg1.fasta \
    --rep-cap ref_rep_cap.fasta \
    --helper ref_helper.fasta \
    --host ref_host.fasta \
    --name pAAV2ST-Pcsk9sg1 \
    --output-dir ./barcode03_results/

# pSM620 packaged genomes (barcode01, wild-type AAV2 rescue)
python aav_genome_qc.py \
    --fastq-dir /path/to/fastq_pass/barcode01/ \
    --cis-plasmid pSM620_ref_rotated.fasta \
    --rep-cap ref_rep_cap.fasta \
    --no-contamination \
    --name pSM620_virus \
    --output-dir ./barcode01_results/

# pBR-AAV2(-) packaged genomes (barcode02, wild-type AAV2 rescue)
python aav_genome_qc.py \
    --fastq-dir /path/to/fastq_pass/barcode02/ \
    --cis-plasmid pAAV2-AAV2-.fasta \
    --rep-cap ref_rep_cap.fasta \
    --no-contamination \
    --name pBR-AAV2 \
    --output-dir ./barcode02_results/
```

---

## 3. PacBio Pileup (`pacbio_pileup.py`)

Analyzes PacBio CLR sequencing data from restriction-digested AAV plasmid fragments. Uses mappy (map-pb preset) for alignment to a linear fragment reference. Produces the same pileup CSV and figure formats as the plasmid workflow.

Note: this script was used for exploratory analysis of the pSM620.7 restriction fragment but the data were not included in the final paper due to homopolymer artifacts at the ITR boundaries inherent to PacBio CLR sequencing of GC-rich palindromic sequences.

### Usage

```bash
python pacbio_pileup.py \
    --bam subreads.bam \
    --fasta fragment_reference.fasta \
    --left-itr 42 185 \
    --right-itr 4583 4724 \
    --name pSM620_fragment \
    --output-dir ./pacbio_results/
```

---

## Data availability

Raw sequencing data are deposited in NCBI SRA under BioProject [PRJNA######]. Processed results and reference files are available at Zenodo (DOI: 10.5281/zenodo.#######).

---

## Citation

If you use these scripts in your work, please cite:

> Jaskula-Ranga V, Li S, Pandya A, Bunz F. *Stable propagation of inverted terminal repeats in adeno-associated virus plasmid vectors.* bioRxiv (2026).

---

## License

MIT License. See LICENSE file for details.
