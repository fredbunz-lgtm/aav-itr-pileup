# Bundled Reference Sequences

This folder contains default reference sequences used by `aav_genome_qc.py`
for contamination analysis. These files must remain in a `references/`
directory alongside the script.

## ref_helper_pHelper.fasta

**pHelper plasmid** (NCBI accession AF369965.1, 11,635 bp)

Contains the adenovirus 5 helper functions required for AAV packaging:
E2A (DNA-binding protein), E4 ORFs, and VA RNA genes. This plasmid is
standard across virtually all helper-free AAV production systems
(Agilent/Stratagene and equivalent).

Source: https://www.ncbi.nlm.nih.gov/nuccore/AF369965.1

## ref_host_ecoli_K12.fasta

**Escherichia coli K-12 substr. MG1655 complete genome**
(NCBI RefSeq accession NC_000913.3, 4,641,652 bp)

Used for detection of encapsidated host DNA. Suitable for detection of
contamination from any common laboratory *E. coli* strain (DH5α, DH10B,
XL10-Gold, etc.) as these strains are all derived from K-12 and share
>99% sequence identity across coding regions.

Source: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3

## Overriding defaults

Both references can be replaced with custom sequences using the
`--helper` and `--host` arguments:

```bash
python aav_genome_qc.py \
    --fastq-dir /path/to/fastq/ \
    --cis-plasmid my_plasmid.fasta \
    --rep-cap my_rep_cap.fasta \
    --helper /path/to/my_helper.fasta \
    --host /path/to/my_host.fasta \
    --name my_sample
```

The rep-cap plasmid reference is not bundled as it varies by AAV serotype
and construct. It must always be provided by the user via `--rep-cap`.
