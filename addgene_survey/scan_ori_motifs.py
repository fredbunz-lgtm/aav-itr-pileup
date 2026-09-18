from Bio import SeqIO
from Bio.Seq import Seq
import sys

MOTIFS = {
    "RNAII_short": "CTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAA",
    "pUC_pMB1_core": "TTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTA",
    "pUC_highcopy": "GCGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAA",
    "colE1_RNAII": "ATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAG",
    "ori_adjacent_pBR322": "TCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTC",
}

def find_all_circular(seq, motif):
    seq = seq.upper()
    motif = motif.upper()
    n = len(seq)
    doubled = seq + seq
    hits = []
    start = 0
    while True:
        i = doubled.find(motif, start)
        if i == -1:
            break
        if i < n:
            hits.append(i)
        start = i + 1
    return hits


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 scan_ori_motifs.py plasmid.fasta")
        sys.exit(1)

    for rec in SeqIO.parse(sys.argv[1], "fasta"):
        seq = str(rec.seq).upper()
        n = len(seq)
        print(f"Record: {rec.id}, length={n}")
        found = False

        for label, motif in MOTIFS.items():
            rc = str(Seq(motif).reverse_complement())
            for strand, query in [("+", motif), ("-", rc)]:
                for pos in find_all_circular(seq, query):
                    found = True
                    end = (pos + len(query)) % n
                    print(f"{label}\tstrand={strand}\tstart_0based={pos}\tend_0based_exclusive={end}\tlength={len(query)}")

        if not found:
            print("No exact ori motif hits found.")


if __name__ == "__main__":
    main()
