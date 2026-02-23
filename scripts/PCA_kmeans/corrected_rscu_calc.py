import os
import argparse
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import numpy as np

# ---------------- CODON TABLE ---------------- #

codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','CAT':'H','CAC':'H',
    'CAA':'Q','CAG':'Q','AAT':'N','AAC':'N',
    'AAA':'K','AAG':'K','GAT':'D','GAC':'D',
    'GAA':'E','GAG':'E','TGT':'C','TGC':'C',
    'TGG':'W','CGT':'R','CGC':'R','CGA':'R',
    'CGG':'R','AGT':'S','AGC':'S','AGA':'R',
    'AGG':'R','GGT':'G','GGC':'G','GGA':'G',
    'GGG':'G'
}

amino_to_codons = defaultdict(list)
for codon, aa in codon_table.items():
    amino_to_codons[aa].append(codon)

# ---------------- FUNCTIONS ---------------- #

def codon_usage(seq):
    counts = defaultdict(int)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in codon_table:
            counts[codon] += 1
    return counts

def calculate_rscu(counts):
    rscu = {}
    for aa, codons in amino_to_codons.items():
        total = sum(counts[c] for c in codons)
        if total == 0:
            for c in codons:
                rscu[c] = 0.0
        else:
            exp = total / len(codons)
            for c in codons:
                rscu[c] = counts[c] / exp if exp > 0 else 0.0
    return rscu

def calculate_gc3(seq):
    third = seq[2::3]
    if not third:
        return 0.0
    return sum(1 for b in third if b in "GC") / len(third)

def gc3_correct_rscu(rscu_df, gc3_series):
    corrected = pd.DataFrame(index=rscu_df.index, columns=rscu_df.columns)
    for codon in rscu_df.columns:
        y = rscu_df[codon].values
        x = gc3_series.values
        if np.std(y) == 0:
            corrected[codon] = 0.0
            continue
        slope, intercept = np.polyfit(x, y, 1)
        predicted = slope * x + intercept
        corrected[codon] = y - predicted
    return corrected

# ---------------- MAIN ---------------- #

def process_directory(fasta_dir, output_dir):

    for file in os.listdir(fasta_dir):
        if not file.lower().endswith((".fa", ".fasta")):
            continue

        base = os.path.splitext(file)[0]

        gc3_rows = {}
        rscu_rows = {}

        path = os.path.join(fasta_dir, file)

        for record in SeqIO.parse(path, "fasta"):
            seq_id = record.id
            seq = str(record.seq).upper().replace("-", "")

            if len(seq) < 3 or len(seq) % 3 != 0:
                continue

            counts = codon_usage(seq)
            rscu = calculate_rscu(counts)
            gc3 = calculate_gc3(seq)

            gc3_rows[seq_id] = gc3
            rscu_rows[seq_id] = rscu

        if not rscu_rows:
            continue

        # ---- DATAFRAMES ---- #
        gc3_df = pd.DataFrame.from_dict(gc3_rows, orient="index", columns=["GC3"])
        rscu_raw_df = pd.DataFrame.from_dict(rscu_rows, orient="index")

        rscu_corr_df = gc3_correct_rscu(
            rscu_raw_df,
            gc3_df.loc[rscu_raw_df.index, "GC3"]
        )

        # ---- WRITE FLAT IN OUTPUT_DIR ---- #
        gc3_df.to_csv(os.path.join(output_dir, f"gc3.tsv"), sep="\t")
        rscu_raw_df.to_csv(os.path.join(output_dir, f"rscu_raw_matrix.tsv"), sep="\t")
        rscu_corr_df.to_csv(os.path.join(output_dir, f"rscu_gc3_corrected_matrix.tsv"), sep="\t")


# ---------------- ENTRY ---------------- #

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    process_directory(args.fasta_dir, args.output_dir)
