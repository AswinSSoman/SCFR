#!/usr/bin/env python3
import sys

def load_sequences(fasta_file):
    sequences = set()
    header = None
    seq_lines = []

    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header and seq_lines:
                    seq = "".join(seq_lines).upper()  # normalize case
                    sequences.add(seq)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            seq = "".join(seq_lines).upper()
            sequences.add(seq)

    return sequences


def main():
    if len(sys.argv) != 3:
        print("Usage: subtract_fasta.py A.fasta B.fasta > filtered_B.fasta")
        sys.exit(1)

    file_a = sys.argv[1]
    file_b = sys.argv[2]

    seqs_a = load_sequences(file_a)
    print(f"Loaded {len(seqs_a)} sequences from {file_a}", file=sys.stderr)

    kept_count = 0
    header = None
    seq_lines = []

    with open(file_b) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header and seq_lines:
                    seq_norm = "".join(seq_lines).upper()
                    if seq_norm not in seqs_a:
                        print(header)
                        print("\n".join(seq_lines))
                        kept_count += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # Final record
        if header and seq_lines:
            seq_norm = "".join(seq_lines).upper()
            if seq_norm not in seqs_a:
                print(header)
                print("\n".join(seq_lines))
                kept_count += 1

    print(f"Kept {kept_count} sequences from {file_b}", file=sys.stderr)

if __name__ == "__main__":
    main()

