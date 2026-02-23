#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Define the set of stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}

def find_stop_codon_free_regions(seq, frame, seq_id, strand="+"):
    """
    Identifies all stop-codon-free regions in a given reading frame.

    Parameters:
    - seq: Biopython Seq object
    - frame: Reading frame (0, 1, or 2)
    - seq_id: FASTA record ID
    - strand: '+' for forward, '-' for reverse

    Returns:
    - List of tuples: (seq_id, start, end, frame_number)
    """
    results = []
    start = None
    seq_len = len(seq)

    for i in range(frame, seq_len - 2, 3):
        codon = str(seq[i:i+3]).upper()
        if codon in STOP_CODONS:
            # End current region if a stop codon is found
            if start is not None:
                end = i
                if end > start:
                    if strand == "+":
                        results.append((seq_id, start, end, frame + 1))
                    else:
                        # Adjust coordinates for reverse strand
                        real_start = seq_len - end
                        real_end = seq_len - start
                        results.append((seq_id, real_start, real_end, -(frame + 1)))
                start = None
        else:
            # Start a new region if not already started
            if start is None:
                start = i

    # Handle tail (remaining sequence after last codon)
    if start is not None and seq_len - start >= 3:
        if strand == "+":
            results.append((seq_id, start, seq_len, frame + 1))
        else:
            real_start = 0
            real_end = seq_len - start
            results.append((seq_id, real_start, real_end, -(frame + 1)))

    return results

def find_gap_regions(seq, seq_id):
    """
    Identifies regions with 'N' or '-' characters in the sequence.

    Parameters:
    - seq: Biopython Seq object
    - seq_id: FASTA record ID

    Returns:
    - List of tuples: (seq_id, start, end) for each ambiguous region
    """
    gap_regions = []
    in_gap = False
    start = None

    for i, base in enumerate(seq.upper()):
        if base in {"N", "-"}:
            if not in_gap:
                in_gap = True
                start = i
        else:
            if in_gap:
                in_gap = False
                gap_regions.append((seq_id, start, i))

    # Handle end of sequence gap
    if in_gap:
        gap_regions.append((seq_id, start, len(seq)))

    return gap_regions

def process_fasta(fasta_path, gap_output_path):
    with open(gap_output_path, "w") as gap_out:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_id = record.id
            seq = record.seq
            rev_seq = seq.reverse_complement()

            # Report gap/N regions
            gap_regions = find_gap_regions(seq, seq_id)
            for chrom, start, end in gap_regions:
                gap_out.write(f"{chrom}\t{start}\t{end}\n")

            # Process forward strand in all 3 frames
            for frame in range(3):
                regions = find_stop_codon_free_regions(seq, frame, seq_id, strand="+")
                for chrom, start, end, frame_num in regions:
                    print(f"{chrom}\t{start}\t{end}\t{frame_num}")

            # Process reverse strand in all 3 frames
            for frame in range(3):
                regions = find_stop_codon_free_regions(rev_seq, frame, seq_id, strand="-")
                for chrom, start, end, frame_num in regions:
                    print(f"{chrom}\t{start}\t{end}\t{frame_num}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find stop codon–free regions and report N/gap regions from both strands of a FASTA file.")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("--gaps", default="gap_regions.bed", help="Output BED file for N or gap regions")
    args = parser.parse_args()

    process_fasta(args.fasta, args.gaps)
