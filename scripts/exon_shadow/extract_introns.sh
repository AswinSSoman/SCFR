#!/usr/bin/env bash

# ==============================
# Usage:
# ./extract_introns.sh exons.tsv genome.fa introns.fa
# ==============================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <exon_table.tsv> <genome.fa> <output_introns.fa>"
    exit 1
fi

EXONS=$1
GENOME=$2
OUTFA=$3

TMP_BED=$(mktemp)

# ------------------------------
# Create intron BED
# ------------------------------
awk 'BEGIN{OFS="\t"}
{
    chrom=$1
    start=$2
    end=$3
    transcript=$5
    strand=$6

    if (transcript == prev_transcript) {
        intron_start = prev_end + 1
        intron_end   = start - 1

        if (intron_start <= intron_end) {
            print chrom, intron_start-1, intron_end, transcript"_intron"++i, ".", strand
        }
    } else {
        i=0
    }

    prev_transcript=transcript
    prev_end=end
}' <(sort -k5,5 -k2,2n "$EXONS") > "$TMP_BED"

# ------------------------------
# Extract FASTA
# ------------------------------
bedtools getfasta \
    -fi "$GENOME" \
    -bed "$TMP_BED" \
    -s \
    -name \
    -fo "$OUTFA"

rm "$TMP_BED"

echo "✅ Intron FASTA written to: $OUTFA"

