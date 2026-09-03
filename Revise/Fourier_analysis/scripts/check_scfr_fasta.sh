#!/bin/bash

# Process and store standard input into memory
INPUT=$(cat)

# Exit early if input is empty
if [ -z "$INPUT" ]; then
    echo "Error: No matching BLAST hits passed to check_scfr_fasta.sh" >&2
    exit 1
fi

# Create temporary file to store all converted BED regions
BED_TMP=$(mktemp)
trap 'rm -f "$BED_TMP"' EXIT

# Read through every line from the grep pipeline
echo "$INPUT" | while read -r LINE; do
    [ -z "$LINE" ] && continue

    # Extract fields from BLAST output ($1=query_id, $5=q.start, $6=q.end)
    Q_ID=$(echo "$LINE" | awk '{print $1}')
    Q_START=$(echo "$LINE" | awk '{print $5}')
    Q_END=$(echo "$LINE" | awk '{print $6}')

    # Parse SCFR region
    SCFR_REGION=$(echo "$Q_ID" | sed 's/^:://g' | tr ":-" "\t" | awk '{print $1":"$2"-"$3}')
    SCFR_CHROM=$(echo "$Q_ID" | sed 's/^:://g' | tr ":-" "\t" | awk '{print $1}')
    SCFR_BASE_START=$(echo "$Q_ID" | sed 's/^:://g' | tr ":-" "\t" | awk '{print $2}')

    # Calculate hit coordinates and strand orientation
    if [ "$Q_START" -le "$Q_END" ]; then
        STRAND="+"
        FRAME="F"
        HIT_START=$((SCFR_BASE_START + Q_START - 1))
        HIT_END=$((SCFR_BASE_START + Q_END))
    else
        STRAND="-"
        FRAME="R"
        HIT_START=$((SCFR_BASE_START + Q_END - 1))
        HIT_END=$((SCFR_BASE_START + Q_START))
    fi

    HIT_REGION="${SCFR_CHROM}:${HIT_START}-${HIT_END}"

    # Print record metadata summary
    echo "=================================================="
    echo "Source SCFR: $SCFR_REGION"
    echo "Hit within SCFR: $HIT_REGION"
    echo "Query Strand: $STRAND (Using Frame: $FRAME)"
    echo "=================================================="

    # Run bedtools and transeq for this individual region
    echo -e "${SCFR_CHROM}\t${HIT_START}\t${HIT_END}" \
      | bedtools getfasta -fi GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed - \
      | transeq -sequence stdin -frame "$FRAME" --auto --stdout \
      | ~/aswin/programmes/myfasta -comb

    echo "" # Add blank line spacing between entries
done
