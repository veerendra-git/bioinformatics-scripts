#!/bin/bash

# Generate pairwise global alignments for FASTA sequences
# Tool: EMBOSS Needle

set -euo pipefail

SEQ_DIR="${1:-sequences}"
OUT_DIR="alignments"

# Check that EMBOSS Needle is installed
if ! command -v needle >/dev/null 2>&1; then
    echo "Error: EMBOSS needle is not installed."
    echo "Install it with: sudo apt install emboss"
    exit 1
fi

if [ ! -d "$SEQ_DIR" ]; then
    echo "Error: sequence directory not found: $SEQ_DIR"
    exit 1
fi

# Find all FASTA files
mapfile -t SEQUENCES < <(
    find "$SEQ_DIR" -maxdepth 1 -type f \
    \( -name "*.fasta" -o -name "*.fa" \) | sort
)

if [ ${#SEQUENCES[@]} -lt 2 ]; then
    echo "Error: at least two FASTA files are required."
    exit 1
fi

mkdir -p "$OUT_DIR"

# Number of unique pairs = n(n-1)/2
total=$((${#SEQUENCES[@]} * (${#SEQUENCES[@]} - 1) / 2))
count=0

echo "=========================================="
echo " Pairwise Protein Alignment"
echo "=========================================="
echo
echo "Sequences found: ${#SEQUENCES[@]}"
echo "Alignments to generate: $total"
echo

for ((i=0; i<${#SEQUENCES[@]}; i++)); do

    for ((j=i+1; j<${#SEQUENCES[@]}; j++)); do

        seq1="${SEQUENCES[$i]}"
        seq2="${SEQUENCES[$j]}"

        species1=$(basename "$seq1")
        species2=$(basename "$seq2")

        species1="${species1%.*}"
        species2="${species2%.*}"

        output="$OUT_DIR/${species1}_vs_${species2}.needle"

        count=$((count + 1))

        echo "[$count/$total] $species1 vs $species2"

        needle \
            -asequence "$seq1" \
            -bsequence "$seq2" \
            -gapopen 10.0 \
            -gapextend 0.5 \
            -outfile "$output" \
            -auto

    done

done

echo
echo "=========================================="
echo " Alignment completed"
echo "=========================================="
echo
echo "Total alignments generated: $count"
echo "Output directory: $OUT_DIR/"
