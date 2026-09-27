#!/bin/bash

# Check FASTA protein sequences before running the alignments

set -euo pipefail

SEQ_DIR="${1:-sequences}"

if [ ! -d "$SEQ_DIR" ]; then
    echo "Error: sequence directory not found: $SEQ_DIR"
    echo "Usage: bash 01_verify_sequences.sh <sequence_directory>"
    exit 1
fi

echo "Checking FASTA files in: $SEQ_DIR"
echo

found_file=false

for file in "$SEQ_DIR"/*.fasta "$SEQ_DIR"/*.fa; do

    if [ ! -f "$file" ]; then
        continue
    fi

    found_file=true
    name=$(basename "$file")

    # Count amino acids excluding the FASTA header
    length=$(grep -v "^>" "$file" | tr -d '\n\r ' | wc -c)

    echo "$name : $length amino acids"

    # Check for characters other than standard amino acid symbols
    invalid=$(grep -v "^>" "$file" | tr -d '[:space:]' | \
        grep -o '[^ACDEFGHIKLMNPQRSTVWY]' || true)

    if [ -n "$invalid" ]; then
        echo "  Warning: unusual amino acid character(s) found"
    else
        echo "  Sequence check: OK"
    fi

    echo
done

if [ "$found_file" = false ]; then
    echo "No FASTA files found in $SEQ_DIR"
    exit 1
fi

echo "Sequence verification completed."
