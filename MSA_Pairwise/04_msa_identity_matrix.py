#!/usr/bin/env python3

# Calculate pairwise percentage identity from a
# multiple sequence alignment.
#
# Columns containing a gap in either sequence are excluded.

from pathlib import Path
import csv
import sys


MSA_FILE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("alignments/alignment.fasta")
RESULT_DIR = Path("results")

RESULT_DIR.mkdir(exist_ok=True)


if not MSA_FILE.is_file():
    raise SystemExit(f"Error: MSA file not found: {MSA_FILE}")


# Read the MSA FASTA file
sequences = {}
current_name = None

with MSA_FILE.open() as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            current_name = line[1:].split("|")[0].strip()

            if current_name in sequences:
                raise SystemExit(
                    f"Error: duplicate sequence name found: {current_name}"
                )

            sequences[current_name] = ""

        else:

            if current_name is None:
                raise SystemExit(
                    "Error: sequence data found before a FASTA header."
                )

            sequences[current_name] += line


if len(sequences) < 2:
    raise SystemExit(
        "Error: the MSA must contain at least two sequences."
    )


# Keep the sequence order from the MSA file
species = list(sequences.keys())


# Check that all aligned sequences have the same length
lengths = {
    sp: len(sequences[sp])
    for sp in species
}

if len(set(lengths.values())) != 1:
    raise SystemExit(
        f"Error: MSA sequences have different lengths: {lengths}"
    )


alignment_length = len(sequences[species[0]])


# Create the identity matrix
matrix = {
    sp1: {sp2: None for sp2 in species}
    for sp1 in species
}


for sp1 in species:

    matrix[sp1][sp1] = 100.0

    for sp2 in species:

        if sp1 == sp2:
            continue

        seq1 = sequences[sp1]
        seq2 = sequences[sp2]

        identical = 0
        comparable = 0

        for aa1, aa2 in zip(seq1, seq2):

            # Ignore columns where either sequence has a gap
            if aa1 == "-" or aa2 == "-":
                continue

            comparable += 1

            if aa1 == aa2:
                identical += 1

        if comparable == 0:
            identity = 0.0
        else:
            identity = (identical / comparable) * 100

        matrix[sp1][sp2] = identity


# Save the matrix
output_file = RESULT_DIR / "MSA_identity_matrix.tsv"

with output_file.open("w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["Species"] + species)

    for sp1 in species:

        row = [sp1]

        for sp2 in species:
            row.append(f"{matrix[sp1][sp2]:.2f}")

        writer.writerow(row)


# Print the results
print()
print("=" * 100)
print("MSA-BASED PAIRWISE PERCENT IDENTITY MATRIX")
print("=" * 100)
print()

print(f"Number of sequences: {len(species)}")
print(f"Alignment length: {alignment_length} columns")
print("Gap-containing columns are excluded from each pairwise calculation.")
print()

header = f"{'Species':<22}"

for sp in species:
    header += f"{sp:>18}"

print(header)
print("-" * 148)

for sp1 in species:

    row = f"{sp1:<22}"

    for sp2 in species:
        row += f"{matrix[sp1][sp2]:>18.2f}"

    print(row)

print()
print(f"Matrix size: {len(species)} x {len(species)}")
print(f"Matrix saved to: {output_file}")
