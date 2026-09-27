#!/usr/bin/env python3

# Extract percentage identity from EMBOSS Needle results
# and create a pairwise identity matrix.

from pathlib import Path
import re
import csv
import sys


ALIGN_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("alignments")
RESULT_DIR = Path("results")

RESULT_DIR.mkdir(exist_ok=True)


if not ALIGN_DIR.is_dir():
    raise SystemExit(f"Error: alignment directory not found: {ALIGN_DIR}")


alignment_files = sorted(ALIGN_DIR.glob("*.needle"))

if not alignment_files:
    raise SystemExit(f"Error: no .needle files found in {ALIGN_DIR}")


# Get the species names from the alignment filenames
species_set = set()

for file in alignment_files:

    pair = file.stem.split("_vs_")

    if len(pair) != 2:
        raise SystemExit(
            f"Error: could not identify species pair in {file.name}"
        )

    species_set.add(pair[0])
    species_set.add(pair[1])


species = sorted(species_set)

if len(species) < 2:
    raise SystemExit("Error: at least two species are required.")


# Number of unique pairwise alignments expected
expected_files = len(species) * (len(species) - 1) // 2

if len(alignment_files) != expected_files:
    raise SystemExit(
        f"Error: expected {expected_files} alignment files for "
        f"{len(species)} species, but found {len(alignment_files)}."
    )


# Create an empty matrix
matrix = {
    sp1: {sp2: None for sp2 in species}
    for sp1 in species
}

# A sequence compared with itself has 100% identity
for sp in species:
    matrix[sp][sp] = 100.0


# EMBOSS Needle identity line:
# # Identity:     318/394 (80.7%)
pattern = re.compile(
    r"# Identity:\s+\d+/\d+\s+\((\d+(?:\.\d+)?)%\)"
)


for file in alignment_files:

    pair = file.stem.split("_vs_")
    sp1, sp2 = pair

    text = file.read_text()

    match = pattern.search(text)

    if not match:
        raise SystemExit(
            f"Error: identity value not found in {file.name}"
        )

    identity = float(match.group(1))

    matrix[sp1][sp2] = identity
    matrix[sp2][sp1] = identity


# Check that the complete matrix was filled
for sp1 in species:
    for sp2 in species:

        if matrix[sp1][sp2] is None:
            raise SystemExit(
                f"Error: missing identity value for "
                f"{sp1} vs {sp2}"
            )


# Save the matrix
output_file = RESULT_DIR / "protein_identity_matrix.tsv"

with output_file.open("w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["Species"] + species)

    for sp1 in species:

        row = [sp1]

        for sp2 in species:
            row.append(f"{matrix[sp1][sp2]:.1f}")

        writer.writerow(row)


# Print the matrix
print()
print("=" * 90)
print("PAIRWISE PERCENT IDENTITY MATRIX")
print("=" * 90)
print()

header = f"{'Species':<22}"

for sp in species:
    header += f"{sp:>18}"

print(header)
print("-" * 148)

for sp1 in species:

    row = f"{sp1:<22}"

    for sp2 in species:
        row += f"{matrix[sp1][sp2]:>18.1f}"

    print(row)

print()
print(f"Number of species: {len(species)}")
print(f"Matrix size: {len(species)} x {len(species)}")
print(f"Matrix saved to: {output_file}")
