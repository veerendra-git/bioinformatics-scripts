#!/usr/bin/env python3

# Generate a heatmap from an MSA-based pairwise
# percentage identity matrix.

from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path("results/MSA_identity_matrix.tsv")
)

OUTPUT_FILE = Path("results/MSA_identity_heatmap.png")

# Create the output directory if it does not exist
OUTPUT_FILE.parent.mkdir(exist_ok=True)


# Check that the input matrix exists
if not INPUT_FILE.is_file():
    raise SystemExit(f"Error: identity matrix not found: {INPUT_FILE}")


# Read the identity matrix
df = pd.read_csv(INPUT_FILE, sep="\t")

if "Species" not in df.columns:
    raise SystemExit("Error: the matrix must contain a 'Species' column.")

df = df.set_index("Species")

if df.empty:
    raise SystemExit("Error: identity matrix is empty.")


# Create heatmap
fig, ax = plt.subplots(figsize=(9, 8))

image = ax.imshow(
    df.values,
    vmin=0,
    vmax=100,
    aspect="equal"
)


# Axis labels
ax.set_xticks(range(len(df.columns)))
ax.set_yticks(range(len(df.index)))

ax.set_xticklabels(
    df.columns,
    rotation=45,
    ha="right"
)

ax.set_yticklabels(df.index)

ax.set_xlabel("Species")
ax.set_ylabel("Species")

ax.set_title(
    "Pairwise Protein Identity\n"
    "Multiple Sequence Alignment"
)


# Add percentage values to the cells
for i in range(len(df.index)):
    for j in range(len(df.columns)):

        value = df.iloc[i, j]

        ax.text(
            j,
            i,
            f"{value:.1f}%",
            ha="center",
            va="center",
            fontsize=9
        )


# Colorbar
cbar = fig.colorbar(image, ax=ax)
cbar.set_label("Pairwise Identity (%)")


fig.tight_layout()

fig.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print()
print("=" * 70)
print("MSA IDENTITY HEATMAP GENERATED")
print("=" * 70)
print()
print(f"Input : {INPUT_FILE}")
print(f"Output: {OUTPUT_FILE}")
