#!/usr/bin/env python3

# Generate a heatmap from the pairwise identity
# matrix calculated from EMBOSS Needle results.

from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path("results/protein_identity_matrix.tsv")
)

OUTPUT_FILE = Path("results/protein_identity_heatmap.png")


if not INPUT_FILE.is_file():
    raise SystemExit(f"Error: identity matrix not found: {INPUT_FILE}")


# Read identity matrix
df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    index_col=0
)

if df.empty:
    raise SystemExit("Error: identity matrix is empty.")


# Create heatmap
fig, ax = plt.subplots(figsize=(9, 7))

image = ax.imshow(
    df.values,
    cmap="viridis",
    vmin=0,
    vmax=100
)


# Axis labels
ax.set_xticks(range(len(df.columns)))
ax.set_yticks(range(len(df.index)))

ax.set_xticklabels(
    df.columns,
    rotation=45,
    ha="right",
    fontsize=10
)

ax.set_yticklabels(
    df.index,
    fontsize=10
)


# Add values to cells
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
cbar.set_label("Global Pairwise Identity (%)")


ax.set_title(
    "Pairwise Protein Identity\n"
    "EMBOSS Needle Global Alignment",
    fontsize=14,
    pad=15
)

ax.set_xlabel("Species")
ax.set_ylabel("Species")


plt.tight_layout()

plt.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print()
print("=" * 70)
print("NEEDLE IDENTITY HEATMAP GENERATED")
print("=" * 70)
print()
print(f"Input : {INPUT_FILE}")
print(f"Output: {OUTPUT_FILE}")
