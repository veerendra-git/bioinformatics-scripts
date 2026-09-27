# Pairwise Protein Sequence Analysis

This folder contains scripts used for protein sequence verification, pairwise alignment, identity calculation and heatmap generation.

## Tools

* Bash
* EMBOSS Needle
* Clustal Omega
* Python 3
* Pandas
* Matplotlib

## Scripts

| Script                          | Purpose                                                |
| ------------------------------- | ------------------------------------------------------ |
| `01_verify_sequences.sh`        | Checks FASTA files and sequence quality                |
| `02_pairwise_alignment.sh`      | Performs pairwise global alignment using EMBOSS Needle |
| `03_identity_matrix.py`         | Creates an identity matrix from Needle results         |
| `04_msa_identity_matrix.py`     | Calculates identity from a multiple sequence alignment |
| `05_identity_heatmap.py`        | Generates a heatmap from the MSA identity matrix       |
| `06_needle_identity_heatmap.py` | Generates a heatmap from the Needle identity matrix    |

## Workflow

### Pairwise alignment

```text
FASTA sequences
      ↓
01_verify_sequences.sh
      ↓
02_pairwise_alignment.sh
      ↓
03_identity_matrix.py
      ↓
06_needle_identity_heatmap.py
```

### Multiple sequence alignment

```text
FASTA sequences
      ↓
Clustal Omega
      ↓
04_msa_identity_matrix.py
      ↓
05_identity_heatmap.py
```

## Input

Protein FASTA files should be placed in a `sequences/` directory.

The scripts can work with different numbers of protein sequences. Pairwise comparisons are generated automatically based on the available FASTA files.

## Main Outputs

```text
alignments/
results/
├── protein_identity_matrix.tsv
├── protein_identity_heatmap.png
├── MSA_identity_matrix.tsv
└── MSA_identity_heatmap.png
```

Input sequences and generated results are not included in this repository.

## Notes

The scripts were developed and used during bioinformatics coursework and project work. They can be modified for other protein sequence comparison studies.
