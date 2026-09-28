# MD Analysis

Scripts used for molecular dynamics trajectory analysis and free energy analysis.

## Tools

* R
* Bio3D
* MASS
* fields
* Plotly

## Scripts

### `PDB_TO_DCD.py`

Combines multiple PDB trajectory frames into:

* a merged PDB file
* a DCD trajectory file

It uses **MDAnalysis** for reading the PDB files and writing the DCD trajectory.

### `FEL_PCA_analysis.R`

Performs analysis of an MD trajectory using a reference PDB structure and DCD trajectory.

The script includes:

* RMSD analysis
* RMSF analysis
* Principal Component Analysis (PCA)
* PCA clustering
* PC1 and PC2 displacement analysis
* Dynamic Cross-Correlation Matrix (DCCM)
* RMSD-based free energy profile
* 2D PCA free energy landscape
* Zoomed 2D free energy landscape
* 3D free energy landscape

## Input

The main R script uses:

```text
Ref.pdb
Output.dcd
```

The input files can also be provided as command-line arguments:

```bash
Rscript FEL_PCA_analysis.R Ref.pdb Output.dcd MD_Results
```

## Output

Results are saved in the specified output directory.

Example outputs include:

```text
RMSD.png
RMSF.png
PCA.png
Clustered_PCA.png
DCCM.png
FEL_RMSD.png
FEL_2D.png
FEL_2D_ZoomedIn.png
3D_FEL_Funnel.html
```

The PCA displacement values are also saved as a CSV file.

MD trajectory and generated result files are not included in this repository.

These scripts were written and used during my bioinformatics practical/project work.
