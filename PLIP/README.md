# PLIP Analysis

Scripts used for protein-ligand interaction analysis after molecular docking.

## Tools

* Python
* AutoDock Vina
* PLIP
* PDB/PDBQT files

## Workflow

The analysis is performed in two steps:

```text
Docking PDBQT files
        ↓
mapk3_plip_pipeline.py
        ↓
Best docking pose
        ↓
Protein-ligand complex PDB
        ↓
run_plip_all.py
        ↓
PLIP interaction analysis
        ↓
Interaction and ligand property tables
```

## Scripts

### `mapk3_plip_pipeline.py`

This script prepares docking results for PLIP.

It:

* reads receptor and ligand PDBQT files
* extracts AutoDock Vina docking scores
* selects the lowest-energy docking pose
* converts the selected pose to PDB format
* creates protein-ligand complex structures
* optionally runs PLIP
* creates a docking pose summary table

Example:

```bash
python3 mapk3_plip_pipeline.py \
    --receptor mapk3.pdbqt \
    --ligands ligand1.pdbqt ligand2.pdbqt \
    --output MAPK3_PLIP_results \
    --run-plip
```

### `run_plip_all.py`

This script runs PLIP on the generated protein-ligand complexes and extracts information from the PLIP XML reports.

It generates:

* interaction details
* interaction summary
* ligand properties

Example:

```bash
python3 run_plip_all.py \
    --input MAPK3_PLIP_results
```

## Output

Example output files include:

```text
docking_pose_summary.tsv
PLIP_interactions.tsv
PLIP_summary.tsv
PLIP_ligand_properties.tsv
```

Individual PLIP reports are stored separately for each compound.

Docking structures, PDB/PDBQT files and generated PLIP results are not included in this repository.

These scripts were written and used during my bioinformatics project work.

