import os
import glob
import sys

import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter


# Input directory and output files
pdb_dir = sys.argv[1] if len(sys.argv) > 1 else "Trj"
output_pdb = sys.argv[2] if len(sys.argv) > 2 else "Merged.pdb"
output_dcd = sys.argv[3] if len(sys.argv) > 3 else "Output.dcd"


# Find PDB files
pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))

if not pdb_files:
    print(f"No PDB files found in {pdb_dir}.")
    sys.exit(1)

print(f"Found {len(pdb_files)} PDB files.")


# Check that all PDB files have the same number of atoms
first_universe = mda.Universe(pdb_files[0])
n_atoms = first_universe.atoms.n_atoms

for pdb_file in pdb_files[1:]:
    universe = mda.Universe(pdb_file)

    if universe.atoms.n_atoms != n_atoms:
        print(f"Atom count mismatch: {pdb_file}")
        print(f"Expected {n_atoms} atoms, found {universe.atoms.n_atoms}.")
        sys.exit(1)


# Merge the PDB files into one file
with open(output_pdb, "w") as output:
    for pdb_file in pdb_files:

        with open(pdb_file, "r") as input_file:
            lines = input_file.readlines()

        if not lines:
            print(f"Skipping empty file: {pdb_file}")
            continue

        # Remove the final END record if present
        if lines[-1].startswith("END"):
            lines = lines[:-1]

        output.writelines(lines)
        output.write("TER\n")

        print(f"Added: {pdb_file}")

print(f"Merged PDB file created: {output_pdb}")


# Convert the PDB frames to DCD
with DCDWriter(output_dcd, n_atoms=n_atoms, overwrite=True) as writer:

    for pdb_file in pdb_files:
        universe = mda.Universe(pdb_file)
        writer.write(universe)

print(f"DCD file created: {output_dcd}")
print(f"Number of frames: {len(pdb_files)}")
