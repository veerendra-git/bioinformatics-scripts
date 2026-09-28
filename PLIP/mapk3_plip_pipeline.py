#!/usr/bin/env python3

"""
Docking pose selection and PLIP preparation pipeline.

Workflow:
1. Read receptor PDBQT
2. Read ligand docking-result PDBQT files
3. Extract Vina docking scores
4. Select the lowest-energy pose for each ligand
5. Convert receptor and selected poses to PDB
6. Create protein-ligand complexes
7. Optionally run PLIP
8. Save a docking pose summary

No docking, minimization, translation, rotation,
or coordinate modification is performed.
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path


VINA_RESULT_PATTERN = re.compile(
    r"REMARK\s+VINA RESULT:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))"
)

ATOM_RECORDS = ("ATOM", "HETATM")


def check_file(path, description):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"{description} does not exist:\n{path}"
        )

    if not path.is_file():
        raise ValueError(
            f"{description} is not a regular file:\n{path}"
        )

    return path.resolve()


def parse_docking_pdbqt(pdbqt_file):
    """Read docking poses and their Vina scores."""

    pdbqt_file = check_file(
        pdbqt_file,
        "Ligand PDBQT"
    )

    with open(
        pdbqt_file,
        "r",
        errors="replace"
    ) as fh:
        lines = [
            line.rstrip("\n")
            for line in fh
        ]

    has_models = any(
        line.startswith("MODEL")
        for line in lines
    )

    # Handle a single-pose PDBQT
    if not has_models:

        score = None

        for line in lines:

            match = VINA_RESULT_PATTERN.search(line)

            if match:
                score = float(match.group(1))
                break

        atom_lines = [
            line
            for line in lines
            if line.startswith(ATOM_RECORDS)
        ]

        if not atom_lines:
            raise ValueError(
                f"No ATOM/HETATM records found in {pdbqt_file}"
            )

        if score is None:
            raise ValueError(
                f"No VINA RESULT found in {pdbqt_file}"
            )

        return [{
            "model": 1,
            "score": score,
            "records": lines
        }]

    # Handle MODEL/ENDMDL docking poses
    models = []

    current_model = None
    current_records = []
    current_score = None

    for line in lines:

        if line.startswith("MODEL"):

            if current_model is not None:

                if current_score is not None:
                    models.append({
                        "model": current_model,
                        "score": current_score,
                        "records": current_records
                    })

            parts = line.split()

            try:
                current_model = int(parts[1])
            except (IndexError, ValueError):
                current_model = len(models) + 1

            current_records = [line]
            current_score = None

            continue

        if line.startswith("ENDMDL"):

            if current_model is None:
                continue

            current_records.append(line)

            if current_score is None:
                raise ValueError(
                    f"MODEL {current_model} in {pdbqt_file} "
                    f"has no VINA RESULT score."
                )

            models.append({
                "model": current_model,
                "score": current_score,
                "records": current_records
            })

            current_model = None
            current_records = []
            current_score = None

            continue

        if current_model is not None:

            current_records.append(line)

            match = VINA_RESULT_PATTERN.search(line)

            if match:
                current_score = float(match.group(1))

    if current_model is not None:

        if current_score is None:
            raise ValueError(
                f"MODEL {current_model} in {pdbqt_file} "
                f"was not closed correctly."
            )

        models.append({
            "model": current_model,
            "score": current_score,
            "records": current_records
        })

    if not models:
        raise ValueError(
            f"No valid docking poses found in {pdbqt_file}"
        )

    return models


def select_best_pose(poses):
    """Select the pose with the lowest Vina score."""

    if not poses:
        raise ValueError(
            "No docking poses available."
        )

    return min(
        poses,
        key=lambda pose: pose["score"]
    )


def infer_element(atom_name):
    """Convert PDBQT atom types to chemical element symbols."""

    name = atom_name.strip().upper()

    atom_type_map = {
        "C": "C",
        "A": "C",
        "N": "N",
        "NA": "N",
        "NS": "N",
        "O": "O",
        "OA": "O",
        "S": "S",
        "SA": "S",
        "P": "P",
        "H": "H",
        "HD": "H",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
        "B": "B",
        "MG": "Mg",
        "CA": "Ca",
        "MN": "Mn",
        "FE": "Fe",
        "ZN": "Zn",
        "CU": "Cu",
        "CO": "Co",
        "NI": "Ni",
        "NA+": "Na"
    }

    if name in atom_type_map:
        return atom_type_map[name]

    cleaned = name.lstrip("0123456789")

    if not cleaned:
        return ""

    two_letter = {
        "CL": "Cl",
        "BR": "Br",
        "MG": "Mg",
        "CA": "Ca",
        "FE": "Fe",
        "ZN": "Zn",
        "MN": "Mn",
        "CU": "Cu",
        "CO": "Co",
        "NI": "Ni"
    }

    if cleaned[:2] in two_letter:
        return two_letter[cleaned[:2]]

    return cleaned[0]


def convert_atom_line(
    line,
    serial,
    ligand=False
):
    """Convert a PDBQT atom record to a PDB atom record."""

    record = line[0:6].strip()

    if record not in ATOM_RECORDS:
        return None

    try:
        atom_name = line[12:16]

        original_resname = line[17:20]

        original_chain = (
            line[21:22]
            if len(line) >= 22
            else "A"
        )

        original_resseq = (
            line[22:26]
            if len(line) >= 26
            else "1"
        )

        insertion_code = (
            line[26:27]
            if len(line) >= 27
            else " "
        )

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

    except (ValueError, IndexError):
        return None

    if ligand:
        resname = "LIG"
        chain = "Z"
        resseq = "1"
    else:
        resname = original_resname
        chain = original_chain
        resseq = original_resseq

    pdbqt_atom_type = ""

    if len(line) >= 79:
        pdbqt_atom_type = line[77:79].strip()

    element = infer_element(pdbqt_atom_type)

    if not element:
        element = infer_element(atom_name)

    pdb_line = (
        f"{record:<6}"
        f"{serial:5d} "
        f"{atom_name:<4}"
        f" "
        f"{resname:>3} "
        f"{chain:1}"
        f"{resseq:>4}"
        f"{insertion_code:1}"
        f"   "
        f"{x:8.3f}"
        f"{y:8.3f}"
        f"{z:8.3f}"
        f"{1.00:6.2f}"
        f"{0.00:6.2f}"
        f"          "
        f"{element:>2}"
    )

    return pdb_line


def convert_receptor_to_pdb(
    receptor_pdbqt,
    output_pdb
):
    """Convert the receptor PDBQT to PDB."""

    receptor_pdbqt = check_file(
        receptor_pdbqt,
        "Receptor PDBQT"
    )

    atom_count = 0

    with open(
        receptor_pdbqt,
        "r",
        errors="replace"
    ) as inp, open(
        output_pdb,
        "w"
    ) as out:

        for line in inp:

            if not line.startswith(ATOM_RECORDS):
                continue

            converted = convert_atom_line(
                line.rstrip("\n"),
                atom_count + 1,
                ligand=False
            )

            if converted:
                out.write(converted + "\n")
                atom_count += 1

        out.write("END\n")

    if atom_count == 0:
        raise ValueError(
            "Receptor contains no valid ATOM/HETATM records."
        )

    return atom_count


def convert_selected_pose_to_pdb(
    pose,
    output_pdb,
    starting_serial=1
):
    """Convert the selected docking pose to PDB."""

    atom_count = 0

    with open(
        output_pdb,
        "w"
    ) as out:

        for line in pose["records"]:

            if not line.startswith(ATOM_RECORDS):
                continue

            converted = convert_atom_line(
                line,
                starting_serial + atom_count,
                ligand=True
            )

            if converted:

                out.write(converted + "\n")
                atom_count += 1

        out.write("END\n")

    if atom_count == 0:
        raise ValueError(
            "Selected ligand pose contains no valid atoms."
        )

    return atom_count


def create_complex_pdb(
    receptor_pdb,
    ligand_pdb,
    output_pdb
):
    """Combine receptor and ligand PDB files."""

    receptor_atoms = []
    ligand_atoms = []

    with open(receptor_pdb, "r") as fh:

        for line in fh:

            if line.startswith(ATOM_RECORDS):
                receptor_atoms.append(line.rstrip())

    with open(ligand_pdb, "r") as fh:

        for line in fh:

            if line.startswith(ATOM_RECORDS):
                ligand_atoms.append(line.rstrip())

    if not receptor_atoms:
        raise ValueError(
            "Receptor PDB contains no atoms."
        )

    if not ligand_atoms:
        raise ValueError(
            "Ligand PDB contains no atoms."
        )

    with open(output_pdb, "w") as out:

        for line in receptor_atoms:
            out.write(line + "\n")

        for line in ligand_atoms:
            out.write(line + "\n")

        out.write("END\n")

    return len(receptor_atoms), len(ligand_atoms)


def run_plip(
    complex_pdb,
    output_dir,
    plip_executable="plip"
):
    """Run PLIP on a protein-ligand complex."""

    if shutil.which(plip_executable) is None:
        raise RuntimeError(
            f"PLIP executable '{plip_executable}' "
            f"was not found in PATH."
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    command = [
        plip_executable,
        "-f",
        str(complex_pdb.resolve()),
        "-o",
        str(output_dir.resolve()),
        "-x",
        "-t"
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(
            "PLIP failed.\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )

    return result


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Select the best docking pose, "
            "create protein-ligand complexes "
            "and optionally run PLIP."
        )
    )

    parser.add_argument(
        "--receptor",
        required=True,
        help="Receptor PDBQT file"
    )

    parser.add_argument(
        "--ligands",
        nargs="+",
        required=True,
        help="Ligand docking-result PDBQT files"
    )

    parser.add_argument(
        "--output",
        default="MAPK3_PLIP_results",
        help="Output directory"
    )

    parser.add_argument(
        "--run-plip",
        action="store_true",
        help="Run PLIP after creating complexes"
    )

    parser.add_argument(
        "--plip",
        default="plip",
        help="PLIP executable"
    )

    args = parser.parse_args()

    output_root = Path(
        args.output
    ).resolve()

    receptor_dir = output_root / "receptor"
    selected_pose_dir = output_root / "selected_poses"
    complex_dir = output_root / "complexes"
    plip_root = output_root / "plip"

    for directory in [
        receptor_dir,
        selected_pose_dir,
        complex_dir,
        plip_root
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True
        )

    print("\n" + "=" * 70)
    print("DOCKING POSE SELECTION AND PLIP PIPELINE")
    print("=" * 70)

    receptor_pdb = (
        receptor_dir / "receptor.pdb"
    )

    receptor_atoms = convert_receptor_to_pdb(
        args.receptor,
        receptor_pdb
    )

    print(
        f"\nReceptor atoms: {receptor_atoms}"
    )

    summary_rows = []

    for index, ligand_path in enumerate(
        args.ligands,
        start=1
    ):

        ligand_path = Path(
            ligand_path
        ).resolve()

        compound_id = (
            f"Compound{index}"
        )

        compound_name = (
            ligand_path.name
        )

        print("\n" + "-" * 70)
        print(
            f"{compound_id}: {compound_name}"
        )

        try:

            poses = parse_docking_pdbqt(
                ligand_path
            )

            best_pose = select_best_pose(
                poses
            )

            print(
                f"Poses detected: {len(poses)}"
            )

            print(
                f"Best pose: MODEL "
                f"{best_pose['model']}"
            )

            print(
                f"Docking score: "
                f"{best_pose['score']:.3f} kcal/mol"
            )

            selected_pose_pdb = (
                selected_pose_dir /
                f"{compound_id}_best_pose.pdb"
            )

            convert_selected_pose_to_pdb(
                best_pose,
                selected_pose_pdb,
                starting_serial=receptor_atoms + 1
            )

            complex_pdb = (
                complex_dir /
                f"MAPK3_{compound_id}_complex.pdb"
            )

            create_complex_pdb(
                receptor_pdb,
                selected_pose_pdb,
                complex_pdb
            )

            print(
                f"Complex created: "
                f"{complex_pdb.name}"
            )

            plip_output = ""

            if args.run_plip:

                compound_plip_dir = (
                    plip_root / compound_id
                )

                run_plip(
                    complex_pdb,
                    compound_plip_dir,
                    args.plip
                )

                plip_output = str(
                    compound_plip_dir
                )

                print(
                    f"PLIP output: "
                    f"{compound_plip_dir}"
                )

            summary_rows.append({
                "compound": compound_name,
                "compound_id": compound_id,
                "source_pdbqt": str(ligand_path),
                "poses_detected": len(poses),
                "best_model": best_pose["model"],
                "best_score_kcal_mol": best_pose["score"],
                "selected_pose_pdb": str(
                    selected_pose_pdb
                ),
                "complex_pdb": str(
                    complex_pdb
                ),
                "plip_output": plip_output
            })

        except Exception as exc:

            print(
                f"ERROR: {exc}",
                file=sys.stderr
            )

    summary_file = (
        output_root /
        "docking_pose_summary.tsv"
    )

    fields = [
        "compound",
        "compound_id",
        "source_pdbqt",
        "poses_detected",
        "best_model",
        "best_score_kcal_mol",
        "selected_pose_pdb",
        "complex_pdb",
        "plip_output"
    ]

    with open(
        summary_file,
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()

        for row in summary_rows:
            writer.writerow(row)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)

    print(
        f"\nDocking summary:\n{summary_file}"
    )


if __name__ == "__main__":
    main()
