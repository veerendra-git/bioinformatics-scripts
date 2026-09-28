#!/usr/bin/env python3

"""
Run PLIP on protein-ligand complexes and extract interactions.

The script:
1. Reads docking scores from docking_pose_summary.tsv
2. Finds MAPK3-ligand complex PDB files
3. Runs PLIP for each complex
4. Reads the PLIP XML reports
5. Extracts interaction details
6. Extracts ligand properties
7. Creates summary TSV files
"""

import csv
import re
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET


def read_docking_scores(summary_file):
    """Read docking scores from the docking summary."""

    scores = {}

    if not summary_file.exists():
        return scores

    with open(
        summary_file,
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            compound_id = row.get(
                "compound_id",
                ""
            )

            if compound_id:

                scores[compound_id] = {
                    "ligand": row.get(
                        "compound",
                        ""
                    ),
                    "score": row.get(
                        "best_score_kcal_mol",
                        ""
                    )
                }

    return scores


def find_xml(folder):
    """Find the PLIP XML report."""

    xml_files = list(
        folder.glob("*_report.xml")
    )

    if not xml_files:
        return None

    return xml_files[0]


def get_text(parent, tag):
    """Get text from an XML element."""

    element = parent.find(tag)

    if element is None or element.text is None:
        return ""

    return element.text.strip()


def run_plip(
    complex_file,
    output_dir,
    plip_executable="plip"
):
    """Run PLIP on one complex."""

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    existing_xml = find_xml(
        output_dir
    )

    if existing_xml:

        print(
            f"  PLIP already completed: "
            f"{existing_xml.name}"
        )

        return existing_xml

    command = [
        plip_executable,
        "-f",
        str(complex_file),
        "-o",
        str(output_dir),
        "-x",
        "-t"
    ]

    print("  Running PLIP...")

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:

        print(
            "  ERROR: PLIP failed"
        )

        print(
            result.stderr
        )

        return None

    xml_file = find_xml(
        output_dir
    )

    if xml_file:

        print(
            f"  PLIP complete: "
            f"{xml_file.name}"
        )

    else:

        print(
            "  ERROR: PLIP completed "
            "but XML report was not found"
        )

    return xml_file


def parse_xml(
    xml_file,
    compound,
    docking_score
):
    """Extract interactions and ligand properties from PLIP XML."""

    tree = ET.parse(
        xml_file
    )

    root = tree.getroot()

    interactions = []
    properties = {}

    binding_site = root.find(
        ".//bindingsite"
    )

    if binding_site is None:
        return interactions, properties

    # -----------------------------
    # Ligand information
    # -----------------------------

    identifiers = binding_site.find(
        "identifiers"
    )

    if identifiers is not None:

        properties["Ligand"] = get_text(
            identifiers,
            "longname"
        )

        properties["Ligand_Type"] = get_text(
            identifiers,
            "ligtype"
        )

        properties["Ligand_Chain"] = get_text(
            identifiers,
            "chain"
        )

        properties["Ligand_Position"] = get_text(
            identifiers,
            "position"
        )

        properties["SMILES"] = get_text(
            identifiers,
            "smiles"
        )

        properties["InChIKey"] = get_text(
            identifiers,
            "inchikey"
        )

    # -----------------------------
    # Ligand properties
    # -----------------------------

    lig_props = binding_site.find(
        "lig_properties"
    )

    if lig_props is not None:

        properties["Heavy_Atoms"] = get_text(
            lig_props,
            "num_heavy_atoms"
        )

        properties["HBD"] = get_text(
            lig_props,
            "num_hbd"
        )

        properties["HBA"] = get_text(
            lig_props,
            "num_hba"
        )

        properties["Aromatic_Rings"] = get_text(
            lig_props,
            "num_aromatic_rings"
        )

        properties["Rotatable_Bonds"] = get_text(
            lig_props,
            "num_rotatable_bonds"
        )

        properties["Molecular_Weight"] = get_text(
            lig_props,
            "molweight"
        )

        properties["LogP"] = get_text(
            lig_props,
            "logp"
        )

    # -----------------------------
    # Interaction information
    # -----------------------------

    interactions_root = binding_site.find(
        "interactions"
    )

    if interactions_root is None:
        return interactions, properties

    interaction_types = [
        "hydrophobic_interactions",
        "hydrogen_bonds",
        "pi_stacks",
        "pi_cation_interactions",
        "salt_bridges",
        "halogen_bonds",
        "water_bridges"
    ]

    for group_name in interaction_types:

        group = interactions_root.find(
            group_name
        )

        if group is None:
            continue

        for interaction in list(group):

            row = {
                "Compound": compound,
                "Docking_Score": docking_score,
                "Interaction_Type": interaction.tag
            }

            for child in interaction:

                key = child.tag

                value = (
                    child.text.strip()
                    if child.text
                    else ""
                )

                row[key] = value

            interactions.append(
                row
            )

    return interactions, properties


def write_tsv(
    filename,
    rows,
    columns
):
    """Write rows to a tab-separated file."""

    with open(
        filename,
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Run PLIP on protein-ligand "
            "complexes and extract interactions."
        )
    )

    parser.add_argument(
        "--input",
        default="MAPK3_PLIP_results",
        help="Directory containing complexes and docking summary"
    )

    parser.add_argument(
        "--plip",
        default="plip",
        help="PLIP executable"
    )

    args = parser.parse_args()

    base = Path(
        args.input
    ).resolve()

    complex_dir = (
        base / "complexes"
    )

    plip_dir = (
        base / "PLIP"
    )

    docking_summary = (
        base / "docking_pose_summary.tsv"
    )

    interactions_tsv = (
        base / "PLIP_interactions.tsv"
    )

    summary_tsv = (
        base / "PLIP_summary.tsv"
    )

    properties_tsv = (
        base / "PLIP_ligand_properties.tsv"
    )

    plip_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    docking_scores = read_docking_scores(
        docking_summary
    )

    complex_files = sorted(
        complex_dir.glob(
            "MAPK3_Compound*_complex.pdb"
        )
    )

    if not complex_files:

        print(
            "ERROR: No complex PDB files found."
        )

        return

    all_interactions = []
    all_properties = []
    summary_rows = []

    print("=" * 70)
    print("MAPK3 - PLIP ANALYSIS")
    print("=" * 70)

    print(
        f"\nFound {len(complex_files)} complexes."
    )

    for complex_file in complex_files:

        match = re.search(
            r"MAPK3_(Compound\d+)_complex\.pdb",
            complex_file.name
        )

        if not match:
            continue

        compound = match.group(1)

        print(
            f"\n[{compound}]"
        )

        print(
            f"  Complex: "
            f"{complex_file.name}"
        )

        output_dir = (
            plip_dir / compound
        )

        xml_file = run_plip(
            complex_file,
            output_dir,
            args.plip
        )

        if xml_file is None:
            continue

        docking_info = docking_scores.get(
            compound,
            {}
        )

        ligand_name = docking_info.get(
            "ligand",
            ""
        )

        docking_score = docking_info.get(
            "score",
            ""
        )

        interactions, properties = parse_xml(
            xml_file,
            compound,
            docking_score
        )

        properties["Compound"] = compound
        properties["Ligand_Name"] = ligand_name
        properties["Docking_Score"] = docking_score

        all_properties.append(
            properties
        )

        all_interactions.extend(
            interactions
        )

        counts = {}

        for row in interactions:

            interaction_type = row[
                "Interaction_Type"
            ]

            counts[interaction_type] = (
                counts.get(
                    interaction_type,
                    0
                ) + 1
            )

        summary_row = {
            "Compound": compound,
            "Docking_Score": docking_score,
            "Total_Interactions": len(
                interactions
            )
        }

        summary_row.update(
            counts
        )

        summary_rows.append(
            summary_row
        )

        print(
            f"  Interactions detected: "
            f"{len(interactions)}"
        )

    # -----------------------------
    # Interaction table
    # -----------------------------

    if all_interactions:

        columns = set()

        for row in all_interactions:
            columns.update(
                row.keys()
            )

        preferred = [
            "Compound",
            "Docking_Score",
            "Interaction_Type",
            "resnr",
            "restype",
            "reschain",
            "resnr_lig",
            "restype_lig",
            "reschain_lig",
            "dist",
            "distance"
        ]

        columns = (
            [
                col
                for col in preferred
                if col in columns
            ]
            +
            sorted(
                col
                for col in columns
                if col not in preferred
            )
        )

        write_tsv(
            interactions_tsv,
            all_interactions,
            columns
        )

    # -----------------------------
    # Interaction summary
    # -----------------------------

    if summary_rows:

        columns = set()

        for row in summary_rows:
            columns.update(
                row.keys()
            )

        first_columns = [
            "Compound",
            "Docking_Score",
            "Total_Interactions"
        ]

        columns = (
            first_columns
            +
            sorted(
                col
                for col in columns
                if col not in first_columns
            )
        )

        write_tsv(
            summary_tsv,
            summary_rows,
            columns
        )

    # -----------------------------
    # Ligand properties
    # -----------------------------

    if all_properties:

        columns = [
            "Compound",
            "Docking_Score",
            "Ligand_Name",
            "Ligand",
            "Ligand_Type",
            "Ligand_Chain",
            "Ligand_Position",
            "Heavy_Atoms",
            "HBD",
            "HBA",
            "Aromatic_Rings",
            "Rotatable_Bonds",
            "Molecular_Weight",
            "LogP",
            "SMILES",
            "InChIKey"
        ]

        write_tsv(
            properties_tsv,
            all_properties,
            columns
        )

    print("\n" + "=" * 70)
    print("PLIP ANALYSIS COMPLETE")
    print("=" * 70)

    print(
        f"\nInteraction table:\n"
        f"{interactions_tsv}"
    )

    print(
        f"\nSummary table:\n"
        f"{summary_tsv}"
    )

    print(
        f"\nLigand properties:\n"
        f"{properties_tsv}"
    )


if __name__ == "__main__":
    main()
