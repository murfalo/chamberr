#!/usr/bin/env python3

import argparse
import os
import parmed
import re
from io import StringIO
from pathlib import Path

###############################################################################
# chamber.py:
#     Convert leaprc (Amber FF format) to CHARMM-readable files (.prm, .rtf,
#     .str)
#
# Outputs:
#     {leaprc_source}.prm, {leaprc_source}.rtf, and {leaprc_source}.str
#
# This script converts AMBER force fields (e.g., RNA.OL3, RNA.Shaw,
# RNA.modrna08) to CHARMM-readable formats for use with, e.g., lambda dynamics.
# It requires AmberTools to be installed and $AMBERHOME environment variable to
# be set. The parameters (.prm) and topologies (.rtf) can be loaded into CHARMM
# either separately or in one go using the created stream (.str) file.
###############################################################################
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Convert Amber leaprc files to CHARMM-readable formats.  "
        "Matches a leaprc of the form "
        "`$AMBERHOME/dat/leap/cmd/leaprc.{leaprc_source}.`"
    )
    parser.add_argument(
        "leaprc_source",
        help="Name of the leaprc file to convert (without the 'leaprc.' prefix)",
    )
    parser.add_argument(
        "--only_unique",
        action="store_true",
        help="Restricts output .prm, .rtf, and .str to only atom types newly "
        "introduced in {leaprc_source}",
    )
    args = parser.parse_args()
    leaprc_source = args.leaprc_source
    
    
    # Identify location of AmberTools installation
    if "AMBERHOME" not in os.environ:
        raise Exception(
            "$AMBERHOME environment variable is not defined. Please install "
            "AmberTools and make sure $AMBERHOME is properly set."
        )
    
    amber_home_path = Path(os.environ["AMBERHOME"]).resolve()
    
    # Create a baseline and combined leaprc file with:
    # 1. ff14SB (baseline, combined)
    # 2. RNA.OL3 (baseline, combined)
    # 3. {leaprc_source} (combined)
    # NOTE(MCA): ParmEd does not have an API for merging leaprcs, so I took a
    #            manual approach and directly combine them before loading the
    #            merged copy.
    print(
        f"Creating combined (protein.ff14SB, RNA.OL3, + {leaprc_source}) leaprcs"
    )
    combined_leaprc_file = StringIO()
    
    amber_leaprc_cmd_path = amber_home_path / "dat" / "leap" / "cmd"
    amber_leaprc_protein_ff14SB_path = (
        amber_leaprc_cmd_path / "leaprc.protein.ff14SB"
    )
    amber_leaprc_rna_OL3_path = amber_leaprc_cmd_path / "leaprc.RNA.OL3"
    amber_leaprc_source_path = amber_leaprc_cmd_path / f"leaprc.{leaprc_source}"
    
    combined_leaprc_paths = [
        amber_leaprc_protein_ff14SB_path,
        amber_leaprc_rna_OL3_path,
        amber_leaprc_source_path,
    ]
    
    # Create combined leaprc
    for leaprc_path in combined_leaprc_paths:
        with open(leaprc_path, "r") as leaprc_file:
            combined_leaprc_file.write(f"# Sourced from {leaprc_path}\n")
            combined_leaprc_file.write(leaprc_file.read())
            combined_leaprc_file.write("\n")
    combined_leaprc_file.seek(0)
    
    
    print("Loading combined parameters...")
    amber_combined_params = parmed.amber.AmberParameterSet.from_leaprc(
        combined_leaprc_file
    )
    combined_atom_types = [
        str(atom_type) for atom_type in amber_combined_params.atom_types
    ]
    
    
    # Remove single-entry, ill-defined atom types to prevent errors when converting
    # to CHARMM
    print(
        "Removing poorly defined atom types... Make sure you don't see anything "
        "you need!"
    )
    for atom_name, atom_type in list(amber_combined_params.atom_types.items()):
        # NOTE(MCA): Feel free to expand this check to other properties if you get
        # errors like "type.<your_property> is None" during conversion.
        if not atom_type.epsilon:
            print(
                f"WARNING: Atom type ({atom_type}) is not fully initialized.  "
                f"Removing it..."
            )
            amber_combined_params.atom_types.pop(atom_name, None)
    
    
    # Convert to CHARMM parameters
    print("Converting to CHARMM parameters...")
    charmm_params = parmed.charmm.CharmmParameterSet.from_parameterset(
        amber_combined_params
    )
    
    
    # Ensure the output directory exists
    output_dir = Path("./toppar").resolve()
    print(f"Creating output directory ({output_dir})...")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    
    # Export CHARMM .prm (parameters), .rtf (topologies), and an all-in-one .str
    # (stream file)
    charmm_prm_path = output_dir / f"{leaprc_source}.prm"
    charmm_rtf_path = output_dir / f"{leaprc_source}.rtf"
    charmm_str_path = output_dir / f"{leaprc_source}.str"
    print(
        f"Saving CHARMM-readable {leaprc_source}.prm, {leaprc_source}.rtf, and "
        f"{leaprc_source}.str in {output_dir}..."
    )
    
    charmm_params.write(
        par=str(charmm_prm_path),
        top=str(charmm_rtf_path),
    )
    
    charmm_params.write(
        stream=str(charmm_str_path),
    )
    
    
    # Logic for identifying new atom types restricted to after this point
    if not args.only_unique:
        exit()
    
    
    # Create baseline leaprc
    # NOTE(MCA): Baseline exists to help us identify which atom types are newly
    # introduced in {leaprc_source}
    print(
        "Only unique selected. Creating baseline (protein.ff14SB and RNA.OL3) "
        "leaprc..."
    )
    baseline_leaprc_file = StringIO()
    baseline_leaprc_paths = [
        leaprc_path
        for leaprc_path in combined_leaprc_paths
        if leaprc_path != amber_leaprc_source_path
    ]
    for leaprc_path in baseline_leaprc_paths:
        with open(leaprc_path, "r") as leaprc_file:
            baseline_leaprc_file.write(f"# Sourced from {leaprc_path}\n")
            baseline_leaprc_file.write(leaprc_file.read())
            baseline_leaprc_file.write("\n")
    baseline_leaprc_file.seek(0)
    
    
    print("Identifying new atom types...")
    amber_baseline_params = parmed.amber.AmberParameterSet.from_leaprc(
        baseline_leaprc_file
    )
    baseline_atom_types = [
        str(atom_type) for atom_type in amber_baseline_params.atom_types
    ]
    amber_source_params = parmed.amber.AmberParameterSet.from_leaprc(
        str(amber_leaprc_source_path)
    )
    new_atom_types = [
        str(atom_type)
        for atom_type in amber_source_params.atom_types
        if atom_type not in baseline_atom_types
    ]
    
    
    # Eliminate all entries that do not contain new atom types
    print("Restricting parameters and topologies to new atom types...")
    charmm_paths = [charmm_str_path, charmm_prm_path, charmm_rtf_path]
    
    combined_atom_types_regex = re.compile(
        r"(?:^|\s)(" + "|".join(map(re.escape, combined_atom_types)) + r")(?:\s|$)"
    )
    new_atom_types_regex = re.compile(
        r"(?:^|\s)(" + "|".join(map(re.escape, new_atom_types)) + r")(?:\s|$)"
    )
    
    print(f"New atom types are: {new_atom_types}")
    
    for charmm_path in charmm_paths:
        with open(charmm_path, "r") as charmm_file:
            lines = charmm_file.readlines()
    
        # Remove any lines that contain exclusively old atom types
        filtered_lines = [
            line
            for line in lines
            if new_atom_types_regex.search(line)
            or not combined_atom_types_regex.search(line)
        ]
    
        with open(charmm_path, "w") as charmm_file:
            charmm_file.writelines(filtered_lines)

if __name__ == "__main__":
    main()
