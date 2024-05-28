#!/usr/bin/env python3

import argparse
import contextlib
import os
import re
from io import StringIO
from pathlib import Path
from typing import Optional

from parmed import AngleType, AtomType
from parmed.amber import AmberParameterSet
from parmed.charmm import CharmmParameterSet


class _Parser:
    """
    Parser for Chamber.
    """

    @staticmethod
    def __build_parser() -> argparse.ArgumentParser:
        """
        Builds a parser for Chamber.
        """
        parser = argparse.ArgumentParser(
            prog="chamber",
            description="Convert Amber leaprc files to CHARMM-readable formats.  "
                        "Matches a leaprc of the form "
                        "`$AMBERHOME/dat/leap/cmd/leaprc.{leaprc_source}.`",
        )
        parser.add_argument(
            "leaprc_source",
            help="name of the leaprc file to convert (without the 'leaprc.' "
                 "prefix)",
            nargs="?",
            default=None,
        )

        return parser

    @staticmethod
    def parse_args() -> argparse.Namespace:
        """
        Parses arguments for Chamber.

        :return: The parsed arguments.
        """
        return _Parser.__build_parser().parse_args()


class _Regex:
    """
    Regular expression patterns for user with Chamber.
    """

    @staticmethod
    def tokenize(regex_pattern: str) -> str:
        """
        Converts a regular expression pattern to match only full tokens.
        On the left, this requires either a word boundary or start of line.
        On the right, this requires either a word boundary or end of line.
        This prevents replacing partial tokens, such as replacing
        "MASS...CLAABC" due to a false-positive matching of "CLM" -> "CLA".

        :param regex_pattern: The original regular expression pattern.
        :return: The tokenized regular expression pattern.
        """
        return r"(?:^|\b)" + regex_pattern + r"(?:\b|$)"

    # Regex for identifying atom type names
    ATOM_TYPE = re.compile(tokenize(r"MASS [ 0-9]{5} ([A-Za-z0-9]{1,6})"))

    # Regex substitutions to conduct on the final `.str` file to ensure it is
    # fully compatible with the official CHARMM FF14SB topologies.  These are
    # applied in order, so earlier matches may effect later matches.
    RE_SUBS = (
        # The following atoms have been renamed to match the names found in the
        # official CHARMM FF14SB ported topology.
        ("CSTR", "CG  "),
        ("NSTR", "NG  "),
        ("CLM", "CLA"),
        ("NAP", "SOD"),
        ("KP", "POT"),
        ("ZN2P", "ZN  "),
        # This dihedral must be reordered to match definitions in
        # the official CHARMM FF14SB port.
        (
            "C5     CB     NG     CT",
            "CB     C5     NG     CT",
        ),
        # Use append for reading topologies, and flex for reading parameters.
        # `append` prevents our skeleton topology, which includes only atom
        # definitions, from overwriting any previously loaded topologies.
        # See: https://academiccharmm.org/documentation/latest/io#RTFFileformat
        # See: https://academiccharmm.org/documentation/latest/io#ParamFiles
        ("read rtf card", "read rtf append card"),
        ("read para card", "read para flex card"),
        # Setting the atom type to -1 allows CHARMM to automatically assign
        # type indices.  This is particularly important when using `append`,
        # as it appears to let CHARMM assign a single type index to all MASS
        # entries with the same type name.
        ("MASS [ 0-9]{5}", "MASS    -1"),
    )


# Compile each RE_SUBS pattern for efficiency
_Regex.RE_SUBS = tuple(
    (re.compile(_Regex.tokenize(find)), replace)
    for find, replace in _Regex.RE_SUBS
)


class Chamber:
    """
    Convert leaprc (Amber FF format) to a combined, CHARMM-readable parameter
    and topology files in `.str` format.  Currently, the output topology only
    includes the required masses and atom types.  To use the output parameters,
    load them into CHARMM using `stream <chamber_output>.str` along with a
    topology file (`.rtf`) defining residues and patches with the
    Chamber-generated atom types.  By default, Chamber loads protein.ff14SB,
    RNA.OL3, and water.tip3p.  An additional forcefield can be included by
    passing :paramref:`extra_leaprc_source`.

    Chamber was built to convert RNA force fields (specifically, RNA.OL3 and
    RNA.modrna08) to CHARMM-readable formats for use with, e.g., λ-dynamics.
    These are currently the only rigorously tested use cases.  Expanding to
    other force fields should be possible, but requires additional care.

    :param extra_leaprc_source: An optional, additional leaprc source to load.
    :ivar __extra_leaprc_source: A path to the input extra leaprc source.
    :ivar __combined_amber_params: An AmberParameterSet containing the
     combined parameters.
    """

    __SLOTS__ = ["__extra_leaprc_source", "__combined_amber_params"]

    # TODO(MCA): Allow Optional[Tuple[Path]]?
    __extra_leaprc_source: Optional[Path]
    __combined_amber_params: AmberParameterSet

    __MAX_TYPE_NAME_LENGTH: int = 6  # Type names can be at most 6 characters
    __ADDITIONAL_ANGLES = [
        # (AngleType Key, AngleType(phi_k, thet_eq))
        # CHARMM requires an angle for TIP3 waters, even though they are
        # rigid (Amber leaves them undefined)
        (("HW", "OW", "HW"), AngleType(100.00, 104.52)),
    ]
    __ADDITIONAL_ATOMS = (
        # (AtomType, (epsilon, rmin, epsilon_14, rmin_14))
        (AtomType("EP", None, 0.0), (0.0, 0.0, 0.0, 0.0)),
    )
    __ADDITIONAL_COMMANDS = (
        # In both the Amber and CHARMM force fields, hydrogen bonds are
        # implicit.  Therefore, explicit hydrogen bonds should be disabled by
        # setting CUTHB to anything less than 1.0.
        # See: https://academiccharmm.org/documentation/latest/hbonds#Function
        "HBOND CUTHB 0.5",
    )
    __HEADER: str = (
        "* Automatically generated with Chamber\n"
        "* https://github.com/murfalo/chamber\n\n"
    )

    def __cleanup(self):
        """
        Helper function for cleaning up Chamber's AmberParameterSet.
        """
        # Remove single-entry, ill-defined atom types to prevent errors when
        # converting to CHARMM.
        print(
            "Removing poorly defined atom types... Make sure you don't see "
            "anything you need!"
        )
        for (
                atom_name,
                atom_type,
        ) in list(self.__combined_amber_params.atom_types.items()):
            if atom_type.epsilon is None:
                print(
                    f"WARNING: Atom type ({atom_type}) is not fully "
                    f"initialized.  Removing it..."
                )
                del self.__combined_amber_params.atom_types[atom_name]

    def __init__(self, extra_leaprc_source: str = None):
        """
        Initializes a Chamber instance.  See :class:`Chamber` documentation.
        """
        assert "AMBERHOME" in os.environ, "AMBERHOME is not set"
        amber_home_path = Path(os.environ["AMBERHOME"]).resolve()

        # Identify force field (protein, RNA, solvent, extra source) source
        # names and paths
        # TODO(MCA): Add commandline options for selecting protein, RNA, and/or
        #  solvent forcefields?
        protein_source = "protein.ff14sb"
        rna_source = "RNA.OL3"
        solvent_source = "water.tip3p"
        combined_sources = (
            extra_leaprc_source,
            rna_source,
            protein_source,
            solvent_source,
        )

        # If extra leaprc source is None, remove it from the combined sources.
        # Otherwise, store a path to it for later use.
        leaprc_cmd_path = amber_home_path / "dat" / "leap" / "cmd"
        if extra_leaprc_source is None:
            self.__extra_leaprc_source = None
            combined_sources = combined_sources[1:]
        else:
            self.__extra_leaprc_source = (
                    leaprc_cmd_path / f"leaprc.{extra_leaprc_source}"
            )

        # In order of ascending priority (parameters defined in lower priority
        # files will overwrite parameters in higher priority files)
        print(f"Creating combined ({', '.join(combined_sources)}) leaprc...")
        combined_leaprc_paths = tuple(
            leaprc_cmd_path / f"leaprc.{source_name}"
            for source_name in combined_sources
        )

        # Create an in-memory leaprc file which combines all of the above
        # sources using the `source` command.
        combined_leaprc_file = StringIO()
        for leaprc_path in combined_leaprc_paths:
            with open(leaprc_path, "r", encoding="utf-8") as leaprc_file:
                combined_leaprc_file.write(f"# Sourced from {leaprc_path}\n")
                combined_leaprc_file.write(leaprc_file.read())
                combined_leaprc_file.write("\n")
        combined_leaprc_file.seek(0)

        # Load combined AmberParameterSet from the combined leaprc
        print("Loading combined parameters...")
        load_leaprc_stderr = StringIO()
        with contextlib.redirect_stderr(load_leaprc_stderr):
            self.__combined_amber_params = AmberParameterSet.from_leaprc(
                combined_leaprc_file
            )

        load_leaprc_stderr = load_leaprc_stderr.getvalue()
        if load_leaprc_stderr != "":
            raise Exception(
                f"Failed to load AmberParameterSet\n{load_leaprc_stderr}"
            )

        self.__cleanup()

    @staticmethod
    def __remove_type_name_decor(type_name: str) -> str:
        """
        ParmEd output Amber type names include P for "+", M for "-", and "LTU"
        ("Lower To Upper") when a type name originally contained lowercase
        characters.  This function removes the "LTU" decorator.

        :param type_name: The type name as a string.
        :return: The type name with "LTU" decorator removed.
        """
        match_len = len(type_name)

        # Removing trailing LTU
        if (
                type_name.endswith("LT")
                and match_len == Chamber.__MAX_TYPE_NAME_LENGTH
        ):
            # Trailing U was clipped off due to name length restriction
            type_name = type_name[:-2]
        elif type_name.endswith("LTU"):
            type_name = type_name[:-3]

        # Left justify with whitespaces to preserve columns/indentation
        return type_name.ljust(match_len, " ")

    @property
    def charmm_str_name(self) -> str:
        """
        Returns the name of the CHARMM stream file that Chamber outputs. This
        will be `ff14sb_chamber.str` if extra_leaprc_source is None, otherwise
        the provided extra leaprc source's name with the suffix `.str`

        :return: The CHARMM .str filename.
        """
        if self.__extra_leaprc_source is None:
            output_name = "ff14sb"
        else:
            # Combines
            split_output_stem = self.__extra_leaprc_source.name.split(".")
            assert (
                    split_output_stem[0] == "leaprc"
            ), "Extra leaprc source must be a leaprc file"
            output_name = "".join(
                self.__extra_leaprc_source.name.split(".")[1:]
            )

        return f"chamber_{output_name}.str"

    def write_charmm_parameters(self):
        """
        Write chamber's combined AmberParameterSet into a CHARMM-readable
        `.str` file.
        """
        # Convert to CHARMM parameters
        print("Converting to CHARMM parameters...")
        charmm_params = CharmmParameterSet.from_parameterset(
            self.__combined_amber_params
        )

        # Add additional angles
        print("Adding additional angles...")
        for new_angle_key, new_angle_type in Chamber.__ADDITIONAL_ANGLES:
            new_angle_key_reverse = new_angle_key[::-1]
            if new_angle_key in charmm_params.angle_types:
                raise KeyError(
                    f"Additional angle ({new_angle_key}) already exists. "
                    f"Please modify Chamber.__ADDITIONAL_ANGLES."
                )
            charmm_params.angle_types[new_angle_key] = new_angle_type
            if (
                    new_angle_key_reverse != new_angle_key
                    and new_angle_key_reverse in charmm_params.angle_types
            ):
                raise KeyError(
                    f"Reverse of additional angle ({new_angle_key_reverse}) "
                    f"already exists. Please modify "
                    f"Chamber.__ADDITIONAL_ANGLES."
                )
            charmm_params.angle_types[new_angle_key_reverse] = new_angle_type

        # Add additional atoms
        for atom_type, lj_params in Chamber.__ADDITIONAL_ATOMS:
            if atom_type not in charmm_params.atom_types:
                atom_type.set_lj_params(
                    lj_params[0], lj_params[1], lj_params[2], lj_params[3]
                )
                charmm_params.atom_types[atom_type.name] = atom_type

        # Create an in-memory CHARMM stream file
        print("Writing in-memory CHARMM stream file...")
        charmm_str = StringIO()
        charmm_params.write(stream=charmm_str)
        charmm_str.seek(0)
        charmm_str = self.__HEADER + charmm_str.getvalue()

        # Substitute atom types
        print("Removing decorators from type names...")
        type_names = _Regex.ATOM_TYPE.findall(charmm_str)
        type_name_map = {
            name: Chamber.__remove_type_name_decor(name) for name in type_names
        }
        for original_name, new_name in type_name_map.items():
            charmm_str = re.sub(
                _Regex.tokenize(re.escape(original_name)), new_name, charmm_str
            )

        # Apply miscellaneous Regex substitutions
        print("Applying regular expression substitutions...")
        for pattern, replace in _Regex.RE_SUBS:
            charmm_str = pattern.sub(replace, charmm_str)

        charmm_str_lines = charmm_str.splitlines()

        # Add additional commands just before "END"
        print("Inserting additional commands...")
        found_end = False
        for index in range(len(charmm_str_lines) - 1, -1, -1):
            if charmm_str_lines[index].strip() == "END":
                charmm_str_lines[index:index] = Chamber.__ADDITIONAL_COMMANDS
                found_end = True
                break

        if not found_end:
            raise IndexError(
                f'Could not find "END" in charmm stream file...  Please '
                f"correct this before continuing."
            )

        # Write CHARMM stream file
        charmm_str = Path(".").resolve() / self.charmm_str_name
        print(f"Saving CHARMM stream file to disk ({charmm_str})...")
        with open(charmm_str, "w", encoding="utf-8") as charmm_str:
            charmm_str.write("\n".join(charmm_str_lines))


def main():
    # Identify location of AmberTools installation
    if "AMBERHOME" not in os.environ:
        raise KeyError(
            "$AMBERHOME environment variable is not defined. Please install "
            "AmberTools and make sure $AMBERHOME is properly set."
        )

    args = _Parser.parse_args()
    chamber = Chamber(extra_leaprc_source=args.leaprc_source)
    chamber.write_charmm_parameters()


if __name__ == "__main__":
    main()
