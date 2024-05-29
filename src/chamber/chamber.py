#!/usr/bin/env python3

import argparse
import contextlib
import os
import re
from io import StringIO
from pathlib import Path
from typing import Optional, List, IO

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
            description="Convert Amber leaprc files to CHARMM-readable "
            "formats.  Input leaprc arguments are expected to be found at "
            "`$AMBERHOME/dat/leap/cmd/leaprc.*`",
        )
        parser.add_argument(
            "-p",
            "--protein_leaprc",
            help="Name of the Amber protein force field leaprc to load "
            "(without the `leaprc.protein.` prefix)",
            default="ff14sb",
        )
        parser.add_argument(
            "-r",
            "--rna_leaprc",
            help="Name of the Amber RNA force field leaprc to load (without "
            "the `leaprc.RNA.` prefix)",
            default="OL3",
        )
        parser.add_argument(
            "-s",
            "--solvent_leaprc",
            help="Name of the Amber solvent force field leaprc to load "
            "(without the `leaprc.water.` prefix)",
            default="tip3p",
        )
        parser.add_argument(
            "extra_leaprc",
            help="Name of an additional Amber force field leaprc to load "
            "(without the `leaprc.` prefix)",
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
        ("ZN2P", "ZN  "),
        # KP -> POT must be replaced in two steps to ensure whitespace/columns
        # are preserved.
        (r"KP\s(\s+)", r"POT\1"),
        (r"KP", "POT"),
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
    passing :paramref:`extra_leaprc`.

    Chamber was built to convert RNA force fields (specifically, RNA.OL3 and
    RNA.modrna08) to CHARMM-readable formats for use with, e.g., λ-dynamics.
    These are currently the only rigorously tested use cases.  Expanding to
    other force fields should be possible, but requires additional care.

    :param protein_leaprc: The Amber protein force field leaprc to load.
     Defaults to "ff14sb".
    :param rna_leaprc: The Amber RNA force field leaprc to load.  Defaults to
     "OL3".
    :param solvent_leaprc: The Amber solvent force field leaprc to load.
     Defaults to "tip3p".
    :param extra_leaprc: An optional, additional leaprc source to load.
    :ivar __protein_leaprc: A path to the input protein leaprc source.
    :ivar __rna_leaprc: A path to the input RNA leaprc source.
    :ivar __solvent_leaprc: A path to the input solvent leaprc source.
    :ivar __extra_leaprc: A path to the input extra leaprc source.
    :ivar __combined_amber_params: An AmberParameterSet containing the
     combined parameters.
    """

    __SLOTS__ = [
        "__protein_leaprc",
        "__rna_leaprc",
        "__solvent_leaprc",
        "__extra_leaprc",
        "__combined_amber_params",
    ]

    __protein_leaprc: Path
    __rna_leaprc: Path
    __solvent_leaprc: Path
    __extra_leaprc: Optional[Path]
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
        f"* Automatically generated with Chamber\n"
        f"* https://github.com/murfalo/chamber\n"
        f"*\n"
        f"\n"
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

    @property
    def leaprc_path(self) -> Path:
        """
        Returns a Path object describing the where Amber's leaprcs are located
        ($AMBERHOME/dat/leap/cmd).

        :return: The Amber leaprc path.
        """
        assert "AMBERHOME" in os.environ, "$AMBERHOME is not set"
        return Path(os.environ["AMBERHOME"]).resolve() / "dat" / "leap" / "cmd"

    @property
    def __leaprc_paths(self) -> List[Path]:
        """
        Returns a tuple of Chamber's leaprc paths in ascending order of
        priority.

        :return: The list of leaprc paths.
        """
        # In order of ascending priority (parameters defined in lower priority
        # files will overwrite parameters in higher priority files)
        leaprcs = [
            self.__extra_leaprc,
            self.__rna_leaprc,
            self.__protein_leaprc,
            self.__solvent_leaprc,
        ]

        if self.__extra_leaprc is None:
            return leaprcs[1:]
        else:
            return leaprcs

    @property
    def combined_leaprc(self) -> IO[str]:
        """
        Returns an in-memory, combined leaprc file which sources all of
        Chamber's leaprc sources in order of ascending priority.

        :return: The in-memory, combined leaprc.
        """
        # Create an in-memory leaprc file which combines all of the above
        # sources using the `source` command.
        combined_leaprc_file = StringIO()
        for leaprc_path in self.__leaprc_paths:
            with open(leaprc_path, "r", encoding="utf-8") as leaprc_file:
                combined_leaprc_file.write(f"# Sourced from {leaprc_path}\n")
                combined_leaprc_file.write(leaprc_file.read())
                combined_leaprc_file.write("\n")
        combined_leaprc_file.seek(0)

        return combined_leaprc_file

    def __init__(
        self,
        protein_leaprc: str = "ff14sb",
        rna_leaprc: str = "OL3",
        solvent_leaprc: str = "tip3p",
        extra_leaprc: str = None,
    ):
        """
        Initializes a Chamber instance.  See :class:`Chamber` documentation.
        """
        leaprc_path = self.leaprc_path

        # Identify force field (protein, RNA, solvent, extra source) source
        # names and paths
        self.__protein_leaprc = (
            leaprc_path / f"leaprc.protein.{protein_leaprc}"
        )
        self.__rna_leaprc = leaprc_path / f"leaprc.RNA.{rna_leaprc}"
        self.__solvent_leaprc = leaprc_path / f"leaprc.water.{solvent_leaprc}"
        self.__extra_leaprc = (
            leaprc_path / f"leaprc.{extra_leaprc}" if extra_leaprc else None
        )

        # Load combined AmberParameterSet from the combined leaprc
        print("Loading combined parameters...")
        load_leaprc_stderr = StringIO()
        with contextlib.redirect_stderr(load_leaprc_stderr):
            self.__combined_amber_params = AmberParameterSet.from_leaprc(
                self.combined_leaprc
            )

        load_leaprc_stderr = load_leaprc_stderr.getvalue()
        if load_leaprc_stderr != "":
            print(
                f"WARNING!  Warnings and/or errors occurred when loading "
                f"AmberParameterSet\n{load_leaprc_stderr}"
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
    def protein_leaprc_name(self) -> str:
        """
        Returns the name of the protein leaprc (without the `leaprc.protein.`
        prefix).
        """
        return self.__protein_leaprc.name.removeprefix("leaprc.protein.")

    @property
    def rna_leaprc_name(self) -> str:
        """
        Returns the name of the RNA leaprc (without the `leaprc.RNA.` prefix).
        """
        return self.__rna_leaprc.name.removeprefix("leaprc.RNA.")

    @property
    def solvent_leaprc_name(self) -> str:
        """
        Returns the name of the solvent leaprc (without the `leaprc.solvent.`
        prefix).
        """
        return self.__solvent_leaprc.name.removeprefix("leaprc.water.")

    @property
    def extra_leaprc_name(self) -> str:
        """
        Returns the name of the extra leaprc (without the `leaprc.` prefix).
        """
        if self.__extra_leaprc:
            return self.__extra_leaprc.name.removeprefix("leaprc.")
        else:
            return "none"

    @property
    def charmm_str_name(self) -> str:
        """
        Returns the name of the CHARMM stream file that Chamber outputs. This
        will be `ff14sb_chamber.str` if extra_leaprc is None, otherwise the
        provided extra leaprc source's name with the suffix `.str`

        :return: The CHARMM .str filename.
        """
        leaprc_names = (
            self.extra_leaprc_name,
            self.rna_leaprc_name,
            self.protein_leaprc_name,
            self.solvent_leaprc_name,
        )

        return f"chamber_{'_'.join(leaprc_names)}.str"

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
                print(
                    f"WARNING! Additional angle ({new_angle_key}) already "
                    f"exists."
                )
            else:
                charmm_params.angle_types[new_angle_key] = new_angle_type
            if (
                new_angle_key_reverse != new_angle_key
                and new_angle_key_reverse in charmm_params.angle_types
            ):
                print(
                    f"WARNING! Reverse of additional angle "
                    f"({new_angle_key_reverse}) already exists."
                )
            else:
                charmm_params.angle_types[new_angle_key_reverse] = (
                    new_angle_type
                )

        # Add additional atoms
        for atom_type, lj_params in Chamber.__ADDITIONAL_ATOMS:
            atom_type.set_lj_params(
                lj_params[0], lj_params[1], lj_params[2], lj_params[3]
            )
            if atom_type not in charmm_params.atom_types:
                charmm_params.atom_types[atom_type.name] = atom_type
            else:
                print(f"WARNING! Additional atom {atom_type} already exists.")

        # Create an in-memory CHARMM stream file
        print("Writing in-memory CHARMM stream file...")
        charmm_str = StringIO()
        charmm_params.write(stream=charmm_str)
        charmm_str.seek(0)
        charmm_str = (
            self.__HEADER
            + f"! Force fields loaded in ascending order of priority\n"
            f"! Extra source: {self.extra_leaprc_name}\n"
            f"! RNA source: {self.rna_leaprc_name}\n"
            f"! Protein source: {self.protein_leaprc_name}\n"
            f"! Solvent source: {self.solvent_leaprc_name}\n\n"
            + charmm_str.getvalue()
        )

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
    chamber = Chamber(
        protein_leaprc=args.protein_leaprc,
        rna_leaprc=args.rna_leaprc,
        solvent_leaprc=args.solvent_leaprc,
        extra_leaprc=args.extra_leaprc,
    )
    chamber.write_charmm_parameters()


if __name__ == "__main__":
    main()
