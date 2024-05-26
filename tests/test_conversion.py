#!/usr/bin/env python3

import math
from typing import Callable, Dict, List, Set, Tuple, Type, Union

import pytest
from parmed import (
    AngleType,
    AtomType,
    BondType,
    CmapType,
    DihedralType,
    DihedralTypeList,
    ImproperType,
    NonbondedExceptionType,
    UreyBradley,
)
from parmed.charmm import CharmmParameterSet

from chamber import Chamber


# NOTE(MCA): With pytest fixtures, these CharmmParameterSet's are only created
#  once, instead of repeatedly creating them in each test they are used.  This
#  saves a little bit of computational overhead.
@pytest.fixture(scope="module")
def charmm_ff14sb_params():
    return CharmmParameterSet.load_set(
        tfile="tests/data/ff14SB_all.rtf", pfile="tests/data/ff14SB_all.prm"
    )


@pytest.fixture(scope="module")
def chamber_ff14sb_params():
    # Run chamber, creating ff14sb_chamber.str
    chamber = Chamber()
    chamber.write_charmm_parameters()
    return CharmmParameterSet.load_set(sfiles=chamber.charmm_str_name)


@pytest.fixture(scope="module")
def chamber_modrna08_params():
    chamber = Chamber(extra_leaprc_source="modrna08")
    chamber.write_charmm_parameters()
    return CharmmParameterSet.load_set(sfiles=chamber.charmm_str_name)


ParameterType = Union[
    Type[AtomType],
    Type[BondType],
    Type[AngleType],
    Type[UreyBradley],
    Type[DihedralType],
    Type[ImproperType],
    Type[CmapType],
    Type[NonbondedExceptionType],
]


# Mapping of parameter types that are allowed to be missing in the Chamber
# output.  Typically, these represent parameters that are either extraneous or
# were erroneously included in the official CHARMM FF14SB port.
__ALLOW_MISSING: Dict[ParameterType, Set] = {
    AtomType: {},
    BondType: {},
    AngleType: {
        # This non-physical angle is not generated in CHARMM simulations,
        # therefore the associated parameters are unnecessary.  The reason for
        # its inclusion is unclear.  Moreover, the associated spring constant
        # is 0.00, so it only serves to disable the relevant angle if it is
        # erroneously generated.
        ("HW", "HW", "OW"),
        # The authors of the CHARMM FF14SB port renamed "C*" -> "CG".  However,
        # this angle appears to be a holdover from "CG" in an older forcefield
        # as no such "HP-CX-C*" or "C*-CX-HP" angle exists in FF14SB.
        ("CG", "CX", "HP"),
    },
    DihedralType: {},
    ImproperType: {},
}

# Add reversed angles
__ALLOW_MISSING[AngleType] |= {
    tuple(reversed(angle)) for angle in __ALLOW_MISSING[AngleType]
}

ParameterKey = Union[
    str,
    Tuple[str, str],
    Tuple[str, str, str],
    Tuple[str, str, str, str],
]


def __create_atom_type(
    name: str,
    mass: float,
    eps: float,
    rmin: float,
    eps14: float = None,
    rmin14: float = None,
) -> AtomType:
    """
    Factory for parmed.AtomType objects with eps and rmin set.

    :param name: The atom name.
    :param mass: The atom mass.
    :param eps: The atom's epsilon parameter.
    :param rmin: The atom's rmin parameter.
    :param eps14: The atom's epsilon 1-4 parameter. If None, defaults to eps.
    :param rmin14: The atom's rmin 1-4 parameter. If None, defaults to rmin.
    """
    atom_type = AtomType(name=name, number=-1, mass=mass)
    atom_type.set_lj_params(eps=eps, rmin=rmin, eps14=eps14, rmin14=rmin14)
    return atom_type


# Corrections to the official CHARMM FF14SB port that have been manually
# verified.
__VERIFY_CORRECTIONS: Dict[
    ParameterType, List[Tuple[ParameterKey, ParameterType]]
] = {
    AtomType: [
        # These atom types have incorrectly specified parameters in the
        # official CHARMM FF14SB port.  The author even noted that they could
        # not find the new values, and just used the older 1994 values...
        # (AtomType, (eps, rmin))
        # Originally:
        ("CLA", __create_atom_type("CLA", 35.45, -0.035591, 2.513, -0.017795)),
        ("POT", __create_atom_type("POT", 39.1, -0.193683, 1.705, -0.096841)),
        ("SOD", __create_atom_type("SOD", 22.99, -0.087439, 1.369, -0.04372)),
        ("ZN", __create_atom_type("ZN", 65.4, -0.003303, 1.271, -0.001651)),
    ],
    BondType: [],
    AngleType: [
        # In the official CHARMM FF14SB port, each of these have incorrect
        # values. In most cases, the corresponding equilibrium angles are
        # correct and the spring constant is off by less than 1.0.  In other
        # cases, the equilibrium angle is different and the spring constant is
        # off by >> 1.0.  Each of these are due to using antiquated values that
        # have been updated in `leaprc.RNA.OL3`.  The correct values identified
        # by Chamber are available in `parm10.dat`.
        # Originally: AngleType(76.0, 125.1),
        (("CC", "NA", "P"), AngleType(76.7, 125.1)),
        # Originally: AngleType(76.0, 125.1),
        (("CR", "NA", "P"), AngleType(76.7, 125.1)),
        # Originally: parmed.AngleType(42.0, 102.38),
        (("NA", "P", "OP"), AngleType(42.9, 102.38)),
        # Originally: AngleType(1.0, 60.0),
        (("OS", "CT", "OS"), AngleType(160.0, 101.0)),
    ],
    DihedralType: [],
    ImproperType: [],
}

# Mapping of parameter types that are allowed to differ from the official
# FF14SB port.  These are included when the official port contained an error,
# Mapping of parameter types that are allowed to differ from the official
# FF14SB port.  These are included when the official port contained an error,
# usually due to relying on older values no longer utilized in FF14SB.  This
# simply creates sets of the keys of the validated solutions from
# __VERIFY_CORRECTIONS for each parameter type. Ex: {
#     ParameterType.BOND: { ("HW", "OW"), ... },
#     ParameterType.ANGLE: { ("HW", "OW", "HW"), ... },
#     ...
# }
__ALLOW_MISMATCH: Dict[ParameterType, Set[Tuple]] = {
    parameter_type: set(key for key, _ in verify_list)
    for parameter_type, verify_list in __VERIFY_CORRECTIONS.items()
}

__ALLOW_MISMATCH[AngleType] |= {
    tuple(reversed(key)) for key in __ALLOW_MISMATCH[AngleType]
}


def __compare_dihedrals(
    charmm_dihedral: Union[DihedralType, DihedralTypeList],
    chamber_dihedral: Union[DihedralType, DihedralTypeList],
) -> bool:
    """
    Compare two dihedrals for equality.  This allows a relative difference of
    0.02% between attributes to allow for very slight differences in k due to
    rounding errors in the official CHARMM FF14SB port.

    :param charmm_dihedral: The official CHARMM FF14SB dihedral.  This can
     either be a single dihedral, or a list of dihedrals for a full Fourier
     series.
    :param chamber_dihedral: The Chamber-generated dihedral.
    :return: True if the two dihedrals are equal, otherwise False.
    """
    # Convert DihedralType to DihedralTypeList to avoid special cases
    if isinstance(charmm_dihedral, DihedralType):
        charmm_dihedral = DihedralTypeList(charmm_dihedral)
    if isinstance(chamber_dihedral, DihedralType):
        chamber_dihedral = DihedralTypeList(chamber_dihedral)

    # Identify which periodicities are present in which dihedral lists
    charmm_per_to_dihedral = {
        dihedral_type.per: dihedral_type for dihedral_type in charmm_dihedral
    }
    chamber_per_to_dihedral = {
        dihedral_type.per: dihedral_type for dihedral_type in chamber_dihedral
    }
    charmm_periodicities = {per for per in charmm_per_to_dihedral.keys()}
    chamber_periodicities = {per for per in chamber_per_to_dihedral.keys()}

    # Return False, immediately Chamber is missing any periodicities
    if charmm_periodicities - chamber_periodicities:
        return False

    # Return False if any extra periodicities exist with non-zero phi_k
    extra_periodicities = chamber_periodicities - charmm_periodicities
    if extra_periodicities:
        if any(
            chamber_per_to_dihedral[per].phi_k != 0.0
            for per in extra_periodicities
        ):
            return False

    # Generate a zipped list of dihedrals with matching periodicities to
    # compare
    matching_periodicities = charmm_periodicities & chamber_periodicities
    dihedrals_to_compare = zip(
        (charmm_per_to_dihedral[per] for per in matching_periodicities),
        (chamber_per_to_dihedral[per] for per in matching_periodicities),
    )

    # Allow up to a 0.2% mismatch in dihedral parameter attributes.  This is
    # necessary in many cases due to rounding errors in the official CHARMM
    # FF14SB port (e.g., mistakenly rounding up when rounding down was
    # appropriate).
    rel_tol = 2e-3

    # Between the two dihedrals, if any pair of `compare_attrs` aren't within
    # `rel_tol` of one another, return False.
    compare_attrs = ("phi_k", "per", "phase", "scee", "scnb")
    for charmm_dihe, chamber_dihe in dihedrals_to_compare:
        for attr in compare_attrs:
            charmm_dihe_attr = getattr(charmm_dihe, attr)
            chamber_dihe_attr = getattr(chamber_dihe, attr)
            if not math.isclose(
                charmm_dihe_attr, chamber_dihe_attr, rel_tol=rel_tol
            ):
                return False

    # Otherwise, return True
    return True


__COMPARE_EQUAL: Dict[
    ParameterType, Callable[[ParameterType, ParameterType], bool]
] = {
    AtomType: lambda atom1, atom2: atom1 == atom2,
    BondType: lambda bond1, bond2: bond1 == bond2,
    AngleType: lambda angle1, angle2: angle1 == angle2,
    DihedralType: __compare_dihedrals,
    ImproperType: lambda im1, im2: im1 == im2,
}


def compare_types(
    parameter_type: ParameterType,
    charmm_params: CharmmParameterSet,
    chamber_params: CharmmParameterSet,
) -> None:
    """
    Compares the types present between two parameter sets.  Intended for
    testing whether the official CHARMM parameter set
    (:param_ref:`charmm_params`) matches the solution generated by Chamber
    (:param_ref:`chamber_params`).  Fails if the CHARMM parameter set contains
    types not found in the Chamber parameter set, or if any matching types have
    different parameters of the specified parameter type.

    :param parameter_type: The type of parameter to compare, e.g.
     parmed.BondType.
    :param charmm_params:  A parameter set for a force field with the known
     CHARMM solution.
    :param chamber_params:  A parameter set for the same force field generated
     by CHAMBER.
    """
    type_name = parameter_type.__name__.removesuffix("Type")
    type_name_snake = "".join(
        "_" + char.lower() if char.isupper() else char for char in type_name
    )[1:]

    # Ensure all allow_mismatch entries are of the correct ParmEd type
    # Ex: If passing ParameterType.BOND, all mismatch entries should be
    #  of type ParmEd.BondType.
    allow_missing = __ALLOW_MISSING[parameter_type]
    allow_mismatch = __ALLOW_MISMATCH[parameter_type]
    compare_equal = __COMPARE_EQUAL[parameter_type]

    # Process the charmm and chamber type dictionaries
    type_dict_attr = f"{type_name_snake}_types"
    charmm_types = getattr(charmm_params, type_dict_attr)
    chamber_types = getattr(chamber_params, type_dict_attr)
    charmm_type_keys = set(charmm_types.keys())
    chamber_type_keys = set(chamber_types.keys())

    # Check for missing type keys
    n_missing = 0
    missing_type_keys = charmm_type_keys - chamber_type_keys
    for key in missing_type_keys:
        charmm_type = charmm_types[key]
        if key not in allow_missing:
            n_missing += 1
            print(
                f"Chamber missing {type_name} type for key {key}.  "
                f"Official CHARMM FF14SB contains {charmm_type!r} but "
                f"Chamber's output does not."
            )

    if n_missing > 0:
        pytest.fail(f"Found {n_missing} missing {type_name} " f"parameters")

    # Ensure that all types present in the official CHARMM port have matching
    # parameters.
    n_mismatches = 0
    for key in charmm_type_keys:
        if key in missing_type_keys:
            assert (
                key not in chamber_types
            ), f"Expected {key} to be absent from Chamber output"
            continue
        charmm_type = charmm_types[key]
        chamber_type = chamber_types[key]
        if (
            not compare_equal(charmm_type, chamber_type)
            and key not in allow_mismatch
        ):
            n_mismatches += 1
            print(
                f"Found {type_name} type mismatch.  For {key}, official "
                f"CHARMM FF14SB port reports {charmm_type!r}. Chamber gives "
                f"{chamber_type!r}."
            )

    if n_mismatches > 0:
        pytest.fail(f"Found {n_mismatches} mismatched {type_name} parameters")

    # Confirm that all manually-verified corrections were made.
    for key, verified_answer in __VERIFY_CORRECTIONS[parameter_type]:
        assert (
            chamber_types[key] == verified_answer
        ), f"Incorrect {type_name} parameter value for key {key}"

    # Identify new types included with Chamber.  In many cases these are
    # desirable, as they represent the types generated when loading parameters
    # not included in the official CHARMM force field (e.g., modified RNA).
    extra_type_keys = chamber_type_keys - charmm_type_keys
    if extra_type_keys:
        print(
            f"Found {type_name} types in Chamber but not the official CHARMM "
            f"FF14SB port:"
        )
        print(
            ", ".join(
                str(chamber_types[extra_type_key])
                for extra_type_key in extra_type_keys
            )
        )
    pass


def test_atom_types(
    charmm_ff14sb_params, chamber_ff14sb_params, chamber_modrna08_params
):
    compare_types(AtomType, charmm_ff14sb_params, chamber_ff14sb_params)
    compare_types(AtomType, charmm_ff14sb_params, chamber_modrna08_params)


def test_bond_types(
    charmm_ff14sb_params, chamber_ff14sb_params, chamber_modrna08_params
):
    compare_types(BondType, charmm_ff14sb_params, chamber_ff14sb_params)
    compare_types(BondType, charmm_ff14sb_params, chamber_modrna08_params)


def test_angle_types(
    charmm_ff14sb_params, chamber_ff14sb_params, chamber_modrna08_params
):
    compare_types(AngleType, charmm_ff14sb_params, chamber_ff14sb_params)
    compare_types(AngleType, charmm_ff14sb_params, chamber_modrna08_params)


def test_dihedral_types(
    charmm_ff14sb_params, chamber_ff14sb_params, chamber_modrna08_params
):
    compare_types(DihedralType, charmm_ff14sb_params, chamber_ff14sb_params)
    compare_types(DihedralType, charmm_ff14sb_params, chamber_modrna08_params)


def test_improper_types(
    charmm_ff14sb_params, chamber_ff14sb_params, chamber_modrna08_params
):
    compare_types(ImproperType, charmm_ff14sb_params, chamber_ff14sb_params)
    compare_types(ImproperType, charmm_ff14sb_params, chamber_modrna08_params)


# TODO(MCA): Ensure other types are empty for all three fixtures
