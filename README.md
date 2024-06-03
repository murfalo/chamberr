# Chamberr: convert Amber force fields for CHARMM ![CI Status](https://github.com/murfalo/chamberr/actions/workflows/ci.yaml/badge.svg)

The reverse of the popular Amber tool `chamber`, Chamberr uses [ParmEd](https://github.com/parmed/parmed) to convert Amber force fields to CHARMM-readable `.str` files.  Chamberr was originally built to port Amber's FF14SB, but should work for other force fields as well.

## Getting Started

Several converted force fields are already available in `ff_stream`, along with the Amber FF14SB parameters and topologies distributed with CHARMM (see `parm14sb_all.{prm|rtf}`).  To convert some yourself, you will need to setup up Chamberr with its required dependencies:

> [!WARNING]
> The officially-distributed Amber FF14SB files (`parm14sb_all.{prm|rtf}`) contain several minor errors.  For more info, see "`Using parm14sb_all`" below.

```bash
> mamba env create -f environment.yaml
> mamba activate chamberr
> pip install -e .
```

Currently, a few custom ParmEd bugfixes are required for Chamberr to work properly.  To install my fork of ParmEd and run `chamberr`:

```bash
> git clone https://github.com/murfalo/parmed
> cd parmed
> # This may result in `pytest` errors, but should successfully install ParmEd
> bash -ex devtools/ci/install.sh
> # Print usage information to ensure `chamberr` is working properly
> chamberr -h
```

## Examples

Only FF14SB and FF14SB + modrna08 have been thoroughly tested.  In all other cases, proceed with care.

```bash
> chamberr  # FF14SB force field
> chamberr modrna08  # FF14SB + modrna08 force field
> chamberr -r Shaw -s tip4pd  # D.E. Shaw RNA + TIP4PD
> chamberr -r Shaw -s opc modrna08  # D.E. Shaw RNA + OPC + modrna08
```

You can now utilize the generated `chamberr_none_OL3_ff14sb_tip3p.str` (FF14SB) and `chamberr_modrna08_OL3_ff14sb_tip3p.str` (FF14SB + modrna08) files in CHARMM!

## Using `parm14sb_all`

`parm14sb_all` is the name for the FF14SB parameters and topologies officially distributed with [CHARMM](https://charmm.chemistry.harvard.edu/request_license.php?version=free).  While developing Chamberr, several issues were noted:

1. `parm14sb_all.prm`:
    - Incorrectly included angle parameter for CG-CX-HP
        - When creating `parm14sb_all.prm` the authors renamed atom type C* to CG.  This parameter appears to be for the original CG (not C\*) which does not exist in FF14SB.
    - Select ions use outdated 1994 Lennard-Jones parameters (CLA, POT, SOD, ZN).  The authors directly noted that they could not locate the newer parameters.
    - Several angle and dihedral parameters are off by ~<1% due to rounding errors.
2. `parm14sb_all.rtf`:
    - `RESI GLU` was missing its `GROUP` definition
    - `RESI`s ADE, GUA, and TRP were missing several improper definitions
    - `RESI`s ADE, GUA, and URA had several misordered improper definitions

For these reasons, using the `Chamberr`-generated parameters may give slight improvements over `parm14sb_all.prm`.  Chamberr does not yet support automatic generation of topologies, so a corrected version of `parm14sb_all.rtf` has been made available in `ff_stream`.  Rudimentary scripts for automatic topology generation are available with [msld-PUMA](https://github.com/murfalo/msld-PUMA/) (currently private, please request access or contact git@murfalo.com).

## Resources

Chamberr makes extensive use of [ParmEd](https://github.com/parmed/parmed) and [Amber](https://ambermd.org/index.php).  The respective home pages are great starting points, along with the [Amber reference manual](https://ambermd.org/doc12/Amber23.pdf).

The force fields themselves are located at `$AMBERHOME/dat/leap/cmd/leaprc.*`.  Parameters are most often found in `$AMBERHOME/dat/leap/parm/<frcmod.*|parm*.dat>`.  To know which files to look at, it's often easiest to start from the respective `leaprc` and then use `find $AMBERHOME -type f -name '<your_file_name>'` to locate any files that pique your interest.

## Testing

Chamberr's tests ensure that it exactly reproduces the officially distributed CHARMM FF14SB force field, excepting a few manually verified corrections.  The included forcefields in this distribution are protein.ff14SB, RNA.OL3, and water.tip3p in descending priority.  When including modrna08, Chamberr's tests ensure that these parameters are maintained.  Thus, the generated forcefield is confirmed to be a superset of FF14SB with the new modrna08 parameters.

These tests are managed automatically through [GitHub actions](https://github.com/murfalo/chamberr/actions/), but you can run them yourself with:

```bash
> pytest
```

## Limitations

- ParmEd does not support exporting full topologies.  You will need to find a source or create your own.  For Amber FF14SB, topologies are distributed with CHARMM's source code in `toppar/non_charmm/parm14sb_all.rtf`.

