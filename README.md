# Chamber: convert Amber force fields for CHARMM ![CI Status](https://github.com/murfalo/chamber/actions/workflows/ci.yaml/badge.svg)

Chamber uses [ParmEd](https://github.com/parmed/parmed) to convert Amber force fields to CHARMM-readable `.str` files.

## Getting Started

Several converted force fields are already available in `ff_stream`, along with the Amber FF14SB parameters and topologies distributed with CHARMM (see `parm14sb_all.{prm|rtf}`).  To convert some yourself, you will need to setup up Chamber with its required dependencies:

```bash
> mamba env create -f environment.yaml
> mamba activate chamber
> pip install -e .
```

Currently, a few custom ParmEd bugfixes are required for Chamber to work properly.  To install my fork of ParmEd and run `chamber`:

```bash
> git clone https://github.com/murfalo/parmed
> cd parmed
> bash -ex devtools/ci/install.sh
> chamber -h
```

## Examples

Only FF14SB and FF14SB + modrna08 have been thoroughly tested.  In all other cases, proceed with care.

```bash
> chamber  # FF14SB force field
> chamber modrna08  # FF14SB + modrna08 force field
> chamber -r Shaw -s tip4pd  # D.E. Shaw RNA + TIP4PD
> chamber -r Shaw -s opc modrna08  # D.E. Shaw RNA + OPC + modrna08
```

You can now utilize the generated `chamber_none_OL3_ff14sb_tip3p.str` (FF14SB) and `chamber_modrna08_OL3_ff14sb_tip3p.str` (FF14SB + modrna08) files in CHARMM!

## Resources

Chamber makes extensive use of [ParmEd](https://github.com/parmed/parmed) and [Amber](https://ambermd.org/index.php).  The respective home pages are great starting points, along with the [Amber reference manual](https://ambermd.org/doc12/Amber23.pdf).

The force fields themselves are located at `$AMBERHOME/dat/leap/cmd/leaprc.*`.  Parameters are most often found in `$AMBERHOME/dat/leap/parm/<frcmod.*|parm*.dat>`.  To know which files to look at, it's often easiest to start from the respective `leaprc` and then use `find $AMBERHOME -type f -name '<your_file_name>'` to locate any files that pique your interest.

## Testing

Chamber's tests ensure that it exactly reproduces the officially distributed CHARMM FF14SB force field, excepting a few manually verified corrections.  The included forcefields in this distribution are protein.ff14SB, RNA.OL3, and water.tip3p in descending priority.  When including modrna08, Chamber's tests ensure that these parameters are maintained.  Thus, the generated forcefield is confirmed to be a superset of FF14SB with the new modrna08 parameters.

These tests are managed automatically through [GitHub actions](https://github.com/murfalo/chamber/actions/), but you can run them yourself with:

```bash
> pytest
```

## Limitations

- ParmEd does not support exporting full topologies.  You will need to find a source or create your own.  For Amber FF14SB, topologies are distributed with CHARMM's source code in `toppar/non_charmm/parm14sb_all.rtf`.

