# Chamber: convert Amber force fields for CHARMM ![CI Status](https://github.com/murfalo/chamber/actions/workflows/ci.yaml/badge.svg)

Chamber uses [ParmEd](https://github.com/parmed/parmed) to convert Amber force fields to CHARMM-readable `.str` files.

## Getting Started

```bash
> mamba env create -f environment.yaml
> mamba activate chamber
> pip install -e .
```

At time of writing, a few custom ParmEd bugfixes are required for Chamber to work properly.  To install my fork of ParmEd:

```bash
> git clone https://github.com/murfalo/parmed
> cd parmed
> bash -ex devtools/ci/install.sh
> chamber -h
```

## Examples

```bash
> chamber  # generate CHARMM-readable FF14SB force field
> chamber modrna08  # generate CHARMM-readable FF14SB + modrna08 force field
> chamber <some_other_force_field>  # requires care, not thoroughly tested
```

You can now utilize the generated `chamber_ff14sb.str` and `chamber_modrna08.str` files in CHARMM!

## Resources

Chamber makes extensive use of [ParmEd](https://github.com/parmed/parmed) and [Amber](https://ambermd.org/index.php).  The respective home pages are great starting points, along with the [Amber reference manual](https://ambermd.org/doc12/Amber23.pdf).

The force fields themselves are located at `$AMBERHOME/dat/leap/cmd/leaprc.*`.  Parameters are most often found in `$AMBERHOME/dat/leap/parm/<frcmod.*|parm*.dat>`.  To know which files to look at, it's often easiest to start from the respective `leaprc` and then use `find $AMBERHOME -type f -name '<your_file_name>'` to locate any files that pique your interest.

## Testing

For FF14SB, Chamber has been extensively tested to reproduce the officially distributed CHARMM FF14SB force field (with the exception of a few corrections discovered via Chamber).  For modrna08, Chamber exactly reproduces all the parameters for FF14SB along with including additional parameters for representing terms with modified RNA atom types.

These tests are managed automatically through [GitHub actions](https://github.com/murfalo/chamber/actions/), but you can run them yourself with:

```bash
> pytest
```

## Limitations

- Chamber has only been extensively tested for converting the Amber FF14SB (protein.ff14sb, RNA.OL3, water.tip3p) forcefield +/- the Amber modified RNA forcefield (modrna08).
- Other protein, RNA, and solvent force field sources can be easily specified but are not offered as commandline arguments.
- ParmEd does not support exporting full topologies.  You will need to create your own.  For Amber FF14SB, topologies are distributed with CHARMM's source code in `toppar/non_charmm/parm14sb_all.rtf`.
