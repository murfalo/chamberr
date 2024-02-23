# chamber: Convert Amber to CHARMM Force Fields

This repository uses [parmed](https://github.com/parmed/parmed) to convert Amber forcefields to CHARMM-readable formats (.prm, .rtf, and .str).

## Getting Started

Install dependencies to a `conda` environment, activate said environment, then run:

``` bash
> mamba env create -f environment.yaml
> mamba activate chamber
> pip install -e .
> chamber -h
```


## Examples

To get the forcefield terms associated with new atom types in the `modrna08` Amber forcefield:

```bash
> chamber --only_unique modrna08
```

You can now utilize the generated `modrna08.prm`, `modrna08.rtf`, and/or `modrna08.str` files to extend the officially distributed `parm14sb_all.prm` parameters and `parm14sb_all.rtf` topologies in your CHARMM scripts!

## Notes

- Currently, protein.ff14SB and RNA.OL3 are used as a "baseline".  Feel free to expand this to suit the needs of your desired forcefield.

## Limitations

- ParmEd does not support exporting full topologies.  You will need to create your own.  In most cases, you can steal existing CHARMM topologies for your molecules and convert their atom types to those used by Amber.
- For unclear reasons, some types (CU, FE, and Se) remain poorly defined even after loading these baseline forcefields.  Their definitions may be found elsewhere, but I do not use them so I simply discard them.  Feel free to fix this if you are in need!

