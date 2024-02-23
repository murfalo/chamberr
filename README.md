# chamber: Convert Amber to CHARMM Force Fields

This repository uses [parmed](https://github.com/parmed/parmed) to convert Amber forcefields to CHARMM-readable formats (.prm, .rtf, and .str).

## Getting Started

Install dependencies to a `conda` environment, activate said environment, then run:

``` bash
> pip install -e .
> chamber -h
```

## Notes

- Currently, protein.ff14SB and RNA.OL3 are used as a "baseline".  Feel free to expand this to suit the needs of your desired forcefield.

## Limitations

- ParmEd does not support exporting full topologies.  You will need to create your own.  In most cases, you can steal existing CHARMM topologies for your molecules and convert their atom types to those used by Amber.
- For unclear reasons, some types (CU, FE, and Se) remain poorly defined even after loading these baseline forcefields.  Their definitions may be found elsewhere, but I do not use them so I simply discard them.  Feel free to fix this if you are in need!

## Dependencies

- AmberTools
- ParmEd


