# g16-isomers

This set of Python scripts will automate calculations of molecules using
Gaussian '16 on a computing cluster. The primary purpose of these scripts is to
perform mass calculations and compare structures and properties across isomers.

There are a few modes of operation, but primarily revolves around
setting up and submitting calculation jobs to a scheduler system
from `xyz` files contained in a `structures` folder.

These `xyz` files are generated primarily from SMILES format using
`obabel`. 

## Requirements

`g16`

