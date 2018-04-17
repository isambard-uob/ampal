# AMPAL
A simple, intuitive and Pythonic framework for representing biomolecular structure.

[![CircleCI](https://circleci.com/gh/isambard-uob/ampal/tree/master.svg?style=shield)](https://circleci.com/gh/isambard-uob/ampal/tree/master)
[![Python Version](https://img.shields.io/badge/python-3.5%2C%203.6-lightgrey.svg)]()
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/isambard-uob/ampal/blob/master/LICENSE)

## Installation

You can install AMPAL from pip:

`pip install ampal`

Or from source by downloading/cloning this repository, navigating to the folder
and typing:

`pip install .`

AMPAL uses Cython, so if you're installing from source make sure you have it
installed.


## Super Quick Start

Load a PDB file into AMPAL:

```Python
my_structure = ampal.load_pdb('3qy1.pdb')
print(my_structure)
# OUT: <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>
```

Select regions of the structure in an intuitive manner:

```Python
my_atom = my_structure['A']['56']['CA']
print(my_structure['A']['56']['CA'])
# OUT: <Carbon Atom (CA). Coordinates: (6.102, -4.287, -29.607)>
```

Then climb all the way back up the hierachy:

```Python
print(my_atom.parent)
# OUT: <Residue containing 9 Atoms. Residue code: GLU>
print(my_atom.parent.parent)
# OUT: <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>
print(my_atom.parent.parent.parent)
# OUT: <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>
```

This is just a quick introduction, AMPAL contain tonnes of tools for making
complex selections and performing analysis. Take a look at the
[docs](https://isambard-uob.github.io/ampal/) to find out more.
