# AMPAL
A simple, intuitive and Pythonic framework for representing biomolecular structure.

## Installation

AMPAL is currently tested with version of Python above 3.8, although it should
still work with earlier versions. You can install AMPAL from pip:

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

## Release Notes

## v1.5.0

* Adds rotamer classification with the `classify_angle_as_rotamer` in
  the `analyse_protein` module. 
* Fixes a bug with parsing structure files

## v1.4.0

* **Adds `get_ss_regions` to `ampal.dssp`.** This function can be used to
  extract all regions of a protein in a particular secondary structure.
* **Fixes bug with DSSP `ss_region` tagging.** End residues used to be missed.

## v1.3.0

* **Adds an interface for NACCESS.** Functions for using NACCESS to calculate
  solvent accessibility.

### v1.2.0

* **Adds an interface for DSSP.** If you have DSSP on your computer and have the
  `mkdssp` command available on your path, you can use the `ampal.tag_dssp_data`
  function to add secondary structure information to the tags dictionary of the
  residues in your structure.
* **Adds the `ampal.align` module.** Contains a simple class for aligning two
  `Polypeptides` using MMC. The simplest interface is the `align_backbones`
  function.
  * This is currently super inefficient and will be reimplemented.

### v1.1.0

* **Adds the centroid property to residues.**
