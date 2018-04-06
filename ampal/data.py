"""Loads data required for running ISAMBARD."""

import json
import pathlib


CE_PATH = pathlib.Path(__file__).parent / 'datafiles' / \
    'chemical_elements.json'
COL_PATH = pathlib.Path(__file__).parent / 'datafiles' / \
    'pdb_atom_column_format.json'

with open(str(CE_PATH), 'r') as inf:
    ELEMENT_DATA = json.loads(inf.read())


with open(str(COL_PATH), 'r') as inf:
    PDB_ATOM_COLUMNS = json.loads(inf.read())


__author__ = "Christopher W. Wood"
