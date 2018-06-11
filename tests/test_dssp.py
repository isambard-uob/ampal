"""Tests basic DSSP functionality."""

import pathlib
import unittest

from ampal import Residue, load_pdb, tag_dssp_data
from ampal.dssp import dssp_available


TEST_FILE_FOLDER = pathlib.Path(__file__).parent / 'testing_files'


@unittest.skipUnless(dssp_available(), "External program not detected.")
class TestDSSP(unittest.TestCase):
    """Tests basic functionality of dssp.py."""

    def test_dssp_tagging(self):
        """Test DSSP tagging of structure."""
        test_file_path = str(TEST_FILE_FOLDER / '3qy1.pdb')
        structure = load_pdb(test_file_path)
        tag_dssp_data(structure)
        for monomer in structure.get_monomers():
            if isinstance(monomer, Residue):
                self.assertTrue('dssp_data' in monomer.tags)
            else:
                self.assertFalse('dssp_data' in monomer.tags)
