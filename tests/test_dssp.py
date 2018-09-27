"""Tests basic DSSP functionality."""

import pathlib
import unittest

from ampal import Residue, load_pdb, tag_dssp_data
from ampal.dssp import dssp_available, get_ss_regions


TEST_FILE_FOLDER = pathlib.Path(__file__).parent / 'testing_files'


@unittest.skipUnless(dssp_available(), "External program not detected.")
class TestDSSP(unittest.TestCase):
    """Tests basic functionality of dssp.py."""

    def setUp(self):
        test_file_path = str(TEST_FILE_FOLDER / '3qy1.pdb')
        self.structure = load_pdb(test_file_path)
        tag_dssp_data(self.structure)

    def test_dssp_tagging(self):
        """Test DSSP tagging of structure."""
        for monomer in self.structure.get_monomers():
            if isinstance(monomer, Residue):
                self.assertTrue('dssp_data' in monomer.tags)
            else:
                self.assertFalse('dssp_data' in monomer.tags)

    def test_get_all_ss_regions(self):
        """Test secondary structure extraction for all labels.

        Note
        ----
        This should return every residue in the `Assembly`.
        """
        ss_regions = get_ss_regions(
            self.structure, ['H', 'B', 'E', 'G', 'I', 'T', 'S', ' '])
        self.assertTrue(''.join(self.structure.sequences),
                        ''.join(ss_regions.sequences))

    def test_get_helices(self):
        """Test helix extraction with `get_ss_regions`."""
        helices = get_ss_regions(self.structure, ['H'])
        self.assertTrue(len(helices) == 20)

    def test_get_strands(self):
        """Test strand extraction with `get_ss_regions`."""
        strands = get_ss_regions(self.structure, ['E'])
        self.assertTrue(len(strands) == 10)
