import unittest
import itertools
import pathlib
import random
import numpy

import ampal

TEST_FILE_FOLDER = pathlib.Path(__file__).parent / 'testing_files'


class PolypeptideStructuralValidationTestCase(unittest.TestCase):
    """ Tests for the valid_backbone_distance method of Polypeptide class. """

    def setUp(self):
        test_files = [str(TEST_FILE_FOLDER / x) for x in [
            '1ek9.pdb', '2ht0.pdb', '3qy1.pdb']]
        test_structures = [ampal.load_pdb(x) for x in test_files]
        self.test_polypeptides = [
            p for p in itertools.chain(*test_structures)
            if isinstance(p, ampal.Polypeptide)]

    def test_valid_distances_testing_files(self):
        # checks for valid bond distances in testing files.
        valid_distances = [p.valid_backbone_bond_lengths()
                           for p in self.test_polypeptides]
        self.assertTrue(all(valid_distances))

    def test_invalid_distance_natural(self):
        polypeptides = self.test_polypeptides[:]
        # translate random monomer
        for p in polypeptides:
            p[random.choice(range(1, len(p)))].translate([0, 0, 10])
            self.assertFalse(p.valid_backbone_bond_lengths())

    def test_valid_angles_testing_files(self):
        # checks for valid bond angles in testing files.
        valid_angles = [p.valid_backbone_bond_angles()
                        for p in self.test_polypeptides]
        self.assertTrue(all(valid_angles))


__author__ = 'Jack W. Heal'
