import pathlib

import ampal
from ampal.interactions import find_salt_bridges

TEST_FILE_FOLDER = pathlib.Path(__file__).parent / "testing_files"


def test_find_salt_bridges():
    test_ampal = ampal.load_pdb(str(TEST_FILE_FOLDER / "3qy1.pdb"))
    salt_bridges = find_salt_bridges(test_ampal.get_atoms())
    # This salt bridges in this structure were manually validated
    assert len(salt_bridges) == 29
    salt_bridges = find_salt_bridges(test_ampal.get_atoms(), min_dist=1, max_dist=1.5)
    assert len(salt_bridges) == 0
    salt_bridges = find_salt_bridges(test_ampal.get_atoms(), positive_labels={})
    assert len(salt_bridges) == 0
    salt_bridges = find_salt_bridges(test_ampal.get_atoms(), negative_labels={})
    assert len(salt_bridges) == 0
    return
