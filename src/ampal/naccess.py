"""This module provides an interface to the program NACCESS."""

import subprocess
import tempfile
import os


def naccess_available():
    """True if naccess is available on the path."""
    available = False
    try:
        subprocess.check_output(["naccess"], stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        available = True
    except FileNotFoundError:
        print(
            "naccess has not been found on your path. If you have already "
            "installed naccess but are unsure how to add it to your path, "
            "check out this: https://stackoverflow.com/a/14638025"
        )
    return available


def run_naccess(
    pdb, mode, path=True, include_hetatms=False, outfile=None, path_to_ex=None
):
    """Uses naccess to run surface accessibility calculations.

    Notes
    -----
    Requires the naccess program, with a path to its executable
    provided in global_settings. For information on the Naccess program,
    see: http://www.bioinf.manchester.ac.uk/naccess/
    This includes information on the licensing, which is not free for
    Industrial and Profit-making instituions.

    Parameters
    ----------
    pdb : str
        Path to pdb file or string.
    mode : str
        Return mode of naccess. One of 'asa', 'rsa' or 'log'.
    path : bool, optional
        Indicates if pdb is a path or a string.
    outfile : str, optional
        Filepath for storing the naccess output.
    path_to_ex : str or None
        Path to the binary for naccess, if none then it is assumed
        that the binary is available on the path as `naccess`.

    Returns
    -------
    naccess_out : str
        naccess output file for given mode as a string.
    """
    if mode not in ["asa", "rsa", "log"]:
        raise ValueError(
            "mode {} not valid. Must be 'asa', 'rsa' or 'log'".format(mode)
        )
    if path_to_ex:
        naccess_exe = path_to_ex
    else:
        naccess_exe = "naccess"

    if not path:
        if type(pdb) == str:
            pdb = pdb.encode()
    else:
        with open(pdb, "r") as inf:
            pdb = inf.read().encode()

    this_dir = os.getcwd()
    # temp pdb file in temp dir.
    temp_dir = tempfile.TemporaryDirectory()
    temp_pdb = tempfile.NamedTemporaryFile(dir=temp_dir.name)
    temp_pdb.write(pdb)
    temp_pdb.seek(0)
    # run naccess in the temp_dir. Files created will be written here.
    os.chdir(temp_dir.name)

    if include_hetatms:
        naccess_args = "-h"
        subprocess.check_output([naccess_exe, naccess_args, temp_pdb.name])
    else:
        subprocess.check_output([naccess_exe, temp_pdb.name])
    temp_pdb.close()
    with open(".{}".format(mode), "r") as inf:
        naccess_out = inf.read()
    # navigate back to initial directory and clean up.
    os.chdir(this_dir)
    if outfile:
        with open(outfile, "w") as inf:
            inf.write(naccess_out)
    temp_dir.cleanup()

    return naccess_out


def total_accessibility(in_rsa, path=True):
    """Parses rsa file for the total surface accessibility data.

    Parameters
    ----------
    in_rsa : str
        Path to naccess rsa file.
    path : bool
        Indicates if in_rsa is a path or a string.

    Returns
    -------
    dssp_residues : 5-tuple(float)
        Total accessibility values for:
        [0] all atoms
        [1] all side-chain atoms
        [2] all main-chain atoms
        [3] all non-polar atoms
        [4] all polar atoms

    """
    if path:
        with open(in_rsa, "r") as inf:
            rsa = inf.read()
    else:
        rsa = in_rsa[:]
    all_atoms, side_chains, main_chain, non_polar, polar = [
        float(x) for x in rsa.splitlines()[-1].split()[1:]
    ]
    return all_atoms, side_chains, main_chain, non_polar, polar


def extract_residue_accessibility(in_rsa, path=True, get_total=False):
    """Parses rsa file for solvent accessibility for each residue.

    Parameters
    ----------
    in_rsa : str
        Path to naccess rsa file
    path : bool
        Indicates if in_rsa is a path or a string
    get_total : bool
        Indicates if the total accessibility from the file needs to
        be extracted. Convenience method for running the
        total_accessibility function but only running NACCESS once

    Returns
    -------
    rel_solv_ac_acc_atoms : list
        Relative solvent accessibility of all atoms in each amino acid
    get_total : float
        Relative solvent accessibility of all atoms in the NACCESS rsa file
    """

    if path:
        with open(in_rsa, "r") as inf:
            rsa = inf.read()
    else:
        rsa = in_rsa[:]

    residue_list = [x for x in rsa.splitlines()]
    rel_solv_acc_all_atoms = [
        float(x[22:28]) for x in residue_list if x[0:3] == "RES" or x[0:3] == "HEM"
    ]

    if get_total:
        (all_atoms, _, _, _, _) = total_accessibility(rsa, path=False)
        return rel_solv_acc_all_atoms, all_atoms
    return rel_solv_acc_all_atoms, None


__author__ = "Jack W. Heal, Gail J. Bartlett"
