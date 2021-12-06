"""This module provides an interface to the program DSSP.

For more information on DSSP see [4]_.

References
----------
.. [4] Kabsch W, Sander C (1983) "Dictionary of protein
   secondary structure: pattern recognition of hydrogen-bonded
   and geometrical features", Biopolymers, 22, 2577-637.
"""

import subprocess
import tempfile

from .assembly import Assembly


def dssp_available():
    """True if mkdssp is available on the path."""
    available = False
    try:
        subprocess.check_output(["mkdssp"], stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        available = True
    except FileNotFoundError:
        print(
            "DSSP has not been found on your path. If you have already "
            "installed DSSP but are unsure how to add it to your path, "
            "check out this: https://stackoverflow.com/a/14638025"
        )
    return available


def run_dssp(pdb, path=True):
    """Uses DSSP to find helices and extracts helices from a pdb file or string.
    Parameters
    ----------
    pdb : str
        Path to pdb file or string.
    path : bool, optional
        Indicates if pdb is a path or a string.

    Returns
    -------
    dssp_out : str
        Std out from DSSP.
    """
    if not path:
        if isinstance(pdb, str):
            pdb = pdb.encode()
        with tempfile.NamedTemporaryFile() as temp_pdb:
            temp_pdb.write(pdb)
            temp_pdb.seek(0)
            dssp_out = subprocess.check_output(["mkdssp", temp_pdb.name])
    else:
        dssp_out = subprocess.check_output(["mkdssp", pdb])
    dssp_out = dssp_out.decode()
    return dssp_out


def extract_all_ss_dssp(in_dssp, path=True):
    """Uses DSSP to extract secondary structure information on every residue.

    Parameters
    ----------
    in_dssp : str
        Path to DSSP file.
    path : bool, optional
        Indicates if pdb is a path or a string.

    Returns
    -------
    dssp_residues : [tuple]
        Each internal list contains:
            [0] int Residue number
            [1] str Secondary structure type
            [2] str Chain identifier
            [3] str Residue type
            [4] float Phi torsion angle
            [5] float Psi torsion angle
            [6] int dssp solvent accessibility
    """

    if path:
        with open(in_dssp, "r") as inf:
            dssp_out = inf.read()
    else:
        dssp_out = in_dssp[:]
    dssp_residues = []
    active = False
    for line in dssp_out.splitlines():
        if active:
            try:
                res_num = int(line[5:10].strip())
                chain = line[10:12].strip()
                residue = line[13]
                ss_type = line[16]
                phi = float(line[103:109].strip())
                psi = float(line[109:116].strip())
                acc = int(line[35:38].strip())
                dssp_residues.append((res_num, ss_type, chain, residue, phi, psi, acc))
            except ValueError:
                pass
        else:
            if line[2] == "#":
                active = True
    return dssp_residues


def find_ss_regions(dssp_residues, loop_assignments=(" ", "B", "S", "T")):
    """Separates parsed DSSP data into groups of secondary structure.

    Notes
    -----
    Example: all residues in a single helix/loop/strand will be gathered
    into a list, then the next secondary structure element will be
    gathered into a separate list, and so on.

    Parameters
    ----------
    dssp_residues : [tuple]
        Each internal list contains:
            [0] int Residue number
            [1] str Secondary structure type
            [2] str Chain identifier
            [3] str Residue type
            [4] float Phi torsion angle
            [5] float Psi torsion angle
            [6] int dssp solvent accessibility

    Returns
    -------
    fragments : [[list]]
        Lists grouped in continuous regions of secondary structure.
        Innermost list has the same format as above.
    """

    loops = loop_assignments
    previous_ele = None
    fragment = []
    fragments = []
    for ele in dssp_residues:
        if previous_ele is None:
            fragment.append(ele)
        elif ele[2] != previous_ele[2]:
            fragments.append(fragment)
            fragment = [ele]
        elif previous_ele[1] in loops:
            if ele[1] in loops:
                fragment.append(ele)
            else:
                fragments.append(fragment)
                fragment = [ele]
        else:
            if ele[1] == previous_ele[1]:
                fragment.append(ele)
            else:
                fragments.append(fragment)
                fragment = [ele]
        previous_ele = ele
    fragments.append(fragment)
    return fragments


def tag_dssp_data(assembly, loop_assignments=(" ", "B", "S", "T")):
    """Adds output data from DSSP to an Assembly.

    A dictionary will be added to the `tags` dictionary of each
    residue called `dssp_data`, which contains the secondary
    structure definition, solvent accessibility phi and psi values
    from DSSP. A list of regions of continuous secondary assignments
    will also be added to each `Polypeptide`.

    The tags are added in place, so nothing is returned from this
    function.

    Parameters
    ----------
    assembly : ampal.Assembly
        An Assembly containing some protein.
    loop_assignments : tuple or list
        A tuple containing the DSSP secondary structure identifiers to
        that are classed as loop regions.
    """
    dssp_out = run_dssp(assembly.pdb, path=False)
    dssp_data = extract_all_ss_dssp(dssp_out, path=False)
    for record in dssp_data:
        rnum, sstype, chid, _, phi, psi, sacc = record
        assembly[chid][str(rnum)].tags["dssp_data"] = {
            "ss_definition": sstype,
            "solvent_accessibility": sacc,
            "phi": phi,
            "psi": psi,
        }
    ss_regions = find_ss_regions(dssp_data, loop_assignments)
    for region in ss_regions:
        chain = region[0][2]
        ss_type = " " if region[0][1] in loop_assignments else region[0][1]
        first_residue = str(region[0][0])
        last_residue = str(region[-1][0])
        if not "ss_regions" in assembly[chain].tags:
            assembly[chain].tags["ss_regions"] = []
        assembly[chain].tags["ss_regions"].append(
            (first_residue, last_residue, ss_type)
        )
    return


def get_ss_regions(assembly, ss_types):
    """Returns an Assembly containing Polymers for each region of structure.

    Parameters
    ----------
    assembly : ampal.Assembly
        `Assembly` object to be searched secondary structure regions.
    ss_types : list
        List of secondary structure tags to be separate i.e. ['H']
        would return helices, ['H', 'E'] would return helices
        and strands.

    Returns
    -------
    fragments : Assembly
        `Assembly` containing a `Polymer` for each region of specified
        secondary structure.
    """
    if not any(map(lambda x: "ss_regions" in x.tags, assembly)):
        raise ValueError(
            "This assembly does not have any tagged secondary structure "
            "regions. Use `ampal.dssp.tag_dssp_data` to add the tags."
        )
    fragments = Assembly()
    for polypeptide in assembly:
        if "ss_regions" in polypeptide.tags:
            for start, end, ss_type in polypeptide.tags["ss_regions"]:
                if ss_type in ss_types:
                    fragment = polypeptide.get_slice_from_res_id(start, end)
                    fragments.append(fragment)
    if not fragments:
        raise ValueError(
            "No regions matching that secondary structure type"
            " have been found. Use standard DSSP labels."
        )
    return fragments


__author__ = "Christopher W. Wood, Gail J. Bartlett"
