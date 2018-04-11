"""AMPAL objects that represent ligands."""

from ampal.base_ampal import Polymer, Monomer


class LigandGroup(Polymer):
    """A container for `Ligand` `Monomers`.

    Parameters
    ----------
    monomers : Monomer or [Monomer], optional
        Monomer or list containing Monomer objects to form the Polymer().
    polymer_id : str, optional
        An ID that the user can use to identify the `Polymer`. This is
        used when generating a pdb file using `Polymer().pdb`.
    parent : ampal.Assembly, optional
        Reference to `Assembly` containing the `Polymer`.
    sl : int, optional
        The default smoothing level used when calculating the
        backbone primitive.
    """

    def __init__(self, monomers=None, polymer_id=' ', parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id, molecule_type='ligands',
            parent=parent, sl=sl)

    def __repr__(self):
        return '<Ligands chain containing {} {}>'.format(
            len(self._monomers),
            'Ligand' if len(self._monomers) == 1 else 'Ligands')


class Ligand(Monomer):
    """`Monomer` that represents a `Ligand`.

    Notes
    -----
    All `Monomers` that do not have dedicated classes are
    represented using the `Ligand` class.

    Parameters
    ----------
    mol_code : str
        PDB molecule code that represents the monomer.
    atoms : OrderedDict, optional
        OrderedDict containing Atoms for the Monomer. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    monomer_id : str, optional
        String used to identify the residue.
    insertion_code : str, optional
        Insertion code of monomer, used if reading from pdb.
    is_hetero : bool, optional
        True if is a hetero atom in pdb. Helps with PDB formatting.

    Attributes
    ----------
    atoms : OrderedDict
        OrderedDict containing Atoms for the Monomer. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    mol_code : str
        PDB molecule code that represents the `Ligand`.
    insertion_code : str
        Insertion code of `Ligand`, used if reading from pdb.
    is_hetero : bool
        True if is a hetero atom in pdb. Helps with PDB formatting.
    self.states : dict
        Contains an `OrderedDicts` containing atom information for each
        state available for the `Ligand`.
    id : str
        String used to identify the residue.
    parent : Polymer or None
        A reference to the `LigandGroup` containing this `Ligand`.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.
    """

    def __init__(self, mol_code, atoms=None, monomer_id=' ', insertion_code=' ',
                 is_hetero=False, parent=None):
        super(Ligand, self).__init__(
            atoms, monomer_id, parent=parent)
        self.mol_code = mol_code
        self.insertion_code = insertion_code
        self.is_hetero = is_hetero

    def __repr__(self):
        return '<Ligand containing {} {}. Ligand code: {}>'.format(
            len(self.atoms), 'Atom' if len(self.atoms) == 1 else 'Atoms',
            self.mol_code)


__author__ = "Christopher W. Wood, Kieran L. Hudson"
