"""AMPAL objects that represent protein."""

from collections import OrderedDict
import warnings

import numpy

from ampal.base_ampal import Polymer, Monomer, Atom
from ampal.pseudo_atoms import Primitive
from ampal.analyse_protein import (
    make_primitive_extrapolate_ends,
    measure_torsion_angles,
    residues_per_turn,
    polymer_to_reference_axis_distances,
    crick_angles,
    alpha_angles,
    sequence_molecular_weight,
    sequence_molar_extinction_280,
    sequence_isoelectric_point,
    measure_sidechain_torsion_angles,
)
from ampal.interactions import (
    generate_covalent_bond_graph,
    generate_bond_subgraphs_from_break,
    find_covalent_bonds,
)
from .amino_acids import (
    get_aa_code,
    get_aa_letter,
    ideal_backbone_bond_lengths,
    ideal_backbone_bond_angles,
)
from .geometry import (
    Quaternion,
    unit_vector,
    dihedral,
    find_transformations,
    distance,
    angle_between_vectors,
)
from .ampal_warnings import MalformedPDBWarning


def flat_list_to_polymer(atom_list, atom_group_s=4):
    """Takes a flat list of atomic coordinates and converts it to a `Polymer`.

    Parameters
    ----------
    atom_list : [Atom]
        Flat list of coordinates.
    atom_group_s : int, optional
        Size of atom groups.

    Returns
    -------
    polymer : Polypeptide
        `Polymer` object containing atom coords converted `Monomers`.

    Raises
    ------
    ValueError
        Raised if `atom_group_s` != 4 or 5
    """
    atom_labels = ["N", "CA", "C", "O", "CB"]
    atom_elements = ["N", "C", "C", "O", "C"]
    atoms_coords = [
        atom_list[x : x + atom_group_s] for x in range(0, len(atom_list), atom_group_s)
    ]
    atoms = [[Atom(x[0], x[1]) for x in zip(y, atom_elements)] for y in atoms_coords]
    if atom_group_s == 5:
        monomers = [Residue(OrderedDict(zip(atom_labels, x)), "ALA") for x in atoms]
    elif atom_group_s == 4:
        monomers = [Residue(OrderedDict(zip(atom_labels, x)), "GLY") for x in atoms]
    else:
        raise ValueError(
            "Parameter atom_group_s must be 4 or 5 so atoms can be labeled correctly."
        )
    polymer = Polypeptide(monomers=monomers)
    return polymer


def flat_list_to_dummy_chain(atom_list, atom_group_s=1):
    """Converts flat list of coordinates into dummy C-alpha carbons

    Parameters
    ----------
    atom_list : [Atom]
        Flat list of co-ordinates.
    atom_group_s : int, optional
        Size of atom groups.

    Returns
    -------
    polymer : Polypeptide
        `Polymer` object containing atom coord converted `Monomers`
        with 'DUM' atom name.

    """
    atom_labels = ["CA"]
    atom_elements = ["C"]
    atoms_coords = [
        atom_list[x : x + atom_group_s] for x in range(0, len(atom_list), atom_group_s)
    ]
    atoms = [[Atom(x[0], x[1]) for x in zip(y, atom_elements)] for y in atoms_coords]
    monomers = [Residue(OrderedDict(zip(atom_labels, x)), "DUM") for x in atoms]
    polymer = Polypeptide(monomers=monomers)
    return polymer


def align(target, mobile, target_i=0, mobile_i=0):
    """Aligns one Polypeptide (mobile) to another (target).

    Notes
    -----
    This function directly modifies atoms of the mobile Polypeptide!
    It does not return a new object.

    Parameters
    ----------
    target : Polypeptide
        Polypeptide to be aligned to.
    mobile : Polypeptide
        Polypeptide to be moved during alignment.
    target_i : int, optional
        Index of `Residue` in target to align to.
    mobile_i : int, optional
        Index of `Residue` in mobile to be aligned.
    """
    # First, align N->CA vectors.
    s1, e1, s2, e2 = [
        x._vector
        for x in [
            mobile[mobile_i]["N"],
            mobile[mobile_i]["CA"],
            target[target_i]["N"],
            target[target_i]["CA"],
        ]
    ]
    translation, angle, axis, point = find_transformations(
        s1, e1, s2, e2, radians=False
    )
    # Rotation first, Then translation.
    mobile.rotate(angle=angle, axis=axis, point=point, radians=False)
    mobile.translate(vector=translation)
    # Second, rotate about N->CA axis to align CA->C vectors.
    angle = dihedral(
        mobile[mobile_i]["C"],
        mobile[mobile_i]["N"],
        mobile[mobile_i]["CA"],
        target[target_i]["C"],
    )
    axis = target[target_i]["CA"] - target[target_i]["N"]
    point = target[target_i]["N"]._vector
    mobile.rotate(angle=angle, axis=axis, point=point)
    return


class Polypeptide(Polymer):
    """Container for `Residues`, inherits from `Polymer`.

    Parameters
    ----------
    monomers : Residue or [Residue], optional
        `Residue` or list containing `Residue` objects to form the
        `Polypeptide`.
    polymer_id : str, optional
        An ID that the user can use to identify the `Polypeptide`. This is
        used when generating a pdb file using `Polypeptide().pdb`.
    parent : ampal.Assembly, optional
        Reference to `Assembly` containing the `Polymer`.
    sl : int, optional
        The default smoothing level used when calculating the
        backbone primitive.

    Attributes
    ----------
    id : str
        `Polypeptide` ID
    parent : ampal.Assembly or None
        Reference to `Assembly` containing the `Polypeptide`
    molecule_type : str
        A description of the type of `Polymer` i.e. Protein, DNA etc.
    ligands : ampal.LigandGroup
        A `LigandGroup` containing all the `Ligands` associated with this
        `Polypeptide` chain.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.
    sl : int
        The default smoothing level used when calculating the
        backbone primitive.

    Raises
    ------
    TypeError
        `Polymer` type objects can only be initialised empty or using
        a `Monomer`.
    """

    def __init__(self, monomers=None, polymer_id=" ", parent=None, sl=2):
        super().__init__(
            monomers=monomers,
            polymer_id=polymer_id,
            molecule_type="protein",
            parent=parent,
            sl=sl,
        )

    def __add__(self, other):
        if isinstance(other, Polymer):
            merged_polymer = self._monomers + other._monomers
        else:
            raise TypeError("Only Polymer objects may be merged with a Polymer.")
        return Polypeptide(monomers=merged_polymer, polymer_id=self.id)

    def __getitem__(self, item):
        if isinstance(item, str):
            id_dict = {str(m.id): m for m in self._monomers}
            return id_dict[item]
        elif isinstance(item, int):
            return self._monomers[item]
        return Polypeptide(self._monomers[item], polymer_id=self.id)

    def __repr__(self):
        if len(self.sequence) > 15:
            seq = self.sequence[:12] + "..."
        else:
            seq = self.sequence
        return "<Polypeptide containing {} {}. Sequence: {}>".format(
            len(self._monomers),
            "Residue" if len(self._monomers) == 1 else "Residues",
            seq,
        )

    def get_slice_from_res_id(self, start, end):
        """Returns a new `Polypeptide` containing the `Residues` in start/end range.

        Parameters
        ----------
        start : str
            string representing start residue id (PDB numbering)
        end : str
            string representing end residue id (PDB numbering)

        Returns
        -------
        slice_polymer : Polymer
            Polymer containing the residue range specified by start-end
        """

        id_dict = {str(m.id): m for m in self._monomers}
        slice_polymer = Polypeptide(
            [id_dict[str(x)] for x in range(int(start), int(end) + 1)], self.id
        )
        return slice_polymer

    @property
    def backbone(self):
        """Returns a new `Polymer` containing only the backbone atoms.

        Notes
        -----
        Metadata is not currently preserved from the parent object.
        Sequence data is retained, but only the main chain atoms are retained.

        Returns
        -------
        bb_poly : Polypeptide
            Polymer containing only the backbone atoms of the original
            Polymer.
        """
        bb_poly = Polypeptide([x.backbone for x in self._monomers], self.id)
        return bb_poly

    @property
    def primitive(self):
        """Primitive of the backbone.

        Notes
        -----
        This is the average of the positions of all the CAs in frames
        of `sl` `Residues`.
        """
        cas = self.get_reference_coords()
        primitive_coords = make_primitive_extrapolate_ends(cas, smoothing_level=self.sl)
        primitive = Primitive.from_coordinates(primitive_coords)
        primitive.relabel_monomers([x.id for x in self])
        primitive.id = self.id
        primitive.parent = self
        return primitive

    @property
    def fasta(self):
        """Generates sequence data for the protein in FASTA format."""
        max_line_length = 79
        fasta_str = ">{0}:{1}|PDBID|CHAIN|SEQUENCE\n".format(
            self.parent.id.upper(), self.id
        )
        seq = self.sequence
        split_seq = [
            seq[i : i + max_line_length] for i in range(0, len(seq), max_line_length)
        ]
        for seq_part in split_seq:
            fasta_str += "{0}\n".format(seq_part)
        return fasta_str

    @property
    def sequence(self):
        """Returns the sequence of the `Polymer` as a string.

        Returns
        -------
        sequence : str
            String of the `Residue` sequence of the `Polypeptide`.
        """
        seq = [x.mol_letter for x in self._monomers]
        return "".join(seq)

    @property
    def molecular_weight(self):
        """Returns the molecular weight of the `Assembly` in Daltons."""
        return sequence_molecular_weight(self.sequence)

    @property
    def molar_extinction_280(self):
        """Returns the extinction co-efficient of the `Assembly` at 280 nm."""
        return sequence_molar_extinction_280(self.sequence)

    @property
    def isoelectric_point(self):
        """Returns the isoelectric point of the `Assembly`."""
        return sequence_isoelectric_point(self.sequence)

    @property
    def backbone_bond_lengths(self):
        """Dictionary containing backbone bond lengths as lists of floats.

        Returns
        -------
        bond_lengths : dict
            Keys are `n_ca`, `ca_c`, `c_o` and `c_n`, referring to the
            N-CA, CA-C, C=O and C-N bonds respectively. Values are
            lists of floats : the bond lengths in Angstroms.
            The lists of n_ca, ca_c and c_o are of length k for
            a Polypeptide containing k Residues. The list of c_n bonds
            is of length k-1 for a Polypeptide containing k Residues
            (C-N formed between successive `Residue` pairs).
        """
        bond_lengths = dict(
            n_ca=[distance(r["N"], r["CA"]) for r in self.get_monomers(ligands=False)],
            ca_c=[distance(r["CA"], r["C"]) for r in self.get_monomers(ligands=False)],
            c_o=[distance(r["C"], r["O"]) for r in self.get_monomers(ligands=False)],
            c_n=[
                distance(r1["C"], r2["N"])
                for r1, r2 in [(self[i], self[i + 1]) for i in range(len(self) - 1)]
            ],
        )
        return bond_lengths

    @property
    def backbone_bond_angles(self):
        """Dictionary containing backbone bond angles as lists of floats.

        Returns
        -------
        bond_angles : dict
            Keys are `n_ca_c`, `ca_c_o`, `ca_c_n` and `c_n_ca`, referring
            to the N-CA-C, CA-C=O, CA-C-N and C-N-CA angles respectively.
            Values are lists of floats : the bond angles in degrees.
            The lists of n_ca_c, ca_c_o are of length k for a `Polypeptide`
            containing k `Residues`. The list of ca_c_n and c_n_ca are of
            length k-1 for a `Polypeptide` containing k `Residues` (These
            angles are across the peptide bond, and are therefore formed
            between successive `Residue` pairs).
        """
        bond_angles = dict(
            n_ca_c=[
                angle_between_vectors(r["N"] - r["CA"], r["C"] - r["CA"])
                for r in self.get_monomers(ligands=False)
            ],
            ca_c_o=[
                angle_between_vectors(r["CA"] - r["C"], r["O"] - r["C"])
                for r in self.get_monomers(ligands=False)
            ],
            ca_c_n=[
                angle_between_vectors(r1["CA"] - r1["C"], r2["N"] - r1["C"])
                for r1, r2 in [(self[i], self[i + 1]) for i in range(len(self) - 1)]
            ],
            c_n_ca=[
                angle_between_vectors(r1["C"] - r2["N"], r2["CA"] - r2["N"])
                for r1, r2 in [(self[i], self[i + 1]) for i in range(len(self) - 1)]
            ],
        )
        return bond_angles

    def c_join(
        self,
        other,
        psi=-40.76,
        omega=-178.25,
        phi=-65.07,
        o_c_n_angle=None,
        c_n_ca_angle=None,
        c_n_length=None,
        relabel=True,
    ):
        """Joins other to self at the C-terminus via a peptide bond.

        Notes
        -----
        This function directly modifies self. It does not return a new object.

        Parameters
        ----------
        other: Residue or Polypeptide
        psi: float, optional
            Psi torsion angle (degrees) between final `Residue` of self
            and first `Residue` of other.
        omega: float, optional
            Omega torsion angle (degrees) between final `Residue` of
            self and first `Residue` of other.
        phi: float, optional
            Phi torsion angle (degrees) between final `Residue` of self
            and first `Residue` of other.
        o_c_n_angle: float or None, optional
            Desired angle between O, C (final `Residue` of self) and N
            (first `Residue` of other) atoms. If `None`, default value is
            taken from `ideal_backbone_bond_angles`.
        c_n_ca_angle: float or None, optional
            Desired angle between C (final `Residue` of self) and N, CA
            (first `Residue` of other) atoms. If `None`, default value is
            taken from `ideal_backbone_bond_angles`.
        c_n_length: float or None, optional
            Desired peptide bond length between final `Residue` of self
            and first `Residue` of other. If `None`, default value is taken
            from `ideal_backbone_bond_lengths`.
        relabel: bool, optional
            If `True`, `relabel_all` is run on self before returning.

        Raises
        ------
        TypeError:
            If other is not a `Residue` or a Polypeptide.
        """
        if isinstance(other, Residue):
            other = Polypeptide([other])
        if not isinstance(other, Polypeptide):
            raise TypeError(
                "Only Polypeptide or Residue objects can be joined to a Polypeptide"
            )
        if abs(omega) >= 90:
            peptide_conformation = "trans"
        else:
            peptide_conformation = "cis"
        if o_c_n_angle is None:
            o_c_n_angle = ideal_backbone_bond_angles[peptide_conformation]["o_c_n"]
        if c_n_ca_angle is None:
            c_n_ca_angle = ideal_backbone_bond_angles[peptide_conformation]["c_n_ca"]
        if c_n_length is None:
            c_n_length = ideal_backbone_bond_lengths["c_n"]
        r1 = self[-1]
        r1_ca = r1["CA"]._vector
        r1_c = r1["C"]._vector
        r1_o = r1["O"]._vector
        # p1 is point that will be used to position the N atom of r2.
        p1 = r1_o[:]
        # rotate p1 by o_c_n_angle, about axis perpendicular to the
        # r1_ca, r1_c, r1_o plane, passing through r1_c.
        axis = numpy.cross((r1_ca - r1_c), (r1_o - r1_c))
        q = Quaternion.angle_and_axis(angle=o_c_n_angle, axis=axis)
        p1 = q.rotate_vector(v=p1, point=r1_c)
        # Ensure p1 is separated from r1_c by the correct distance.
        p1 = r1_c + (c_n_length * unit_vector(p1 - r1_c))
        # rotate p1 and r1['O'] by to obtain desired psi value at the join.
        measured_psi = dihedral(r1["N"], r1["CA"], r1["C"], p1)
        q = Quaternion.angle_and_axis(angle=(psi - measured_psi), axis=(r1_c - r1_ca))
        p1 = q.rotate_vector(v=p1, point=r1_c)
        r1["O"]._vector = q.rotate_vector(v=r1_o, point=r1_c)
        # translate other so that its first N atom is at p1
        other.translate(vector=(p1 - other[0]["N"]._vector))
        # rotate other so that c_n_ca angle is correct.
        v1 = r1_c - other[0]["N"]._vector
        v2 = other[0]["CA"]._vector - other[0]["N"]._vector
        measured_c_n_ca = angle_between_vectors(v1, v2)
        axis = numpy.cross(v1, v2)
        other.rotate(
            angle=(c_n_ca_angle - measured_c_n_ca),
            axis=axis,
            point=other[0]["N"]._vector,
        )
        # rotate other to obtain desired omega and phi values at the join
        measured_omega = dihedral(r1["CA"], r1["C"], other[0]["N"], other[0]["CA"])
        other.rotate(
            angle=(omega - measured_omega),
            axis=(other[0]["N"] - r1["C"]),
            point=other[0]["N"]._vector,
        )
        measured_phi = dihedral(r1["C"], other[0]["N"], other[0]["CA"], other[0]["C"])
        other.rotate(
            angle=(phi - measured_phi),
            axis=(other[0]["CA"] - other[0]["N"]),
            point=other[0]["CA"]._vector,
        )
        self.extend(other)
        if relabel:
            self.relabel_all()
        self.tags["assigned_ff"] = False
        return

    def n_join(
        self,
        other,
        psi=-40.76,
        omega=-178.25,
        phi=-65.07,
        o_c_n_angle=None,
        c_n_ca_angle=None,
        c_n_length=None,
        relabel=True,
    ):
        """Joins other to self at the N-terminus via a peptide bond.

        Notes
        -----
        This function directly modifies self. It does not return a new object.

        Parameters
        ----------
        other: Residue or Polypeptide
        psi: float
            Psi torsion angle (degrees) between final `Residue` of other
            and first `Residue` of self.
        omega: float
            Omega torsion angle (degrees) between final `Residue` of
            other and first `Residue` of self.
        phi: float
            Phi torsion angle (degrees) between final `Residue` of other
            and first `Residue` of self.
        o_c_n_angle: float or None
            Desired angle between O, C (final `Residue` of other) and N
            (first `Residue` of self) atoms. If `None`, default value is
            taken from `ideal_backbone_bond_angles`.
        c_n_ca_angle: float or None
            Desired angle between C (final `Residue` of other) and N, CA
            (first `Residue` of self) atoms. If `None`, default value is taken
            from `ideal_backbone_bond_angles`.
        c_n_length: float or None
            Desired peptide bond length between final `Residue` of other
            and first `Residue` of self. If None, default value is taken
            from ideal_backbone_bond_lengths.
        relabel: bool
            If True, relabel_all is run on self before returning.

        Raises
        ------
        TypeError:
            If other is not a `Residue` or a `Polypeptide`
        """
        if isinstance(other, Residue):
            other = Polypeptide([other])
        if not isinstance(other, Polypeptide):
            raise TypeError(
                "Only Polypeptide or Residue objects can be joined to a Polypeptide"
            )
        if abs(omega) >= 90:
            peptide_conformation = "trans"
        else:
            peptide_conformation = "cis"
        if o_c_n_angle is None:
            o_c_n_angle = ideal_backbone_bond_angles[peptide_conformation]["o_c_n"]
        if c_n_ca_angle is None:
            c_n_ca_angle = ideal_backbone_bond_angles[peptide_conformation]["c_n_ca"]
        if c_n_length is None:
            c_n_length = ideal_backbone_bond_lengths["c_n"]
        r1 = self[0]
        r1_n = r1["N"]._vector
        r1_ca = r1["CA"]._vector
        r1_c = r1["C"]._vector
        # p1 is point that will be used to position the C atom of r2.
        p1 = r1_ca[:]
        # rotate p1 by c_n_ca_angle, about axis perpendicular to the
        # r1_n, r1_ca, r1_c plane, passing through r1_ca.
        axis = numpy.cross((r1_ca - r1_n), (r1_c - r1_n))
        q = Quaternion.angle_and_axis(angle=c_n_ca_angle, axis=axis)
        p1 = q.rotate_vector(v=p1, point=r1_n)
        # Ensure p1 is separated from r1_n by the correct distance.
        p1 = r1_n + (c_n_length * unit_vector(p1 - r1_n))
        # translate other so that its final C atom is at p1
        other.translate(vector=(p1 - other[-1]["C"]._vector))
        # Force CA-C=O-N to be in a plane, and fix O=C-N angle accordingly
        measured_dihedral = dihedral(
            other[-1]["CA"], other[-1]["C"], other[-1]["O"], r1["N"]
        )
        desired_dihedral = 180.0
        axis = other[-1]["O"] - other[-1]["C"]
        other.rotate(
            angle=(measured_dihedral - desired_dihedral),
            axis=axis,
            point=other[-1]["C"]._vector,
        )
        axis = numpy.cross(other[-1]["O"] - other[-1]["C"], r1["N"] - other[-1]["C"])
        measured_o_c_n = angle_between_vectors(
            other[-1]["O"] - other[-1]["C"], r1["N"] - other[-1]["C"]
        )
        other.rotate(
            angle=(measured_o_c_n - o_c_n_angle),
            axis=axis,
            point=other[-1]["C"]._vector,
        )
        # rotate other to obtain desired phi, omega, psi values at the join.
        measured_phi = dihedral(other[-1]["C"], r1["N"], r1["CA"], r1["C"])
        other.rotate(angle=(phi - measured_phi), axis=(r1_n - r1_ca), point=r1_ca)
        measured_omega = dihedral(other[-1]["CA"], other[-1]["C"], r1["N"], r1["CA"])
        other.rotate(
            angle=(measured_omega - omega), axis=(r1["N"] - other[-1]["C"]), point=r1_n
        )
        measured_psi = dihedral(
            other[-1]["N"], other[-1]["CA"], other[-1]["C"], r1["N"]
        )
        other.rotate(
            angle=-(measured_psi - psi),
            axis=(other[-1]["CA"] - other[-1]["C"]),
            point=other[-1]["CA"]._vector,
        )
        self._monomers = other._monomers + self._monomers
        if relabel:
            self.relabel_all()
        self.tags["assigned_ff"] = False
        return

    def tag_sidechain_dihedrals(self, force=False, tag_rotamers=False):
        """Tags each monomer with side-chain dihedral angles

        force: bool, optional
            If `True` the tag will be run even if `Residues` are
            already tagged.
        rotamers: bool, optional
            If `True` it will tag rotamers as tags.
        """
        tagged = ["chi_angles" in x.tags.keys() for x in self._monomers]
        if tag_rotamers:
            tagged = ["rotamers" in x.tags.keys() for x in self._monomers]
        if (not all(tagged)) or force:
            for monomer in self._monomers:
                if tag_rotamers:
                    chi_angles, rotamer = measure_sidechain_torsion_angles(
                        monomer, verbose=False, return_rotamers=True
                    )
                    monomer.tags["rotamers"] = rotamer
                else:
                    chi_angles = measure_sidechain_torsion_angles(
                        monomer, verbose=False
                    )
                monomer.tags["chi_angles"] = chi_angles
        return

    def tag_torsion_angles(self, force=False):
        """Tags each Monomer of the Polymer with its omega, phi and psi torsion angle.

        Parameters
        ----------
        force : bool, optional
            If `True` the tag will be run even if `Residues` are
            already tagged.
        """
        tagged = ["omega" in x.tags.keys() for x in self._monomers]
        if (not all(tagged)) or force:
            tas = measure_torsion_angles(self._monomers)
            for monomer, (omega, phi, psi) in zip(self._monomers, tas):
                monomer.tags["omega"] = omega
                monomer.tags["phi"] = phi
                monomer.tags["psi"] = psi
                monomer.tags["tas"] = (omega, phi, psi)
        return

    def rise_per_residue(self):
        """List of rise per residue values along the `Polypeptide`.
        Notes
        -----
        Calculated from `Polypeptide.primitive`."""
        return self.primitive.rise_per_residue()

    def radii_of_curvature(self):
        """ List of radius of curvature values along the `Polypeptide`."""
        return self.primitive.radii_of_curvature()

    def tag_ca_geometry(
        self, force=False, reference_axis=None, reference_axis_name="ref_axis"
    ):
        """Tags each `Residue` with rise_per_residue, radius_of_curvature and residues_per_turn.

        Parameters
        ----------
        force : bool, optional
            If `True` the tag will be run even if `Residues` are already
            tagged.
        reference_axis : list(numpy.array or tuple or list), optional
            Coordinates to feed to geometry functions that depend on
            having a reference axis.
        reference_axis_name : str, optional
            Used to name the keys in tags at `Polypeptide` and `Residue` level.
        """
        tagged = ["rise_per_residue" in x.tags.keys() for x in self._monomers]
        if (not all(tagged)) or force:
            # Assign tags None if Polymer is too short to have a primitive.
            if len(self) < 7:
                rprs = [None] * len(self)
                rocs = [None] * len(self)
                rpts = [None] * len(self)
            else:
                rprs = self.rise_per_residue()
                rocs = self.radii_of_curvature()
                rpts = residues_per_turn(self)
            for monomer, rpr, roc, rpt in zip(self._monomers, rprs, rocs, rpts):
                monomer.tags["rise_per_residue"] = rpr
                monomer.tags["radius_of_curvature"] = roc
                monomer.tags["residues_per_turn"] = rpt
        # Functions that require a reference_axis.
        if (reference_axis is not None) and (len(reference_axis) == len(self)):
            # Set up arguments to pass to functions.
            ref_axis_args = dict(
                p=self,
                reference_axis=reference_axis,
                tag=True,
                reference_axis_name=reference_axis_name,
            )
            # Run the functions.
            polymer_to_reference_axis_distances(**ref_axis_args)
            crick_angles(**ref_axis_args)
            alpha_angles(**ref_axis_args)
        return

    def valid_backbone_bond_lengths(self, atol=0.1):
        """True if all backbone bonds are within atol Angstroms of the expected distance.

        Notes
        -----
        Ideal bond lengths taken from [1].

        References
        ----------
        .. [1] Schulz, G. E, and R. Heiner Schirmer. Principles Of
           Protein Structure. New York: Springer-Verlag, 1979.

        Parameters
        ----------
        atol : float, optional
            Tolerance value in Angstoms for the absolute deviation
            away from ideal backbone bond lengths.
        """
        bond_lengths = self.backbone_bond_lengths
        a1 = numpy.allclose(
            bond_lengths["n_ca"],
            [ideal_backbone_bond_lengths["n_ca"]] * len(self),
            atol=atol,
        )
        a2 = numpy.allclose(
            bond_lengths["ca_c"],
            [ideal_backbone_bond_lengths["ca_c"]] * len(self),
            atol=atol,
        )
        a3 = numpy.allclose(
            bond_lengths["c_o"],
            [ideal_backbone_bond_lengths["c_o"]] * len(self),
            atol=atol,
        )
        a4 = numpy.allclose(
            bond_lengths["c_n"],
            [ideal_backbone_bond_lengths["c_n"]] * (len(self) - 1),
            atol=atol,
        )
        return all([a1, a2, a3, a4])

    def valid_backbone_bond_angles(self, atol=20):
        """True if all backbone bond angles are within atol degrees of their expected values.

        Notes
        -----
        Ideal bond angles taken from [1].

        References
        ----------
        .. [1] Schulz, G. E, and R. Heiner Schirmer. Principles Of
           Protein Structure. New York: Springer-Verlag, 1979.

        Parameters
        ----------
        atol : float, optional
            Tolerance value in degrees for the absolute deviation
            away from ideal backbone bond angles.
        """
        bond_angles = self.backbone_bond_angles
        omegas = [x[0] for x in measure_torsion_angles(self)]
        trans = [
            "trans" if (omega is None) or (abs(omega) >= 90) else "cis"
            for omega in omegas
        ]
        ideal_n_ca_c = [ideal_backbone_bond_angles[x]["n_ca_c"] for x in trans]
        ideal_ca_c_o = [
            ideal_backbone_bond_angles[trans[i + 1]]["ca_c_o"]
            for i in range(len(trans) - 1)
        ]
        ideal_ca_c_o.append(ideal_backbone_bond_angles["trans"]["ca_c_o"])
        ideal_ca_c_n = [ideal_backbone_bond_angles[x]["ca_c_n"] for x in trans[1:]]
        ideal_c_n_ca = [ideal_backbone_bond_angles[x]["c_n_ca"] for x in trans[1:]]
        a1 = numpy.allclose(bond_angles["n_ca_c"], [ideal_n_ca_c], atol=atol)
        a2 = numpy.allclose(bond_angles["ca_c_o"], [ideal_ca_c_o], atol=atol)
        a3 = numpy.allclose(bond_angles["ca_c_n"], [ideal_ca_c_n], atol=atol)
        a4 = numpy.allclose(bond_angles["c_n_ca"], [ideal_c_n_ca], atol=atol)
        return all([a1, a2, a3, a4])


class Residue(Monomer):
    """Represents a amino acid `Residue`.

    Parameters
    ----------
    atoms : OrderedDict, optional
        OrderedDict containing Atoms for the Monomer. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    mol_code : str, optional
        One or three letter code that represents the monomer.
    monomer_id : str, optional
        String used to identify the residue.
    insertion_code : str, optional
        Insertion code of monomer, used if reading from pdb.
    is_hetero : bool, optional
        True if is a hetero atom in pdb. Helps with PDB formatting.
    parent : ampal.Polypeptide, optional
        Reference to `Polypeptide` containing the `Residue`.

    Attributes
    ----------
    mol_code : str
        PDB molecule code that represents the `Residue`.
    insertion_code : str
        Insertion code of `Residue`, used if reading from pdb.
    is_hetero : bool
        True if is a hetero atom in pdb. Helps with PDB formatting.
    states : dict
        Contains an `OrderedDicts` containing atom information for each
        state available for the `Residue`.
    id : str
        String used to identify the residue.
    reference_atom : str
        The key that corresponds to the reference atom. This is used
        by various functions, for example backbone primitives are
        calculated using the atom defined using this key.
    parent : Polypeptide or None
        A reference to the `Polypeptide` containing this `Residue`.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.

    Raises
    ------
    ValueError
        Raised if `mol_code` is not length 1 or 3.
    """

    def __init__(
        self,
        atoms=None,
        mol_code="UNK",
        monomer_id=" ",
        insertion_code=" ",
        is_hetero=False,
        parent=None,
    ):
        super(Residue, self).__init__(atoms, monomer_id, parent=parent)
        if len(mol_code) == 3:
            self.mol_code = mol_code
            self.mol_letter = get_aa_letter(mol_code)
        elif len(mol_code) == 1:
            self.mol_code = get_aa_code(mol_code)
            self.mol_letter = mol_code
        else:
            raise ValueError(
                "Monomer requires either a 1-letter or a 3-letter "
                "amino acid code ({})".format(mol_code)
            )
        self.insertion_code = insertion_code
        self.is_hetero = is_hetero
        self.reference_atom = "CA"

    def __repr__(self):
        return "<Residue containing {} {}. Residue code: {}>".format(
            len(self.atoms), "Atom" if len(self.atoms) == 1 else "Atoms", self.mol_code
        )

    @property
    def backbone(self):
        """Returns a new `Residue` containing only the backbone atoms.

        Returns
        -------
        bb_monomer : Residue
            `Residue` containing only the backbone atoms of the original
            `Monomer`.

        Raises
        ------
        IndexError
            Raise if the `atoms` dict does not contain the backbone
            atoms (N, CA, C, O).
        """
        try:
            backbone = OrderedDict(
                [
                    ("N", self.atoms["N"]),
                    ("CA", self.atoms["CA"]),
                    ("C", self.atoms["C"]),
                    ("O", self.atoms["O"]),
                ]
            )
        except KeyError:
            missing_atoms = filter(
                lambda x: x not in self.atoms.keys(), ("N", "CA", "C", "O")
            )
            raise KeyError(
                "Error in residue {} {} {}, missing ({}) atoms. "
                "`atoms` must be an `OrderedDict` with coordinates "
                "defined for the backbone (N, CA, C, O) atoms.".format(
                    self.parent.id, self.mol_code, self.id, ", ".join(missing_atoms)
                )
            )
        bb_monomer = Residue(
            backbone,
            self.mol_code,
            monomer_id=self.id,
            insertion_code=self.insertion_code,
            is_hetero=self.is_hetero,
        )
        return bb_monomer

    @property
    def unique_id(self):
        """Generates a tuple that uniquely identifies a `Monomer` in an `Assembly`.

        Notes
        -----
        The unique_id will uniquely identify each monomer within a polymer.
        If each polymer in an assembly has a distinct id, it will uniquely
        identify each monomer within the assembly.

        The hetero-flag is defined as in Biopython as a string that is
        either a single whitespace in the case of a non-hetero atom,
        or 'H_' plus the name of the hetero-residue (e.g. 'H_GLC' in
        the case of a glucose molecule), or 'W' in the case of a water
        molecule.

        For more information, see the Biopython documentation or this
        Biopython wiki page:
        http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

        Returns
        -------
        unique_id : tuple
            unique_id[0] is the polymer_id unique_id[1] is a triple
            of the hetero-flag, the monomer id (residue number) and the
            insertion code.
        """
        if self.is_hetero:
            if self.mol_code == "HOH":
                hetero_flag = "W"
            else:
                hetero_flag = "H_{0}".format(self.mol_code)
        else:
            hetero_flag = " "
        return self.parent.id, (hetero_flag, self.id, self.insertion_code)

    @property
    def side_chain(self):
        """List of the side-chain atoms (R-group).

        Notes
        -----
        Returns empty list for glycine.

        Returns
        -------
        side_chain_atoms: list(`Atoms`)
        """
        side_chain_atoms = []
        if self.mol_code != "GLY":
            covalent_bond_graph = generate_covalent_bond_graph(
                find_covalent_bonds(self)
            )
            try:
                subgraphs = generate_bond_subgraphs_from_break(
                    covalent_bond_graph, self["CA"], self["CB"]
                )
                if len(subgraphs) == 1:
                    subgraphs = generate_bond_subgraphs_from_break(
                        subgraphs[0], self["CD"], self["N"]
                    )
                if len(subgraphs) == 2:
                    for g in subgraphs:
                        if self["CB"] in g:
                            side_chain_atoms = g.nodes()
                            break
            except:
                warning_message = "Malformed PDB for Residue {0}: {1}.".format(
                    self.id, self
                )
                if "CB" in self.atoms.keys():
                    side_chain_atoms.append(self["CB"])
                    warning_message += " Side-chain is just the CB atom."
                else:
                    warning_message += " Empty side-chain."
                warnings.warn(warning_message, MalformedPDBWarning)
        return side_chain_atoms

    @property
    def centroid(self):
        """Calculates the centroid of the residue.

        Returns
        -------
        centroid : numpy.array or None
            Returns a 3D coordinate for the residue unless a CB
            atom is not available, in which case `None` is
            returned.

        Notes
        -----
        Uses the definition of the centroid from Huang *et al* [2]_.

        References
        ----------
        .. [2] Huang ES, Subbiah S and Levitt M (1995) Recognizing Native
           Folds by the Arrangement of Hydrophobic and Polar Residues, J. Mol.
           Biol return., **252**, 709-720.
        """
        if "CB" in self.atoms:
            cb_unit_vector = unit_vector(self["CB"]._vector - self["CA"]._vector)
            return self["CA"]._vector + (3.0 * cb_unit_vector)
        return None


__author__ = (
    "Jack W. Heal, Christopher W. Wood, Gail J. Bartlett, "
    "Andrew R. Thomson, Kieran L. Hudson"
)
