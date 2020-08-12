"""Contains code for analysing chemical interactions in AMPAL objects."""

import itertools
import networkx
from networkx import Graph
import typing as t

from .base_ampal import Atom, BaseAmpal, Monomer
from .data import ELEMENT_DATA
from .geometry import distance, gen_sectors  # pylint: disable=no-name-in-module

core_components = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HOH",
]


class Interaction:
    """A container for all types of interaction with donor and acceptor.

    Parameters
    ----------
    a : ampal.Atom
        A member of a pairwise interaction.
    b : ampal.Atom
        A member of a pairwise interaction.
    dist : float
        The distance between `a` and `b`.

    Attributes
    ----------
    a : ampal.Atom
        A member of a pairwise interaction.
    b : ampal.Atom
        A member of a pairwise interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    """

    def __init__(self, a: Atom, b: Atom, dist: float):
        self._a = a
        self._b = b
        self.dist = dist

    def __hash__(self):
        return hash((self._a, self._b))

    def __eq__(self, other):
        return (type(self), self._a, self._b) == (type(other), other._a, other._b)

    def __repr__(self):
        am = self._a.parent
        ac = am.parent
        bm = self._b.parent
        bc = bm.parent
        return "<Interaction between {} {}{} and {} {}{}>".format(
            self._a.res_label, am.id, ac.id, self._b.res_label, bm.id, bc.id
        )


class CovalentBond(Interaction):
    """Defines a covalent bond."""

    @property
    def a(self) -> Atom:
        """One `Atom` involved in the covalent bond."""
        return self._a

    @property
    def b(self) -> Atom:
        """One `Atom` involved in the covalent bond."""
        return self._b

    def __repr__(self):
        am = self._a.parent
        ac = am.parent
        bm = self._b.parent
        bc = bm.parent
        return "<Covalent bond between {}{} {} {} --- {} {} {}{}>".format(
            ac.id,
            am.id,
            am.mol_code,
            self._a.res_label,
            self._b.res_label,
            bm.mol_code,
            bc.id,
            bm.id,
        )


def find_covalent_bonds(
    atoms: t.List[Atom],
    max_range: float = 2.2,
    threshold: float = 1.1,
    tag: bool = True,
):
    """Finds all covalent bonds in the AMPAL object.

    Parameters
    ----------
    atoms : [Atoms]
        A list of atoms.
    max_range : float, optional
        Used to define the sector size, so interactions at longer ranges
        will not be found.
    threshold : float, optional
        Allows deviation from ideal covalent bond distance to be included.
        For example, a value of 1.1 would allow interactions up to 10% further
        from the ideal distance to be included.
    tag : bool, optional
        If `True`, will add the covalent bond to the tags dictionary of
        each `Atom` involved in the interaction under the `covalent_bonds`
        key.

    Returns
    -------
    bonds : [CovalentBond]
        A list of `CovalentBond` objects.
    """
    atom_pairs = itertools.combinations(atoms, 2)
    bonds = []
    for a, b in atom_pairs:
        bond_distance = (
            ELEMENT_DATA[a.element.title()]["atomic radius"]
            + ELEMENT_DATA[b.element.title()]["atomic radius"]
        ) / 100
        dist = distance(a._vector, b._vector)
        if dist <= bond_distance * threshold:
            bonds.append(CovalentBond(a, b, dist))
    if tag:
        for bond in bonds:
            a, b = bond.a, bond.b
            if "covalent_bonds" not in a.tags:
                a.tags["covalent_bonds"] = [b]
            else:
                a.tags["covalent_bonds"].append(b)
            if "covalent_bonds" not in b.tags:
                b.tags["covalent_bonds"] = [a]
            else:
                b.tags["covalent_bonds"].append(a)
    return bonds


def generate_covalent_bond_graph(covalent_bonds):
    """Generates a graph of the covalent bond network described by the interactions.

    Parameters
    ----------
    covalent_bonds: [CovalentBond]
        List of `CovalentBond`.

    Returns
    -------
    bond_graph: networkx.Graph
        A graph of the covalent bond network.
    """
    bond_graph = networkx.Graph()
    for inter in covalent_bonds:
        bond_graph.add_edge(inter.a, inter.b)
    return bond_graph


def generate_bond_subgraphs_from_break(bond_graph: Graph, covalent_bond: CovalentBond):
    """Splits the bond graph between two atoms to producing subgraphs.

    Notes
    -----
    This will not work if there are cycles in the bond graph.

    Parameters
    ----------
    bond_graph: networkx.Graph
        Graph of covalent bond network
    covalent_bond: CovalentBond
        Covalent bond to break to produce subgraphs.

    Returns
    -------
    subgraphs: [networkx.Graph]
        A list of subgraphs generated when a bond is broken in the covalent
        bond network.
    """
    bond_graph.remove_edge(covalent_bond.a, covalent_bond.b)
    try:
        # pylint: disable=no-member
        subgraphs = list(networkx.connected_component_subgraphs(bond_graph, copy=False))
    finally:
        # Add edge in order to leave the base graph unchanged
        bond_graph.add_edge(covalent_bond.a, covalent_bond.b)
    return subgraphs


class NonCovalentInteraction(Interaction):
    """ A container for all non-covalent interaction.

    Parameters
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.

    Attributes
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    """

    def __init__(self, donor: Atom, acceptor: Atom, dist: float):
        super().__init__(donor, acceptor, dist)

    @property
    def donor(self) -> Atom:
        """The donor `Atom` in the interaction."""
        return self._a

    @property
    def acceptor(self) -> Atom:
        """The acceptor in the interaction."""
        return self._b

    def __repr__(self):
        return "<Interaction between {} {}{} (donor) and {} {}{} (acceptor)>".format(
            self.donor.mol_code,
            self.donor.parent.id,
            self.donor.id,
            self.acceptor.mol_code,
            self.acceptor.parent.id,
            self.acceptor.id,
        )


class HydrogenBond(NonCovalentInteraction):
    """Defines a hydrogen bond in terms of a donor and an acceptor.

    Parameters
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    ang_a : float
        Angle between the acceptor and the interaction vector.
    ang_d : float
        Angle between the donor and the interaction vector.

    Attributes
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    ang_a : float
        Angle between the acceptor and the interaction vector.
    ang_d : float
        Angle between the donor and the interaction vector.
    """

    def __init__(
        self, donor: Atom, acceptor: Atom, dist: float, ang_d: float, ang_a: float
    ):
        super().__init__(donor, acceptor, dist)
        self.ang_d = ang_d
        self.ang_a = ang_a

    @property
    def donor_monomer(self) -> Monomer:
        """The donor `Monomer` in the interaction."""
        return self._a.parent

    @property
    def acceptor_monomer(self) -> Monomer:
        """The acceptor `Monomer` in the interaction."""
        return self._b.parent

    def __repr__(self):
        dm = self.donor.parent
        dc = dm.parent
        am = self.acceptor.parent
        ac = am.parent
        return "<Hydrogen Bond between ({}{}) {}-{} ||||| {}-{} ({}{})>".format(
            dc.id,
            dm.id,
            dm.mol_code,
            self.donor.res_label,
            self.acceptor.res_label,
            am.mol_code,
            ac.id,
            am.id,
        )


class SaltBridge(NonCovalentInteraction):
    """Defines a salt bridge in terms of a negative and positive atom."""

    def __init__(self, donor: Atom, acceptor: Atom, dist: float):
        super().__init__(donor, acceptor, dist)

    @property
    def pos_monomer(self) -> Monomer:
        return self._a.parent

    @property
    def neg_monomer(self) -> Monomer:
        return self._b.parent

    def __repr__(self):
        dm = self.donor.parent
        dc = dm.parent
        am = self.acceptor.parent
        ac = am.parent

        return "<Salt Bridge between ({}{}) {}-{} ||||| {}-{} ({}{})>".format(
            dc.id,
            dm.id,
            dm.mol_code,
            self.donor.res_label,
            self.acceptor.res_label,
            am.mol_code,
            ac.id,
            am.id,
        )


PROTEIN_SALT_BRIDGE_POS = {"ARG": ["NH2", "NH1", "NE"], "LYS": ["NZ"]}
PROTEIN_SALT_BRIDGE_NEG = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}


def find_salt_bridges(
    atoms: t.List[Atom],
    positive_labels: t.Optional[t.Dict[str, t.List[str]]] = None,
    negative_labels: t.Optional[t.Dict[str, t.List[str]]] = None,
    min_dist: float = 2.5,
    max_dist: float = 4.0,
):
    """Defines salt bridges as between positively and negatively charged residues.

    Parameters
    ----------
    atoms : [Atom]
        A list of Atoms that will be searched for salt bridges.
    positive_labels : {str: [str]}, optional
        Labels that define positive atoms. If no labels are provided,
        `{"ARG": ["NH2", "NH1", "NE"], "LYS": ["NZ"]}` are used.
    negative_labels : {str: [str]}, optional
        Labels that define negative atoms. If no labels are provided,
        `{"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}` are used.
    min_dist : float
        The minimum distance the interaction will be considered a salt bridge.
    max_dist : float
        The maximum distance the interaction will be considered a salt bridge.

    Returns
    -------
    salt_bridges : [SaltBridge]
        A list of SaltBridge objects.
    """

    pos = []
    neg = []
    positive_labels = (
        PROTEIN_SALT_BRIDGE_POS if positive_labels is None else positive_labels
    )
    negative_labels = (
        PROTEIN_SALT_BRIDGE_NEG if negative_labels is None else negative_labels
    )

    for atom in atoms:
        if atom.parent.mol_code in negative_labels:
            if atom.res_label in negative_labels[atom.parent.mol_code]:
                neg.append(atom)

        elif atom.parent.mol_code in positive_labels:
            if atom.res_label in positive_labels[atom.parent.mol_code]:
                pos.append(atom)

    salt_bridges = []
    for p in pos:
        for n in neg:
            dist = distance(p._vector, n._vector)
            if min_dist <= dist <= max_dist:
                sb = SaltBridge(p, n, dist)
                salt_bridges.append(sb)

    return salt_bridges


__author__ = "Kieran L. Hudson, Christopher W. Wood, Gail J. Bartlett"

