"""Contains code for analysing chemical interactions in AMPAL objects."""

import itertools
import networkx

# from .ampal_databases import element_data
# from tools.components import (ch_bond_dict, component_pi_systems,
#     get_hbond_dicts, get_hbond_acceptors, get_hbond_donors, known_pi_systems,
#     get_chbond_donors, get_chbond_dict)
from .geometry import (  # type: ignore
    distance,
    angle_between_vectors,
    find_foot_on_plane,  
    centre_of_mass,
    gen_sectors,
    dihedral,
)

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
    "ZN",
]

# hbond_donors, hbond_acceptors = get_hbond_dicts(mol_codes=core_components)
# chbond_donors = get_chbond_dict(mol_codes=core_components)
# ch_bonds = ch_bond_dict(mol_codes=core_components)
# all_pi_systems = known_pi_systems()
npi_dict = {
    "ASN": ["OD1"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "GLN": ["OE1"],
}


class Interaction(object):
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

    def __init__(self, a, b, dist):
        self._a = a
        self._b = b
        self.dist = dist

    def __hash__(self):
        return hash((self._a, self._b))

    def __eq__(self, other):
        return (type(self), self._a, self._b) == (type(other), other._a, other._b)

    def __repr__(self):
        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent
        return "<Interaction between {} {}{} and {} {}{}>".format(
            self._a.res_label, am.id, ac.id, self._b.res_label, bm.id, bc.id
        )


class CovalentBond(Interaction):
    """Defines a covalent bond."""

    @property
    def a(self):
        """One `Atom` involved in the covalent bond."""
        return self._a

    @property
    def b(self):
        """One `Atom` involved in the covalent bond."""
        return self._b

    def __repr__(self):
        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent
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


class NonCovalentInteraction(Interaction):
    """A container for all non-covalent interaction.
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

    def __init__(self, donor, acceptor, dist):
        super().__init__(donor, acceptor, dist)

    @property
    def donor(self):
        """The donor `Atom` in the interaction."""
        return self._a

    @property
    def acceptor(self):
        """The acceptor in the interaction."""
        return self._b

    def __repr__(self):
        return "<Interaction between {} {}{} (donor) " "and {} {}{} (acceptor)>".format(
            self.donor.mol_code,
            self.donor.id,
            self.donor.ampal_parent.id,
            self.acceptor.mol_code,
            self.acceptor.id,
            self.acceptor.ampal_parent.id,
        )


class C5HydrogenBond(Interaction):
    def __init__(self, a, b, dist):
        super().__init__(a, b, dist)

    @property
    def donor(self):
        return self._a.ampal_parent

    @property
    def acceptor(self):
        return self._b.ampal_parent

    def __repr__(self):
        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent

        return "<C5 Hydrogen bond in residue {}-{}-{}>".format(
            ac.id, am.mol_code, am.id
        )

    @property
    def distance(self):
        return distance(self._a, self._b)

    @property
    def phi(self):

        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent

        if not "phi" in am.tags:
            ac.tag_torsion_angles()
        return am.tags["phi"]

    @property
    def psi(self):

        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent

        if not "psi" in am.tags:
            ac.tag_torsion_angles()
        return am.tags["psi"]


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

    def __init__(self, donor, acceptor, dist, ang_d, ang_a):
        super().__init__(donor, acceptor, dist)
        self.ang_d = ang_d
        self.ang_a = ang_a

    @property
    def donor_monomer(self):
        """The donor `Monomer` in the interaction."""
        return self._a.ampal_parent

    @property
    def acceptor_monomer(self):
        """The acceptor `Monomer` in the interaction."""
        return self._b.ampal_parent

    def __repr__(self):
        dm = self.donor.ampal_parent
        dc = dm.ampal_parent
        am = self.acceptor.ampal_parent
        ac = am.ampal_parent
        return "<Hydrogen Bond between ({}{}) {}-{} ||||| {}-{} ({}{})>".format(
            dm.id,
            dc.id,
            dm.mol_code,
            self.donor.res_label,
            self.acceptor.res_label,
            am.mol_code,
            am.id,
            ac.id,
        )

    @property
    def is_mainchain(self):

        if self.donor.res_label == "H" and self.acceptor.res_label == "O":
            return True
        else:
            return False

    @property
    def is_sidechain_sidechain(self):

        if self.donor.res_label != "H" and self.acceptor.res_label != "O":
            return True
        else:
            return False


class CHydrogenBond(HydrogenBond):
    def __init__(self, donor, acceptor, dist, ang_d, ang_a):
        super().__init__(donor, acceptor, dist, ang_d, ang_a)

    def __repr__(self):
        dm = self.donor.ampal_parent
        dc = dm.ampal_parent
        am = self.acceptor.ampal_parent
        ac = am.ampal_parent
        return "<C-Hydrogen Bond between ({}{}) {}-{} ||||| {}-{} ({}{})>".format(
            dm.id,
            dc.id,
            dm.mol_code,
            self.donor.res_label,
            self.acceptor.res_label,
            am.mol_code,
            am.id,
            ac.id,
        )

    @property
    def is_calpha(self):
        calpha_h_atoms = ["HA", "HA1", "HA2", "HA3"]

        if self.donor.res_label in calpha_h_atoms:
            return True

        else:
            return False


class NPiStarInteraction(object):
    """Defines an n-->pi* interaction in terms of a donor and acceptor CARBONYL BONDS."""

    def __init__(self, carbonyl_donor, carbonyl_acceptor):
        self.carbonyl_donor = carbonyl_donor
        self.carbonyl_acceptor = carbonyl_acceptor
        self.c1 = self.carbonyl_donor.a
        self.o1 = self.carbonyl_donor.b
        self.c2 = self.carbonyl_acceptor.a
        self.o2 = self.carbonyl_acceptor.b

    def __repr__(self):
        c1 = self.carbonyl_donor.a
        o1 = self.carbonyl_donor.b
        c2 = self.carbonyl_acceptor.a
        o2 = self.carbonyl_acceptor.b

        return "<n-->pi* interaction between ({}-{}) {} ||||| {} ({}-{})>".format(
            c1.ampal_parent.mol_code,
            c1.ampal_parent.id,
            o1.res_label,
            o2.res_label,
            c2.ampal_parent.mol_code,
            c2.ampal_parent.id,
        )

    @property
    def distance(self):
        """Distance between donor O atom and acceptor C atom"""
        return distance(self.o1, self.c2)

    @property
    def angle(self):
        """Angle between the O-C n->pi* bond and the acceptor C=O bond"""

        oc_vector = self.c2._vector - self.o1._vector
        co_vector = self.c2._vector - self.o2._vector
        return angle_between_vectors(oc_vector, co_vector)

    @property
    def carbonyl_dihedral(self):

        aaa = None

        if self.c1.res_label == "C":
            aaa = self.c1.ampal_parent["CA"]
        elif self.c1.res_label == "CG":
            aaa = self.c1.ampal_parent["CB"]
        elif self.c1.res_label == "CD":
            aaa = self.c1.ampal_parent["CG"]

        return dihedral(aaa, self.c1, self.o1, self.c2)

    def parameters(
        self, dist_cutoff=3.22, angle_min=95, angle_max=125, dihedral_min=120
    ):
        """Returns all N-pi* measurements, and whether these consistute a N-pi* interaction with defined parameters.
        Parameters
        ----------
        dist_cutoff : float
            Maximum distance between donor oxygen and acceptor carbonyl
        angle_min : float
            Minimum angle between O-C n-pi* interaction and acceptor C=O bond.
        angle_max : float
            Maximum angle between O-C n-pi* interaction and acceptor C=O bond.
        dihedral_min : float
            Minimum dihedral angle for CA-C-O-C (position of pi cloud).
        Returns
        -------
        answer : bool
            Whether it constitutes a CH-pi intetaction.
        measurements : dict
            Calculated measurements.
        """
        if self.distance is None or self.distance > dist_cutoff:
            return False, {
                "distance": self.distance,
                "angle": "Not calculated",
                "dihedral": "Not calculated",
            }
        if self.angle is None or self.angle < angle_min or self.angle > angle_max:
            return False, {
                "distance": self.distance,
                "angle": self.angle,
                "dihedral": "Not calculated",
            }
        if self.carbonyl_dihedral is None or abs(self.carbonyl_dihedral) > dihedral_min:
            return False, {
                "distance": self.distance,
                "angle": self.angle,
                "dihedral": self.carbonyl_dihedral,
            }
        return True, {
            "distance": self.distance,
            "angle": self.angle,
            "dihedral": self.carbonyl_dihedral,
        }


class SaltBridge(NonCovalentInteraction):
    """Defines a salt bridge in terms of a negative and positive atom"""

    def __init__(self, donor, acceptor, dist):
        super().__init__(donor, acceptor, dist)

    @property
    def pos_monomer(self):
        return self._a.ampal_parent

    @property
    def neg_monomer(self):
        return self._b.ampal_parent

    def __repr__(self):
        dm = self.donor.ampal_parent
        dc = dm.ampal_parent
        am = self.acceptor.ampal_parent
        ac = am.ampal_parent

        return "<Salt Bridge between ({}{}) {}-{} ||||| {}-{} ({}{})".format(
            dm.id,
            dc.id,
            dm.mol_code,
            self.donor.res_label,
            self.acceptor.res_label,
            am.mol_code,
            am.id,
            ac.id,
        )


class PiBase(object):
    """A container for all types of interaction, with donor and acceptor as AMPAL Monomers."""

    def __init__(self, donor, acceptor):
        self.donor = donor
        self.acceptor = acceptor

    def __repr__(self):
        return "<Interaction between {} {}{} (donor) and {} {}{} (acceptor)>".format(
            self.donor.mol_code,
            self.donor.id,
            self.donor.ampal_parent.id,
            self.acceptor.mol_code,
            self.acceptor.id,
            self.acceptor.ampal_parent.id,
        )

    @property
    def donor_monomer(self):
        return self.donor

    @property
    def acceptor_monomer(self):
        return self.acceptor


class Cation_pi(PiBase):
    def __init__(self, donor, acceptor, pi_system=None):
        super(Cation_pi, self).__init__(donor, acceptor)
        if self.donor.mol_code == "LYS":
            self.cation = self.donor["NZ"]
        elif self.donor.mol_code == "ARG":
            self.cation = self.donor["CZ"]
        else:
            raise AttributeError(
                "{0} is not a recognized cationic amino acid.".format(
                    self.donor.mol_code
                )
            )

        if pi_system:
            self.pi_system = pi_system
        elif self.acceptor_monomer.mol_code not in all_pi_systems:
            raise AttributeError(
                "{0} has no defined pi systems - it cannot act as an acceptor.".format(
                    self.acceptor_monomer.mol_code
                )
            )
        elif len(all_pi_systems[self.acceptor_monomer.mol_code]) > 1:
            raise NameError(
                "{0} has multiple pi systems - pi_system argument must be defined from {1}.".format(
                    self.acceptor_monomer.mol_code,
                    all_pi_systems[self.acceptor_monomer.mol_code].keys(),
                )
            )
        else:
            self.pi_system = list(
                all_pi_systems[self.acceptor_monomer.mol_code].keys()
            )[0]

    def __repr__(self):
        return "<Cation-pi interaction ({0}{1}) {2} ||||| {3} ({4}{5})>".format(
            self.donor.mol_code,
            self.donor.id,
            self.cation,
            self.pi_system,
            self.acceptor_monomer.mol_code,
            self.acceptor.id,
        )

    @property
    def pi_atoms(self):
        """List of AMPAL Atoms making up acceptor pi system"""
        pi_system_atoms = all_pi_systems[self.acceptor_monomer.mol_code][self.pi_system]
        return [
            self.acceptor_monomer[x]
            for x in pi_system_atoms
            if x in self.acceptor_monomer.atoms
        ]

    @property
    def pi_centre(self):
        """Coordinates as array of centre of pi system"""
        return centre_of_mass([x._vector for x in self.pi_atoms])

    def distance(self):
        """Distance between H atom and centre of pi system"""
        if self.pi_atoms:
            return distance(self.pi_centre, self.cation._vector)
        else:
            return None

    @property
    def cation_proj(self):
        """Coordinates of projection of cationic atom onto plane of pi system."""
        if len(self.pi_atoms) > 2:
            pi1 = self.pi_atoms[0]
            pi2 = self.pi_atoms[1]
            pi3 = self.pi_atoms[2]
            return find_foot_on_plane(
                pi1._vector, pi2._vector, pi3._vector, self.cation._vector
            )
        else:
            print(
                "Cationic projection cannot be defined for {0} - fewer than three atoms in the pi-system".format(
                    self
                )
            )
            return None

    @property
    def angle(self):
        """Angle between C-H bond and normal to plane of pi system"""
        if not self.cation_proj is None:

            centre_pi_cation_vector = self.pi_centre - self.cation._vector
            pi_cation_vector = self.cation_proj - self.cation._vector
            return angle_between_vectors(pi_cation_vector, centre_pi_cation_vector)
        else:
            print(
                "Angle cannot be measured for {0} - no S projection defined.".format(
                    self
                )
            )
            return None


class Met_pi(PiBase):
    def __init__(self, donor, acceptor, pi_system=None):
        super(Met_pi, self).__init__(donor, acceptor)
        self.s = donor["SD"]

        if pi_system:
            self.pi_system = pi_system
        elif self.acceptor_monomer.mol_code not in all_pi_systems:
            raise AttributeError(
                "{0} has no defined pi systems - it cannot act as an acceptor.".format(
                    self.acceptor_monomer.mol_code
                )
            )
        elif len(all_pi_systems[self.acceptor_monomer.mol_code]) > 1:
            raise NameError(
                "{0} has multiple pi systems - pi_system argument must be defined from {1}.".format(
                    self.acceptor_monomer.mol_code,
                    all_pi_systems[self.acceptor_monomer.mol_code].keys(),
                )
            )
        else:
            self.pi_system = list(
                all_pi_systems[self.acceptor_monomer.mol_code].keys()
            )[0]

    def __repr__(self):
        return "<Met-pi interaction ({0}{1}) {2} ||||| {3} ({4}{5})>".format(
            self.donor.mol_code,
            self.donor.id,
            self.s,
            self.pi_system,
            self.acceptor.mol_code,
            self.acceptor.id,
        )

    @property
    def s_atom(self):
        """Donor C atom as AMPAL Atom"""
        return self.donor["SD"]

    @property
    def pi_atoms(self):
        """List of AMPAL Atoms making up acceptor pi system"""
        pi_system_atoms = all_pi_systems[self.acceptor_monomer.mol_code][self.pi_system]
        return [
            self.acceptor_monomer[x]
            for x in pi_system_atoms
            if x in self.acceptor_monomer.atoms
        ]

    @property
    def pi_centre(self):
        """Coordinates as array of centre of pi system"""
        return centre_of_mass([x._vector for x in self.pi_atoms])

    @property
    def distance(self):
        """Distance between H atom and centre of pi system"""
        if self.pi_atoms:
            return distance(self.pi_centre, self.s_atom._vector)
        else:
            return None

    @property
    def s_proj(self):
        """Coordinates of projection of S atom onto plane of pi system."""
        if len(self.pi_atoms) > 2:
            pi1 = self.pi_atoms[0]
            pi2 = self.pi_atoms[1]
            pi3 = self.pi_atoms[2]
            return find_foot_on_plane(
                pi1._vector, pi2._vector, pi3._vector, self.s_atom._vector
            )
        else:
            print(
                "S projection cannot be defined for {0} - fewer than three atoms in the pi-system".format(
                    self
                )
            )
            return None

    @property
    def angle(self):
        """Angle between C-H bond and normal to plane of pi system"""
        if not self.s_proj is None:

            centre_pi_S_vector = self.pi_centre - self.s_atom._vector
            pi_S_vector = self.s_proj - self.s_atom._vector
            return angle_between_vectors(pi_S_vector, centre_pi_S_vector)
        else:
            print(
                "Angle cannot be measured for {0} - no S projection defined.".format(
                    self
                )
            )
            return None

    def parameters(self, dist_cutoff=6.0, angle_cutoff=75):
        """Returns all Met-pi measurements, and whether these consistute a Met-pi interaction with defined parameters.
        Parameters
        ----------
        dist_cutoff : float
            Maximum distance between proton and centre of pi system.
        angle_cutoff : float
            Maximum angle between C-H bond and normal to plane of pi system.
        proj_cutoff : float
            Maximum distance between projection of H onto plane of pi system and centre of pi system.
        Returns
        -------
        answer : bool
            Whether it constitutes a CH-pi intetaction.
        measurements : dict
            Calculated measurements.
        """
        if self.distance is None or self.distance > dist_cutoff:
            return False, {"distance": self.distance, "angle": "Not calculated"}
        if self.angle is None or self.angle > angle_cutoff:
            return False, {"distance": self.distance, "angle": self.angle}
        return True, {"distance": self.distance, "angle": self.angle}


class Pi_pi(PiBase):
    def __init__(self, donor, acceptor, pi_system1=None, pi_system2=None):
        super(Pi_pi, self).__init__(donor, acceptor)

        if pi_system1:
            self.pi_system1 = pi_system1
        elif self.donor_monomer.mol_code not in all_pi_systems:
            raise AttributeError(
                "{0} has no identified pi systems - it cannot take part in a pi-pi interaction.".format(
                    self.donor_monomer.mol_code
                )
            )
        elif len(all_pi_systems[self.donor_monomer.mol_code]) > 1:
            raise NameError(
                "{0} has multiple pi systems - pi_system argument must be defined from {1}".format(
                    self.donor_monomer.mol_code,
                    all_pi_systems[self.donor_monomer.mol_code].keys(),
                )
            )
        else:
            self.pi_system1 = list(all_pi_systems[self.donor_monomer.mol_code].keys())[
                0
            ]

        if pi_system2:
            self.pi_system2 = pi_system2
        elif self.acceptor_monomer.mol_code not in all_pi_systems:
            raise AttributeError(
                "{0} has no identified pi systems - it cannot take part in a pi-pi interaction.".format(
                    self.acceptor_monomer.mol_code
                )
            )
        elif len(all_pi_systems[self.acceptor_monomer.mol_code]) > 1:
            raise NameError(
                "{0] has multiple pi systems - pi system argument must be defined from {1}".format(
                    self.acceptor_monomer.mol_code,
                    all_pi_systems[self.acceptor_monomer.mol_code].keys(),
                )
            )

        else:
            self.pi_system2 = list(
                all_pi_systems[self.acceptor_monomer.mol_code].keys()
            )[0]

    def __repr__(self):
        return "<Pi-pi interaction ({0} {1} ||||| {2} {3}".format(
            self.donor_monomer.mol_code,
            self.donor_monomer.id,
            self.acceptor_monomer.id,
            self.acceptor_monomer.mol_code,
        )

    @property
    def pi_atoms1(self):
        pi_system_atoms = all_pi_systems[self.donor_monomer.mol_code][self.pi_system1]
        return [
            self.donor_monomer[x]
            for x in pi_system_atoms
            if x in self.donor_monomer.atoms
        ]

    @property
    def pi_atoms2(self):
        pi_system_atoms = all_pi_systems[self.acceptor_monomer.mol_code][
            self.pi_system2
        ]
        return [
            self.acceptor_monomer[x]
            for x in pi_system_atoms
            if x in self.acceptor_monomer.atoms
        ]

    @property
    def pi_centre1(self):
        """Coordinates as array of CoM of pi system"""
        return centre_of_mass([x._vector for x in self.pi_atoms1])

    @property
    def pi_centre2(self):
        return centre_of_mass([x._vector for x in self.pi_atoms2])

    @property
    def pi_1_proj(self):
        """projection of first pi system onto plane of second"""

        if len(self.pi_atoms2) > 2:
            pi1 = self.pi_atoms2[0]
            pi2 = self.pi_atoms2[1]
            pi3 = self.pi_atoms2[2]
            return find_foot_on_plane(
                pi1._vector, pi2._vector, pi3._vector, self.pi_centre1
            )
        else:
            print(
                "Projection cannot be defined for {0} - fewer than three atoms in the pi-system".format(
                    self
                )
            )
            return None

    @property
    def pi_2_proj(self):
        """projection of second pi system onto plane of first"""
        if len(self.pi_atoms1) > 2:
            pi1 = self.pi_atoms1[0]
            pi2 = self.pi_atoms1[1]
            pi3 = self.pi_atoms1[2]
            return find_foot_on_plane(
                pi1._vector, pi2._vector, pi3._vector, self.pi_centre2
            )
        else:
            print(
                "Projection cannot be defined for {0} - fewer than three atoms in the pi-system".format(
                    self
                )
            )
            return None

    @property
    def distance(self):
        if self.pi_atoms1 and self.pi_atoms2:
            return distance(self.pi_centre1, self.pi_centre2)

    @property
    def planar_angle(self):
        vec1 = self.pi_centre1 - self.pi_1_proj
        vec2 = self.pi_centre2 - self.pi_2_proj

        return angle_between_vectors(vec1, vec2)

    @property
    def orientational_angle(self):
        vec1 = self.donor_monomer["CB"] - self.pi_centre1
        vec2 = self.acceptor_monomer["CB"] - self.pi_centre2

        return angle_between_vectors(vec1, vec2)


class CH_pi(PiBase):
    """Defines a CH-pi interaction in terms of donor C and H atoms and acceptor pi-system.
    Parameters
    ----------
    donor : AMPAL Monomer
    donor_atoms : dict
        Dictionary containing donor C and H atoms, as strings corresponding to PDB atom labels.
    acceptor : AMPAL Monomer
    pi_system : str
        Code for pi system corresponding to pi_systems table in components database.
        Required if acceptor has multiple pi systems (e.g., Trp).
    """

    def __init__(self, donor, donor_atoms, acceptor, pi_system=None):
        super(CH_pi, self).__init__(donor, acceptor)
        self.c = donor_atoms["C"]
        self.h = donor_atoms["H"]
        if pi_system:
            self.pi_system = pi_system
        elif self.acceptor_monomer.mol_code not in all_pi_systems:
            raise AttributeError(
                "{0} has no defined pi systems - it cannot act as an acceptor.".format(
                    self.acceptor_monomer.mol_code
                )
            )
        elif len(all_pi_systems[self.acceptor_monomer.mol_code]) > 1:
            raise NameError(
                "{0} has multiple pi systems - pi_system argument must be defined from {1}.".format(
                    self.acceptor_monomer.mol_code,
                    all_pi_systems[self.acceptor_monomer.mol_code].keys(),
                )
            )
        else:
            self.pi_system = list(
                all_pi_systems[self.acceptor_monomer.mol_code].keys()
            )[0]

    def __repr__(self):
        return "<CH-pi interaction ({0}{1}) {2}-{3} ||||| {4} ({5}{6})>".format(
            self.donor.mol_code,
            self.donor.id,
            self.c,
            self.h,
            self.pi_system,
            self.acceptor.mol_code,
            self.acceptor.id,
        )

    @property
    def c_atom(self):
        """Donor C atom as AMPAL Atom"""
        return self.donor[self.c]

    @property
    def h_atom(self):
        """Donor H atom as AMPAL Atom"""
        return self.donor[self.h]

    @property
    def pi_atoms(self):
        """List of AMPAL Atoms making up acceptor pi system"""
        pi_system_atoms = all_pi_systems[self.acceptor_monomer.mol_code][self.pi_system]
        return [
            self.acceptor_monomer[x]
            for x in pi_system_atoms
            if x in self.acceptor_monomer.atoms
        ]

    @property
    def h_proj(self):
        """Coordinates of projection of H atom onto plane of pi system."""
        if len(self.pi_atoms) > 2:
            pi1 = self.pi_atoms[0]
            pi2 = self.pi_atoms[1]
            pi3 = self.pi_atoms[2]
            return find_foot_on_plane(
                pi1._vector, pi2._vector, pi3._vector, self.h_atom._vector
            )
        else:
            print(
                "H projection cannot be defined for {0} - fewer than three atoms in the pi-system".format(
                    self
                )
            )
            return None

    @property
    def pi_centre(self):
        """Coordinates as array of centre of pi system"""
        return centre_of_mass([x._vector for x in self.pi_atoms])

    @property
    def distance(self):
        """Distance between H atom and centre of pi system"""
        if self.pi_atoms:
            return distance(self.pi_centre, self.h_atom._vector)
        else:
            return None

    @property
    def proj_dist(self):
        """Distance between projection of H atom onto plane of pi system and centre of pi system"""
        if not self.h_proj is None:
            return distance(self.h_proj, self.pi_centre)
        else:
            print(
                "Projection distance cannot be measured for {0} - no H projection defined.".format(
                    self
                )
            )
            return None

    @property
    def angle(self):
        """Angle between C-H bond and normal to plane of pi system"""
        if not self.h_proj is None:
            pi_H_vector = self.h_proj - self.h_atom._vector
            CH_vector = self.h_atom._vector - self.c_atom._vector
            return angle_between_vectors(pi_H_vector, CH_vector)
        else:
            print(
                "Angle cannot be measured for {0} - no H projection defined.".format(
                    self
                )
            )
            return None

    def parameters(self, dist_cutoff=3.5, angle_cutoff=55, proj_cutoff=2):
        """Returns all CH-pi measurements, and whether these consistute a CH-pi interaction with defined parameters.
        Parameters
        ----------
        dist_cutoff : float
            Maximum distance between proton and centre of pi system.
        angle_cutoff : float
            Maximum angle between C-H bond and normal to plane of pi system.
        proj_cutoff : float
            Maximum distance between projection of H onto plane of pi system and centre of pi system.
        Returns
        -------
        answer : bool
            Whether it constitutes a CH-pi intetaction.
        measurements : dict
            Calculated measurements.
        """
        if self.distance is None or self.distance > dist_cutoff:
            return False, {
                "distance": self.distance,
                "angle": "Not calculated",
                "proj_distance": "Not calculated",
            }
        if self.angle is None or self.angle > angle_cutoff:
            return False, {
                "distance": self.distance,
                "angle": self.angle,
                "proj_distance": "Not calculated",
            }
        if self.proj_dist is None or self.proj_dist > proj_cutoff:
            return False, {
                "distance": self.distance,
                "angle": self.angle,
                "proj_distance": self.proj_dist,
            }
        return True, {
            "distance": self.distance,
            "angle": self.angle,
            "proj_distance": self.proj_dist,
        }


def covalent_bonds(atoms, threshold=1.1):
    """Returns all the covalent bonds in a list of `Atom` pairs.
    Notes
    -----
    Uses information `element_data`, which can be accessed directly
    through this module i.e. `isambard.ampal.interactions.element_data`.
    Parameters
    ----------
    atoms : [(`Atom`, `Atom`)]
        List of pairs of `Atoms`.
    threshold : float, optional
        Allows deviation from ideal covalent bond distance to be included.
        For example, a value of 1.1 would allow interactions up to 10% further
        from the ideal distance to be included.
    """
    bonds = []
    for a, b in atoms:
        bond_distance = (
            element_data[a.element.title()]["atomic radius"]
            + element_data[b.element.title()]["atomic radius"]
        ) / 100
        dist = distance(a._vector, b._vector)
        if dist <= bond_distance * threshold:
            bonds.append(CovalentBond(a, b, dist))
    return bonds


def find_covalent_bonds(ampal, max_range=2.2, threshold=1.1, tag=True):
    """Finds all covalent bonds in the AMPAL object.
    Parameters
    ----------
    ampal : AMPAL Object
        Any AMPAL object with a `get_atoms` method.
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
    """
    sectors = gen_sectors(ampal.get_atoms(), max_range * 1.1)
    bonds = []
    for sector in sectors.values():
        atoms = itertools.combinations(sector, 2)
        bonds.extend(covalent_bonds(atoms, threshold=threshold))
    bond_set = list(set(bonds))
    if tag:
        for bond in bond_set:
            a, b = bond.a, bond.b
            if "covalent_bonds" not in a.tags:
                a.tags["covalent_bonds"] = [b]
            else:
                a.tags["covalent_bonds"].append(b)
            if "covalent_bonds" not in b.tags:
                b.tags["covalent_bonds"] = [a]
            else:
                b.tags["covalent_bonds"].append(a)
    return bond_set


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


def generate_bond_subgraphs_from_break(bond_graph, atom1, atom2):
    """Splits the bond graph between two atoms to producing subgraphs.
    Notes
    -----
    This will not work if there are cycles in the bond graph.
    Parameters
    ----------
    bond_graph: networkx.Graph
        Graph of covalent bond network
    atom1: isambard.ampal.Atom
        First atom in the bond.
    atom2: isambard.ampal.Atom
        Second atom in the bond.
    Returns
    -------
    subgraphs: [networkx.Graph]
        A list of subgraphs generated when a bond is broken in the covalent
        bond network.
    """
    bond_graph.remove_edge(atom1, atom2)
    try:
        subgraphs = list(networkx.connected_component_subgraphs(bond_graph, copy=False))
    finally:
        # Add edge
        bond_graph.add_edge(atom1, atom2)
    return subgraphs


def hydrogen_bonds(
    atoms, dist_range=(1.5, 2.7), angular_cutoff=90.0, water_donors=True
):
    """Returns a list of hydrogen bonds found in the input structure.
    Parameters
    ----------
    atoms: [Atom]
        Atoms to be analysed.
    dist_range: (float, float)
        Minimum and maximum distance for interaction.
    angular_cutoff: float
        Minimum interaction angle.
    water_donors : bool
        True to include waters as donors if water protons are missing, e.g. for crystal structures but not models.
    Returns
    -------
    hbonds: [HydrogenBonds]
        A list of the hydrogen bonds found in the ampal object.
    """
    donors = []
    acceptors = []
    donor_waters = []
    for atom in atoms:
        component = atom.ampal_parent.mol_code
        if component not in hbond_donors:
            hbond_donors[component] = get_hbond_donors(component)
            hbond_acceptors[component] = get_hbond_acceptors(component)
        if atom.res_label in hbond_donors[component]:
            donors.append(atom)
        elif atom.res_label in hbond_acceptors[component]:
            acceptors.append(atom)
        if component == "HOH" and atom.res_label == "O":
            donor_waters.append(atom)
    potential_hbonds = []
    for d in donors:
        for a in acceptors:
            dist = distance(d._vector, a._vector)
            if dist_range[0] < dist < dist_range[1]:
                potential_hbonds.append(((d, a), dist))
    if water_donors:
        for d in donor_waters:
            for a in acceptors:
                if a.ampal_parent.mol_code == "HOH":
                    continue
                dist = distance(d._vector, a._vector)
                if dist_range[0] + 0.96 < dist < dist_range[1] + 0.96:
                    potential_hbonds.append(((d, a), dist - 0.96))
    hbonds = []
    for ((da, aa), dist) in potential_hbonds:
        da_mc = da.ampal_parent.mol_code
        aa_mc = aa.ampal_parent.mol_code
        bv = aa._vector - da._vector  # hbond vector
        #  Donor and hydrogen
        if da.ampal_parent.mol_code == "HOH" and da.res_label == "O":
            dba = None
        else:
            dv = (
                da.ampal_parent.atoms[hbond_donors[da_mc][da.res_label]]._vector
                - da._vector
            )
            dba = angle_between_vectors(dv, bv)
            if dba < angular_cutoff:
                continue
        # Acceptor and acceptor-antecedent
        baa = None
        fail = False
        for aaa in hbond_acceptors[aa_mc][aa.res_label]:
            if aaa not in aa.ampal_parent.atoms:
                continue
            av = aa.ampal_parent.atoms[aaa]._vector - aa._vector
            baa = angle_between_vectors(-bv, av)
            if baa < angular_cutoff:
                fail = True
        if not fail:
            hbond = HydrogenBond(da, aa, dist, dba, baa)
            hbonds.append(hbond)
    return hbonds


def C_hydrogen_bonds(atoms, dist_range=(1.5, 2.7), angular_cutoff=90.0):
    """Returns a list of C_hydrogen bonds found in the input structure.
    Parameters
    ----------
    atoms: [Atom]
        Atoms to be analysed.
    dist_range: (float, float)
        Minimum and maximum distance for interaction.
    angular_cutoff: float
        Minimum interaction angle.
    water_donors : False:
        True to include waters as donors if water protons are missing, e.g. for crystal structures but not models.
    Returns
    -------
    chbonds: [HydrogenBonds]
        A list of the chydrogen bonds found in the ampal object.
    """
    donors = []
    acceptors = []
    donor_waters = []
    for atom in atoms:
        component = atom.ampal_parent.mol_code
        if component not in chbond_donors:
            chbond_donors[component] = get_chbond_donors(component)
            hbond_acceptors[component] = get_hbond_acceptors(component)
        if atom.res_label in chbond_donors[component]:
            donors.append(atom)
        elif atom.res_label in hbond_acceptors[component]:
            acceptors.append(atom)
        if component == "HOH" and atom.res_label == "O":
            donor_waters.append(atom)
    potential_chbonds = []
    for d in donors:
        for a in acceptors:
            dist = distance(d._vector, a._vector)
            if dist_range[0] < dist < dist_range[1]:
                potential_chbonds.append(((d, a), dist))
    chbonds = []
    for ((da, aa), dist) in potential_chbonds:
        da_mc = da.ampal_parent.mol_code
        aa_mc = aa.ampal_parent.mol_code
        bv = aa._vector - da._vector  # hbond vector
        #  Donor and hydrogen

        dv = (
            da.ampal_parent.atoms[chbond_donors[da_mc][da.res_label]]._vector
            - da._vector
        )
        dba = angle_between_vectors(dv, bv)
        if dba < angular_cutoff:
            continue
        # Acceptor and acceptor-antecedent
        baa = None
        fail = False
        for aaa in hbond_acceptors[aa_mc][aa.res_label]:
            if aaa not in aa.ampal_parent.atoms:
                continue
            av = aa.ampal_parent.atoms[aaa]._vector - aa._vector
            baa = angle_between_vectors(-bv, av)
            if baa < angular_cutoff:
                fail = True
        if not fail:
            chbond = CHydrogenBond(da, aa, dist, dba, baa)
            chbonds.append(chbond)
    return chbonds


def find_C5HydrogenBonds(ampal, dist_range=(1.5, 2.61), phi=140, psi=140):
    """Returns a list of C5 hydrogen bonds found in the input structure involving
    the 20 proteinogenic amino acids. Parameters taken from Newberry & Raines, Nat Chem Bio 2016
    Parameters
    ----------
    ampal : BaseAmpal or subclass
        AMPAL object to be analysed.
    dist_range : (float,float)
        Minimum and maximum distance for interaction
    phi : int
        Max absolute phi angle
    psi : int
        Max absolute psi angle
    Returns
    -------
    c5_hbonds : list
        List of C5 hydrogen bonds
    """
    c5_hbonds = []
    potential_c5s = [
        C5HydrogenBond(x["N"], x["C"], distance(x["N"], x["C"]))
        for x in ampal.get_monomers(ligands=False)
    ]
    for c in potential_c5s:
        if c.phi and c.psi:
            if (
                c.distance > dist_range[0]
                and c.distance < dist_range[1]
                and abs(c.phi) > phi
                and abs(c.psi) > psi
            ):
                c5_hbonds.append(c)
    return c5_hbonds


# This should be expanded to a dictionary with values for each hbond pair
# check figure 1 in Gail's latest ProSci paper


def find_hydrogen_bonds(
    ampal, dist_range=(1.5, 2.7), angular_cutoff=90.0, water_donors=True
):
    """Returns a list of hydrogen bonds found in the input structure involving the 20 proteinogenic amino acids.
    Parameters
    ----------
    ampal: BaseAmpal or subclass
        AMPAL object to be analysed.
    dist_range: (float, float)
        Minimum and maximum distance for interaction.
    angular_cutoff: float
        Minimum interaction angle.
    Returns
    -------
    hbonds: [HydrogenBonds]
        A list of the hydrogen bonds found in the ampal object.
    """
    sectors = gen_sectors(ampal.get_atoms(), dist_range[1] * 1.1)
    hbonds = []
    for sector in sectors.values():
        hbonds.extend(
            hydrogen_bonds(
                sector, dist_range, angular_cutoff, water_donors=water_donors
            )
        )
    return list(set(hbonds))


def find_C_hydrogen_bonds(ampal, dist_range=(1.5, 2.7), angular_cutoff=90.0):
    """Returns a list of hydrogen bonds found in the input structure involving the 20 proteinogenic amino acids.
    Parameters
    ----------
    ampal: BaseAmpal or subclass
        AMPAL object to be analysed.
    dist_range: (float, float)
        Minimum and maximum distance for interaction.
    angular_cutoff: float
        Minimum interaction angle.
    Returns
    -------
    hbonds: [HydrogenBonds]
        A list of the hydrogen bonds found in the ampal object.
    """
    sectors = gen_sectors(ampal.get_atoms(), dist_range[1] * 1.1)
    chbonds = []
    for sector in sectors.values():
        chbonds.extend(C_hydrogen_bonds(sector, dist_range, angular_cutoff))
    return list(set(chbonds))


def find_Met_pi_interactions(
    polymer, acceptor_codes=None, dist_cutoff=6.0, angle_cutoff=75, inter_chain=True
):
    """Finds all Met-Aromatic interactions based on defined parameters
    Parameters
    ----------
    polymer: Ampal object
    acceptor_codes : list or None
        optional list of mol codes of residues that will be considered as acceptors
    dist_cutoff: float
        max accepted distance between S and center of pi system that constitutes an interaction
    angle_cutoff:float
        max accepted angle between S and normal to pi system to constitute an interaction.
    inter_chain : bool
        If false, includes only Met-Pi interactions where the acceptor is in the same chain as the monomer donor
    Returns
    interactions : list
        list of Met-Pi objects that fall within the bounds
    """

    interactions = []

    if acceptor_codes:
        pi_systems = {}
        for acceptor in acceptor_codes:
            if acceptor in all_pi_systems:
                pi_systems[acceptor] = all_pi_systems[acceptor]

    monomers = []

    for monomer in polymer:
        if monomer.mol_code == "MET":
            monomers.append(monomer)

    for monomer in monomers:
        for residue in monomer.environment(include_self=False, inter_chain=inter_chain):
            if residue.mol_code not in pi_systems:
                continue
            pi_codes = pi_systems[residue.mol_code]
            for system in pi_codes:
                possible_interaction = Met_pi(
                    donor=monomer, acceptor=residue, pi_system=system
                )
                within_parameters, parameters = possible_interaction.parameters(
                    dist_cutoff=dist_cutoff, angle_cutoff=angle_cutoff
                )
                if within_parameters:
                    interactions.append(possible_interaction)
    return interactions


def find_pi_pi_interactions(
    polymer, dist_cutoff=(4.4, 5.5), angle_cutoff=(30, 60, 120)
):
    """Finds pi-pi stacking in structures, categorized into face-on and edge-on
    Parameters
    ----------
    polymer : AMPAL
        AMPAL object
    dist_cutoff: (4.4,5.5)
        Two distances for defining face-on and edge-on distances
    angle_cutoff: (30,60,120)
        planar angle definitions: Face on should be < 30 degrees, edge-on should be between 60 and 120.
    Returns
    -------
    face_interactions : []
        list of Pi_pi interaction objects that are face-on
    edge_interactions : []
        list of Pi_pi interaction objects that are edge-on
    """
    allowed_pi = ["PHE", "TRP", "TYR"]
    pi_systems = {}
    for pi in allowed_pi:
        if pi in all_pi_systems:
            pi_systems[pi] = all_pi_systems[pi]

    pis = []

    for m in polymer.get_monomers():
        if m.mol_code in allowed_pi:
            pis.append(m)

    poss_ints = []
    face_interactions = []
    edge_interactions = []

    for i in range(0, len(pis) - 1):
        for j in range(i + 1, len(pis)):
            pi_codes1 = pi_systems[pis[i].mol_code]
            pi_codes2 = pi_systems[pis[j].mol_code]

            for pi1, pi2 in zip(pi_codes1, pi_codes2):
                pipi = Pi_pi(pis[i], pis[j], pi_system1=pi1, pi_system2=pi2)
                poss_ints.append(pipi)

    for poss_int in poss_ints:

        if (
            poss_int.distance <= dist_cutoff[0]
            and poss_int.planar_angle <= angle_cutoff[0]
        ):

            face_interactions.append(poss_int)

        elif (
            poss_int.distance <= dist_cutoff[1]
            and poss_int.planar_angle >= angle_cutoff[1]
            and poss_int.planar_angle <= angle_cutoff[2]
        ):
            edge_interactions.append(poss_int)

    return face_interactions, edge_interactions


def find_cation_pi_interactions(ampal, dist_cutoff=(2.8, 6.6), angle_cutoff=(0, 30)):
    """Finds cation pi interactions between arg, lys, phe, tyr and trp usin metrics taken from ligand-interaction
    diagram from Schrodinger
    Parameters
    ----------
    ampal : AMPAL object
    dist_cutoff : (2.8,6.6)
        min/max bounds for distance
    angle_cutoff: (0,30)
        min/max bounds for angle between cation atom and normal to plane of aromatic
    Returns
    -------
    interactions : []
        list of CationPi interaction objects
    """
    allowed_donors = ["LYS", "ARG"]
    allowed_acceptors = ["PHE", "TRP", "TYR"]

    pi_systems = {}
    for acceptor in allowed_acceptors:
        if acceptor in all_pi_systems:
            pi_systems[acceptor] = all_pi_systems[acceptor]

    donors = []
    acceptors = []

    cationic_atoms = {"LYS": "NZ", "ARG": "CZ"}

    for m in ampal.get_monomers():

        if m.mol_code in allowed_donors and cationic_atoms[m.mol_code] in m.atoms:

            donors.append(m)

        if m.mol_code in allowed_acceptors:
            acceptors.append(m)

    poss_ints = []
    interactions = []

    for donor in donors:
        for acceptor in acceptors:
            pi_codes = pi_systems[acceptor.mol_code]

            for pi_code in pi_codes:
                cpi = Cation_pi(donor, acceptor, pi_system=pi_code)
                poss_ints.append(cpi)

    for cpi in poss_ints:
        if cpi.distance() < dist_cutoff[1] and cpi.distance() > dist_cutoff[0]:
            if cpi.angle <= angle_cutoff[1] and cpi.angle >= angle_cutoff[0]:
                interactions.append(cpi)

    return interactions


def find_CH_pi_interactions(
    monomer,
    acceptor_codes=None,
    dist_cutoff=3.5,
    angle_cutoff=55,
    proj_cutoff=2,
    inter_chain=True,
):
    """Finds all CH-pi interactions where an AMPAL monomer is the CH donor based on defined parameters.
    Notes
    -----
    Currently only finds interactions with acceptors that are part of Chains (i.e., not ligands) due to constraint
    of environment method.
    Parameters
    ----------
    monomer : AMPAL Monomer
        The monomer for which CH-pi interactions will be found.
    acceptor_codes : list or None
        Optional list of mol_codes of Residues that will be considered as acceptors, else all considered.
    dist_cutoff : float
        Maximum accepted distance between proton and centre of pi system to constitute an interaction.
    angle_cutoff : float
        Maximum accepted angle between CH bond and normal pi system to constitute an interaction.
    proj_cutoff : float
        Maximum accepted distance between projection of CH proton and centre of pi system to constitute an interaction.
    inter_chain : bool
        If false, only includes CH-pi interactions where the acceptor is in the same chain as the monomer donor.
    Returns
    -------
    interactions : list
        All of the interactions (as CH_pi class objects) that meet the parameters for the monomer.
    """
    interactions = []
    if monomer.mol_code in ch_bonds:
        comp_ch_bonds = ch_bonds[monomer.mol_code]
    else:
        comp_ch_bonds = monomer.comp_covalent_bonds("C", "H")
        ch_bonds[monomer.mol_code] = comp_ch_bonds
    if not comp_ch_bonds:
        return interactions
    hc_bonds = {}
    for comp_c in comp_ch_bonds:
        try:
            c_atom = monomer.atoms[comp_c]
        except KeyError:
            continue
        for comp_h in comp_ch_bonds[comp_c]:
            try:
                h_atom = monomer.atoms[comp_h]
                hc_bonds[h_atom.res_label] = c_atom.res_label
            except KeyError:
                continue
    if not hc_bonds:
        return interactions
    if acceptor_codes:
        pi_systems = {}
        for acceptor in acceptor_codes:
            if acceptor in all_pi_systems:
                pi_systems[acceptor] = all_pi_systems[acceptor]
    else:
        pi_systems = all_pi_systems
    for residue in monomer.environment(include_self=False, inter_chain=inter_chain):
        if residue.mol_code not in pi_systems:
            continue
        pi_codes = pi_systems[residue.mol_code]
        for system in pi_codes:
            if len(pi_codes[system]) > 5:
                new_proj_cutoff = proj_cutoff
            else:
                new_proj_cutoff = proj_cutoff - 0.4
            for bound_h in hc_bonds:
                possible_interaction = CH_pi(
                    donor=monomer,
                    donor_atoms={"C": hc_bonds[bound_h], "H": bound_h},
                    acceptor=residue,
                    pi_system=system,
                )
                within_parameters, parameters = possible_interaction.parameters(
                    dist_cutoff=dist_cutoff,
                    angle_cutoff=angle_cutoff,
                    proj_cutoff=new_proj_cutoff,
                )
                if within_parameters:
                    interactions.append(possible_interaction)
    return interactions


def find_CH_pi_acceptor(
    monomer,
    donor_codes=None,
    dist_cutoff=3.5,
    angle_cutoff=55,
    proj_cutoff=2,
    inter_chain=True,
    include_ligands=True,
):
    """Finds all CH-pi interactions where an AMPAL monomer is the CH acceptor based on defined parameters.
    Notes
    -----
    Currently only finds interactions with donors that are part of Chains (i.e., not ligands) due to constraint
    of environment method.
    Parameters
    ----------
    monomer : AMPAL Monomer
        The monomer for which CH-pi interactions will be found.
    donor_codes : list or None
        Optional list of mol_codes of Residues that will be considered as donors, else all considered.
    dist_cutoff : float
        Maximum accepted distance between proton and centre of pi system to constitute an interaction.
    angle_cutoff : float
        Maximum accepted angle between CH bond and normal pi system to constitute an interaction.
    proj_cutoff : float
        Maximum accepted distance between projection of CH proton and centre of pi system to constitute
        an interaction.
    inter_chain : bool
        If false, only includes CH-pi interactions where the donor is in the same chain as the monomer acceptor.
    Returns
    -------
    interactions : list
        All of the interactions (as CH_pi class objects) that meet the parameters for the monomer.
    """
    interactions = []
    if monomer.mol_code not in all_pi_systems:
        return interactions
    else:
        pi_codes = all_pi_systems[monomer.mol_code]
    for residue in monomer.environment(
        inter_chain=inter_chain, include_self=False, include_ligands=include_ligands
    ):
        if donor_codes and residue.mol_code not in donor_codes:
            continue
        if residue.mol_code in ch_bonds:
            comp_ch_bonds = ch_bonds[residue.mol_code]
        else:
            comp_ch_bonds = residue.comp_covalent_bonds("C", "H")
            ch_bonds[residue.mol_code] = comp_ch_bonds
        if not comp_ch_bonds:
            continue
        hc_bonds = {}
        for comp_c in comp_ch_bonds:
            try:
                c_atom = residue.atoms[comp_c]
            except KeyError:
                continue
            for comp_h in comp_ch_bonds[comp_c]:
                try:
                    h_atom = residue.atoms[comp_h]
                    hc_bonds[h_atom.res_label] = c_atom.res_label
                except KeyError:
                    continue
        if not hc_bonds:
            continue
        for system in pi_codes:
            if len(pi_codes[system]) > 5:
                new_proj_cutoff = proj_cutoff
            else:
                new_proj_cutoff = proj_cutoff - 0.4
            for bound_h in hc_bonds:
                possible_interaction = CH_pi(
                    donor=residue,
                    donor_atoms={"C": hc_bonds[bound_h], "H": bound_h},
                    acceptor=monomer,
                    pi_system=system,
                )
                within_parameters, parameters = possible_interaction.parameters(
                    dist_cutoff=dist_cutoff,
                    angle_cutoff=angle_cutoff,
                    proj_cutoff=new_proj_cutoff,
                )
                if within_parameters:
                    interactions.append(possible_interaction)
    return interactions


def find_CH_pis_in_list(
    monomer_list,
    donor_codes=None,
    donor_categories=None,
    acceptor_codes=None,
    dist_cutoff=3.5,
    angle_cutoff=55,
    proj_cutoff=2,
    inter_chain=True,
):
    """Finds CH-pi interactions for all monomers in a chain that fit the specified parameters.
    Parameters
    ----------
    monomer_list : list
        List of AMPAL Monomers to be investigated.
    donor_codes : list or None
        Optional list of mol_codes of donor Monomers to be included.
    donor_categories : list or None
        Optional list of categories of donor Monomers to be included.
    acceptor_codes : list or None
        Optional list of mol_codes of Monomers to be included as acceptors.
    dist_cutoff : float
        Maximum accepted distance between proton and centre of pi system to constitute an interaction.
    angle_cutoff : float
        Maximum accepted angle between CH bond and normal pi system to constitute an interaction.
    proj_cutoff : float
        Maximum accepted distance between projection of CH proton and centre of pi system to constitute an interaction.
    inter_chain : bool
        If false, only includes CH-pi interactions where the acceptor is in the same chain as the monomer donor.
    Returns
    -------
    all_interactions : dict
        A dictionary of lists of found CH_pis for each Monomer with one or more interactions.
    """
    all_interactions = {}
    for monomer in monomer_list:
        if monomer.is_solvent:
            continue
        if donor_codes and monomer.mol_code not in donor_codes:
            continue
        if donor_categories and monomer.category not in donor_categories:
            continue
        interactions = find_CH_pi_interactions(
            monomer,
            acceptor_codes=acceptor_codes,
            dist_cutoff=dist_cutoff,
            angle_cutoff=angle_cutoff,
            proj_cutoff=proj_cutoff,
            inter_chain=inter_chain,
        )
        if interactions:
            all_interactions[monomer] = interactions
    return all_interactions


def find_calpha_carbonyl_carbons(ampal):
    """Find all calpha-carbonyl-carbon bonds in an ampal object and return them
    Parameters
    ----------
    ampal : AMPAL object
    Returns
    -------
    calpha_carbonyls : list
        list of CovalentBond objects"""

    cvs = find_covalent_bonds(ampal)
    calpha_carbonyls = []

    for cv in cvs:
        if cv.a.res_label == "CA" and cv.b.res_label == "C":
            calpha_carbonyls.append(cv)

    return calpha_carbonyls


def find_carbonyls(ampal, with_ca=False):
    """Find all carbonyl bonds in an ampal object and return them
    Parameters
    ----------
    ampal : AMPAL object
    Returns
    -------
    carbonyls : list
        list of CovalentBond objects
    """

    cvs = find_covalent_bonds(ampal)
    carbonyls = []

    calpha_keys = []
    if with_ca:
        calpha_carbonyls = find_calpha_carbonyl_carbons(ampal)
        calpha_keys = [
            bond.a.ampal_parent.ampal_parent.id + bond.a.ampal_parent.id
            for bond in calpha_carbonyls
        ]

    for cv in cvs:
        if cv.a.res_label == "C" and cv.b.res_label == "O":
            if with_ca:
                if (
                    cv.a.ampal_parent.ampal_parent.id + cv.a.ampal_parent.id
                    in calpha_keys
                ):
                    carbonyls.append(cv)
                else:
                    continue
            else:
                carbonyls.append(cv)
        elif (
            cv.a.res_label == "CG"
            and cv.b.res_label == "OD1"
            or cv.b.res_label == "OD2"
        ):
            carbonyls.append(cv)
        elif (
            cv.a.res_label == "CD"
            and cv.b.res_label == "OE1"
            or cv.b.res_label == "OE2"
        ):
            carbonyls.append(cv)

    return carbonyls


def find_N_pis(
    polymer, dist_cutoff=3.22, angle_max=125, angle_min=95, dihedral_min=120
):
    """Finds n-->pi* interactions for all backbone carbonyls in a chain that fit the specified parameters
    Will not currently find n->pi* within residues (e.g. Glu OE1 to Glu C=O)
    Parameters
    ----------
    polymer : Polymer
        AMPAL polymer to be investigated
    dist_cutoff : float
        O...C distance cutoff
    angle_max : int
        Maximum O...Ci+1...Oi+1 angle
    angle_min : int
        Minimum O...Ci+1...Oi+1 angle
    dihedral_min: int
        Minimum CA---C---O---Ci+1 dihedral angle
    Returns
    -------
    interactions : list
        A list of N_PiStar_Interaction objects
    """
    interactions = []

    poss_interactions = []

    carbonyls = find_carbonyls(polymer, with_ca=True)
    for i in range(0, len(carbonyls) - 1):
        for j in range(i + 1, len(carbonyls)):
            if carbonyls[i].a.ampal_parent.id != carbonyls[j].a.ampal_parent.id:
                npistar = NPiStarInteraction(carbonyls[i], carbonyls[j])
                poss_interactions.append(npistar)
                npistar2 = NPiStarInteraction(carbonyls[j], carbonyls[i])
                poss_interactions.append(npistar2)

    for int in poss_interactions:

        if (
            int.distance <= dist_cutoff
            and int.angle >= angle_min
            and int.angle <= angle_max
            and abs(int.carbonyl_dihedral) >= dihedral_min
        ):
            interactions.append(int)

    return interactions


def salt_bridges(atoms, dist_range=(2.5, 4.0)):
    """Defines salt bridges as between positively and negatively charged residues
    Parameters
    ----------
    atoms : []
        list of AMPAL atom objects to investigate
    dist_range: ()
        min,max values for salt bridge
    Returns
    -------
    sbs : []
        list of SaltBridge objects
    """
    salt_bridge_pos = {"ARG": ["NH2", "NH1", "NE"], "LYS": ["NZ"]}
    salt_bridge_neg = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}

    pos = []
    neg = []

    for atom in atoms:
        if atom.ampal_parent.mol_code in salt_bridge_neg:
            if atom.res_label in salt_bridge_neg[atom.ampal_parent.mol_code]:
                neg.append(atom)

        elif atom.ampal_parent.mol_code in salt_bridge_pos:
            if atom.res_label in salt_bridge_pos[atom.ampal_parent.mol_code]:
                pos.append(atom)

    sbs = []
    for p in pos:
        for n in neg:
            dist = distance(p._vector, n._vector)
            if dist_range[0] < dist < dist_range[1]:
                sb = SaltBridge(p, n, dist)
                sbs.append(sb)

    return sbs


def find_salt_bridges(ampal, dist_range=(2.5, 4.0)):
    """Finds salt bridges in an ampal object
    Parameters
    ----------
    ampal : AMPAL object
    dist_range : []
        distance range for salt bridge
    Returns
    -------
    sbs : []
        list of salt bridges in ampal object
    """

    sectors = gen_sectors(ampal.get_atoms(), dist_range[1] * 1.1)
    sbs = []
    for sector in sectors.values():
        sbs.extend(salt_bridges(sector, dist_range))
    return list(set(sbs))


__author__ = "Kieran L. Hudson, Christopher W. Wood, Gail J. Bartlett"
