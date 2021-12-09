"""Contains code for parsing PDB files."""

from collections import OrderedDict
import itertools
import pathlib

from .base_ampal import Atom
from .assembly import AmpalContainer, Assembly
from .protein import Polypeptide, Residue
from .nucleic_acid import Polynucleotide, Nucleotide
from .ligands import Ligand, LigandGroup
from .amino_acids import standard_amino_acids
from .data import PDB_ATOM_COLUMNS


def load_pdb(pdb, path=True, pdb_id="", ignore_end=False):
    """Converts a PDB file into an AMPAL object.

    Parameters
    ----------
    pdb : str
        Either a path to a PDB file or a string containing PDB
        format structural data.
    path : bool, optional
        If `true`, flags `pdb` as a path and not a PDB string.
    pdb_id : str, optional
        Identifier for the `Assembly`.
    ignore_end : bool, optional
        If `false`, parsing of the file will stop when an "END"
        record is encountered.

    Returns
    -------
    ampal : ampal.Assembly or ampal.AmpalContainer
        AMPAL object that contains the structural information from
        the PDB file provided. If the PDB file has a single state
        then an `Assembly` will be returned, otherwise an
        `AmpalContainer` will be returned.
    """
    pdb_p = PdbParser(pdb, path=path, pdb_id=pdb_id, ignore_end=ignore_end)
    return pdb_p.make_ampal()


class PdbParser(object):
    """Parses a PDB file and produces and AMPAL Assembly.

    Parameters
    ----------
    pdb : str
        Either a path to a PDB file or a string containing PDB
        format structural data.
    path : bool, optional
        If `true`, flags `pdb` as a path and not a PDB string.
    pdb_id : str, optional
        Identifier for the `Assembly`.
    ignore_end : bool, optional
        If `false`, parsing of the file will stop when an "END"
        record is encountered.

    Attributes
    ----------
    proc_functions : dict
        Keys are PDB labels i.e. "ATOM" or "END", values are the
        function that processes that specific line.
    id : str
        Identifier for the `Assembly`.
    pdb_lines : [str]
        Input PDB split into line.
    new_labels : bool
        Indicates if new atom or residue labels have been found while
        parsing the PDB file.
    state : int
        Current state being appended to. This is used on multi-state
        files like NMR structures.
    pdb_parse_tree : dict
        This dictionary represents the parse tree of the PDB file.
        Each line of the structure is broken down into a key, the
        entry label, and a value, the data.
    current_line : int
        The line that is currently being parsed.
    ignore_end : bool, optional
        If `false`, parsing of the file will stop when an "END"
        record is encountered.
    """

    def __init__(self, pdb, path=True, pdb_id="", ignore_end=False):
        self.proc_functions = {
            "ATOM": self.proc_atom,
            "HETATM": self.proc_atom,
            "ENDMDL": self.change_state,
            "END": self.end,
        }
        if path:
            pdb_path = pathlib.PurePath(pdb)
            with open(str(pdb_path), "r") as inf:
                pdb_str = inf.read()
            self.id = pdb_path.stem
        else:
            pdb_str = pdb
            self.id = pdb_id
        self.pdb_lines = pdb_str.splitlines()
        self.new_labels = False
        self.state = 0
        self.pdb_parse_tree = None
        self.current_line = None
        self.ignore_end = ignore_end
        self.parse_pdb_file()

    def parse_pdb_file(self):
        """Runs the PDB parser."""
        self.pdb_parse_tree = {"info": {}, "data": {self.state: {}}}
        try:
            for line in self.pdb_lines:
                self.current_line = line
                record_name = line[:6].strip()
                if record_name in self.proc_functions:
                    self.proc_functions[record_name]()
                else:
                    if record_name not in self.pdb_parse_tree["info"]:
                        self.pdb_parse_tree["info"][record_name] = []
                    self.pdb_parse_tree["info"][record_name].append(line)
        except EOFError:
            # Raised by END record
            pass
        return

    def proc_atom(self):
        """Processes an "ATOM" or "HETATM" record."""
        atom_data = self.proc_line_coordinate(self.current_line)
        (
            at_type,
            at_ser,
            at_name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
            charge,
        ) = atom_data
        # currently active state
        a_state = self.pdb_parse_tree["data"][self.state]
        res_id = (res_seq, i_code)
        if chain_id not in a_state:
            a_state[chain_id] = (set(), OrderedDict())
        if res_id not in a_state[chain_id][1]:
            a_state[chain_id][1][res_id] = (set(), OrderedDict())
        if at_type == "ATOM":
            if res_name in standard_amino_acids.values():
                poly = "P"
            else:
                poly = "N"
        else:
            poly = "H"
        a_state[chain_id][0].add((chain_id, at_type, poly))
        a_state[chain_id][1][res_id][0].add((at_type, res_seq, res_name, i_code))
        if at_ser not in a_state[chain_id][1][res_id][1]:
            a_state[chain_id][1][res_id][1][at_ser] = [atom_data]
        else:
            a_state[chain_id][1][res_id][1][at_ser].append(atom_data)
        return

    def change_state(self):
        """Increments current state and adds a new dict to the parse tree."""
        self.state += 1
        self.pdb_parse_tree["data"][self.state] = {}
        return

    def end(self):
        """Processes an "END" record."""
        if not self.ignore_end:
            raise EOFError
        else:
            return

    def proc_line_coordinate(self, line):
        """Extracts data from columns in ATOM/HETATM record."""
        at_type = line[0:6].strip()  # 0
        at_ser = int(line[6:11].strip())  # 1
        at_name = line[12:16].strip()  # 2
        alt_loc = line[16].strip()  # 3
        res_name = line[17:20].strip()  # 4
        chain_id = line[21].strip()  # 5
        res_seq = int(line[22:26].strip())  # 6
        i_code = line[26].strip()  # 7
        x = float(line[30:38].strip())  # 8
        y = float(line[38:46].strip())  # 9
        z = float(line[46:54].strip())  # 10
        occupancy = float(line[54:60].strip())  # 11
        temp_factor = float(line[60:66].strip())  # 12
        element = line[76:78].strip()  # 13
        charge = line[78:80].strip()  # 14
        if at_name not in PDB_ATOM_COLUMNS:
            PDB_ATOM_COLUMNS[at_name] = line[12:16]
            self.new_labels = True
        return (
            at_type,
            at_ser,
            at_name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
            charge,
        )

    # Generate PDB from parse tree
    def make_ampal(self):
        """Generates an AMPAL object from the parse tree.

        Notes
        -----
        Will create an `Assembly` if there is a single state in the
        parese tree or an `AmpalContainer` if there is more than one.
        """
        data = self.pdb_parse_tree["data"]
        if len(data) > 1:
            ac = AmpalContainer(id=self.id)
            for state, chains in sorted(data.items()):
                if chains:
                    ac.append(
                        self.proc_state(chains, self.id + "_state_{}".format(state + 1))
                    )
            return ac
        elif len(data) == 1:
            return self.proc_state(data[0], self.id)
        else:
            raise ValueError("Empty parse tree, check input PDB format.")

    def proc_state(self, state_data, state_id):
        """Processes a state into an `Assembly`.

        Parameters
        ----------
        state_data : dict
            Contains information about the state, including all
            the per line structural data.
        state_id : str
            ID given to `Assembly` that represents the state.
        """
        assembly = Assembly(assembly_id=state_id)
        for k, chain in sorted(state_data.items()):
            assembly._molecules.append(self.proc_chain(chain, assembly))
        return assembly

    def proc_chain(self, chain_info, parent):
        """Converts a chain into a `Polymer` type object.

        Parameters
        ----------
        chain_info : (set, OrderedDict)
            Contains a set of chain labels and atom records.
        parent : ampal.Assembly
            `Assembly` used to assign `parent` on created
            `Polymer`.

        Raises
        ------
        ValueError
            Raised if multiple or unknown atom types found
            within the same chain.
        AttributeError
            Raised if unknown `Monomer` type encountered.
        """
        hetatom_filters = {"nc_aas": self.check_for_non_canonical}

        polymer = False
        chain_labels, chain_data = chain_info
        chain_label = list(chain_labels)[0]
        monomer_types = {x[2] for x in chain_labels if x[2]}
        if ("P" in monomer_types) and ("N" in monomer_types):
            raise ValueError('Malformed PDB, multiple "ATOM" types in a single chain.')
        # Changes Polymer type based on chain composition
        if "P" in monomer_types:
            polymer_class = Polypeptide
            polymer = True
        elif "N" in monomer_types:
            polymer_class = Polynucleotide
            polymer = True
        elif "H" in monomer_types:
            polymer_class = LigandGroup
        else:
            raise AttributeError("Malformed parse tree, check inout PDB.")
        chain = polymer_class(polymer_id=chain_label[0], parent=parent)
        # Changes where the ligands should go based on the chain composition
        if polymer:
            chain.ligands = LigandGroup(polymer_id=chain_label[0], parent=parent)
            ligands = chain.ligands
        else:
            ligands = chain

        for residue in chain_data.values():
            res_info = list(residue[0])[0]
            if res_info[0] == "ATOM":
                chain._monomers.append(self.proc_monomer(residue, chain))
            elif res_info[0] == "HETATM":
                mon_cls = None
                on_chain = False
                for filt_func in hetatom_filters.values():
                    filt_res = filt_func(residue)
                    if filt_res:
                        mon_cls, on_chain = filt_res
                        break
                    mon_cls = Ligand
                if on_chain:
                    chain._monomers.append(
                        self.proc_monomer(residue, chain, mon_cls=mon_cls)
                    )
                else:
                    ligands._monomers.append(
                        self.proc_monomer(residue, chain, mon_cls=mon_cls)
                    )
            else:
                raise ValueError("Malformed PDB, unknown record type for data")
        return chain

    def proc_monomer(self, monomer_info, parent, mon_cls=False):
        """Processes a records into a `Monomer`.

        Parameters
        ----------
        monomer_info : (set, OrderedDict)
            Labels and data for a monomer.
        parent : ampal.Polymer
            `Polymer` used to assign `parent` on created
            `Monomer`.
        mon_cls : `Monomer class or subclass`, optional
            A `Monomer` class can be defined explicitly.
        """
        monomer_labels, monomer_data = monomer_info
        if len(monomer_labels) > 1:
            raise ValueError(
                "Malformed PDB, single monomer id with "
                "multiple labels. {}".format(monomer_labels)
            )
        else:
            monomer_label = list(monomer_labels)[0]
        if mon_cls:
            monomer_class = mon_cls
            het = True
        elif monomer_label[0] == "ATOM":
            if monomer_label[2] in standard_amino_acids.values():
                monomer_class = Residue
            else:
                monomer_class = Nucleotide
            het = False
        else:
            raise ValueError("Unknown Monomer type.")
        monomer = monomer_class(
            atoms=None,
            mol_code=monomer_label[2],
            monomer_id=monomer_label[1],
            insertion_code=monomer_label[3],
            is_hetero=het,
            parent=parent,
        )
        monomer.states = self.gen_states(monomer_data.values(), monomer)
        monomer._active_state = sorted(monomer.states.keys())[0]
        return monomer

    def gen_states(self, monomer_data, parent):
        """Generates the `states` dictionary for a `Monomer`.

        monomer_data : list
            A list of atom data parsed from the input PDB.
        parent : ampal.Monomer
            `Monomer` used to assign `parent` on created
            `Atoms`.
        """
        states = {}
        for atoms in monomer_data:
            for atom in atoms:
                state = "A" if not atom[3] else atom[3]
                if state not in states:
                    states[state] = OrderedDict()
                states[state][atom[2]] = Atom(
                    tuple(atom[8:11]),
                    atom[13],
                    atom_id=atom[1],
                    res_label=atom[2],
                    occupancy=atom[11],
                    bfactor=atom[12],
                    charge=atom[14],
                    state=state,
                    parent=parent,
                )

        # This code is to check if there are alternate states and populate any
        # both states with the full complement of atoms
        states_len = [(k, len(x)) for k, x in states.items()]
        if (len(states) > 1) and (len(set([x[1] for x in states_len])) > 1):
            for t_state, t_state_d in states.items():
                new_s_dict = OrderedDict()
                for k, v in states[
                    sorted(states_len, key=lambda x: x[0])[0][0]
                ].items():
                    if k not in t_state_d:
                        c_atom = Atom(
                            v._vector,
                            v.element,
                            atom_id=v.id,
                            res_label=v.res_label,
                            occupancy=v.tags["occupancy"],
                            bfactor=v.tags["bfactor"],
                            charge=v.tags["charge"],
                            state=t_state[0],
                            parent=v.parent,
                        )
                        new_s_dict[k] = c_atom
                    else:
                        new_s_dict[k] = t_state_d[k]
                states[t_state] = new_s_dict
        return states

    # HETATM filters
    @staticmethod
    def check_for_non_canonical(residue):
        """Checks to see if the residue is non-canonical."""
        res_label = list(residue[0])[0][2]
        atom_labels = {
            x[2] for x in itertools.chain(*residue[1].values())
        }  # Used to find unnatural aas
        if (all(x in atom_labels for x in ["N", "CA", "C", "O"])) and (
            len(res_label) == 3
        ):
            return Residue, True
        return None


__author__ = "Christopher W. Wood"
