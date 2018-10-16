import re

to_canonical = {
    'BOND': ['bond'],
    'ANGLE': ['angle'],
    'UB': ['angle', 'urey_bradley'],
    'DIHED': ['proper'],
    'IMP': ['improper'],
    'CMAP': ['cmap'],
    'HBOND': ['h-bond'],
    'VDWAALS': ['vdw_total'],
    '1-4 VDW': ['vdw_14'],
    'EEL': ['coulomb_total'],
    '1-4 EEL': ['coulomb-14'],
    'ENERGY': ['potential']
}

# Possible size keys to look for in the header
size_keys = [
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH", "MPHIA",
    "NHPARM", "NPARM", "NNB", "NRES", "NBONA", "NTHETA", "NPHIA", "NUMBND",
    "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT", "NBPER", "NGPER", "NDPER",
    "MBPER", "MGPER", "MDPER", "IFBOX", "NMXRS", "IFCAP", "NUMEXTRA"
]

data_labels = {
    "POINTERS": [31, "%FORMAT(10I8)"],
    "ATOM_NAME": ["NATOM", "%FORMAT(20a4)"],
    "CHARGE": ["NATOM", "%FORMAT(5E16.8)"],
    "ATOMIC_NUMBER": ["NATOM", "%FORMAT(10I8)"],
    "MASS": ["NATOM", "%FORMAT(5E16.8)"],
    "ATOM_TYPE_INDEX": ["NATOM", "%FORMAT(10I8)"],
    "NUMBER_EXCLUDED_ATOMS": ["NATOM", "%FORMAT(10I8)"],
    "NONBONDED_PARM_INDEX": ["NTYPES ** 2", "%FORMAT(10I8)"],
    "RESIDUE_LABEL": ["NRES", "%FORMAT(20a4)"],
    "RESIDUE_POINTER": ["NRES", "%FORMAT(10I8)"],
    "BOND_FORCE_CONSTANT": ["NUMBND", "%FORMAT(5E16.8)"],
    "BOND_EQUIL_VALUE": ["NUMBND", "%FORMAT(5E16.8)"],
    "ANGLE_FORCE_CONSTANT": ["NUMANG", "%FORMAT(5E16.8)"],
    "ANGLE_EQUIL_VALUE": ["NUMANG", "%FORMAT(5E16.8)"],
    "DIHEDRAL_FORCE_CONSTANT": ["NPTRA", "%FORMAT(5E16.8)"],
    "DIHEDRAL_PERIODICITY": ["NPTRA", "%FORMAT(5E16.8)"],
    "DIHEDRAL_PHASE": ["NPTRA", "%FORMAT(5E16.8)"],
    "SCEE_SCALE_FACTOR": ["NPTRA", "%FORMAT(5E16.8)"],
    "SCNB_SCALE_FACTOR": ["NPTRA", "%FORMAT(5E16.8)"],
    "SOLTY": ["NATYP", "%FORMAT(5E16.8)"],
    "LENNARD_JONES_ACOEF": ["(NTYPES * (NTYPES + 1)) / 2", "%FORMAT(5E16.8)"],
    "LENNARD_JONES_BCOEF": ["(NTYPES * (NTYPES + 1)) / 2", "%FORMAT(5E16.8)"],
    "BONDS_INC_HYDROGEN": ["3 * NBONH", "%FORMAT(10I8)"],
    "BONDS_WITHOUT_HYDROGEN": ["3 * NBONA", "%FORMAT(10I8)"],
    "ANGLES_INC_HYDROGEN": ["4 * NTHETH", "%FORMAT(10I8)"],
    "ANGLES_WITHOUT_HYDROGEN": ["4 * NTHETA", "%FORMAT(10I8)"],
    "DIHEDRALS_INC_HYDROGEN": ["5 * NPHIH", "%FORMAT(10I8)"],
    "DIHEDRALS_WITHOUT_HYDROGEN": ["5 * NPHIA", "%FORMAT(10I8)"],
    "EXCLUDED_ATOMS_LIST": ["NNB", "%FORMAT(10I8)"],
    "HBOND_ACOEF": ["NPHB", "%FORMAT(5E16.8)"],
    "HBOND_BCOEF": ["NPHB", "%FORMAT(5E16.8)"],
    "HBCUT": ["NPHB", "%FORMAT(5E16.8)"],
    "AMBER_ATOM_TYPE": ["NATOM", "%FORMAT(20a4)"],
    "TREE_CHAIN_CLASSIFICATION": ["NATOM", "%FORMAT(20a4)"],
    "JOIN_ARRAY": [
        "NATOM", "%FORMAT(10I8)"
    ],  # This section is filled with zeros, but is unused. We should not store it.
    "IROTAT": ["NATOM", "%FORMAT(10I8)"],
    "SOLVENT_POINTERS": ["3 if IFBOX else 0", "%FORMAT(3I8)"],
    "ATOMS_PER_MOLECULE": ["NATOM if IFBOX else 0", "%FORMAT(10I8)"],
    # "ATOMS_PER_MOLECULE": ["SOLVENT_POINTERS[1] if IFBOX else 0"], # SOLVENT_POINTERS[1] == NPSM
    "BOX_DIMENSIONS": [4, "%FORMAT(5E16.8)"],
    "CAP_INFO": ["1 if IFCAP else 0", "%FORMAT(10I8)"],
    "CAP_INFO2": ["4 if IFCAP else 0", "%FORMAT(5E16.8)"],
    "RADIUS_SET": [1, "%FORMAT(1a80)"],
    "RADII": ["NATOM", "%FORMAT(5E16.8)"],
    "SCREEN": ["NATOM", "%FORMAT(5E16.8)"],
    "IPOL": [1, "%FORMAT(1I8)"],
    # "POLARIZABILITY": [0, "%FORMAT(5E16.8)"]
    # "POLARIZABILITY": ["NATOM if IPOL else 0"]
}

atom_data_units = {
    "charge": "e / 18.2223",  # Internal units
    "mass": "g * mol ** -1",
    "length": "angstrom",
}

box_units = {
    "length": "angstrom",
    "angle": "degree",
    "center": ["a/2", "b/2", "c/2"],
}

# Box boundaries are always periodic if there is a box (i.e. IFBOX > 0).
# Else, this is a vacuum simulation where there are no boundaries. This
# is most similar to the 'shrink-wrapped' option in lammps and the EEX
# data layer
box_boundaries = {
    "x": ["'periodic' if IFBOX>0, else 'shrink-wrapped'"],
    "y": ["'periodic' if IFBOX>0, else 'shrink-wrapped'"],
    "z": ["'periodic' if IFBOX>0, else 'shrink-wrapped'"],
}

atom_property_names = {
    "ATOM_NAME": "atom_name",
    "CHARGE": "charge",
    "MASS": "mass",
    "ATOM_TYPE_INDEX": "atom_type",
    "ATOMIC_NUMBER": "atomic_number",
    # "AMBER_ATOM_TYPE": "atom_type_name",
    # "RADII" : "implicit_solvent_radius",
}

residue_store_names = ["RESIDUE_LABEL", "RESIDUE_POINTER"]
molecule_store_names = ["ATOMS_PER_MOLECULE"]

topology_store_names = [
    "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN",
    "ANGLES_WITHOUT_HYDROGEN", "DIHEDRALS_INC_HYDROGEN",
    "DIHEDRALS_WITHOUT_HYDROGEN"
]

forcefield_parameters = {
    "bond": {
        "order": 2,
        "form": "harmonic",
        "units": {
            "K": "kcal * mol ** -1 angstrom ** -2",
            "R0": "angstrom"
        },
        "column_names": {
            "BOND_FORCE_CONSTANT": "K",
            "BOND_EQUIL_VALUE": "R0"
        }
    },
    "angle": {
        "order": 3,
        "form": "harmonic",
        "units": {
            "K": "kcal * mol ** -1 radian ** -2",
            "theta0": "radian"
        },
        "column_names": {
            "ANGLE_FORCE_CONSTANT": "K",
            "ANGLE_EQUIL_VALUE": "theta0"
        }
    },
    "dihedral": {
        "order": 4,
        "form": "charmmfsw",
        "units": {
            "K": "kcal * mol ** -1",
            "n": "phase",
            "d": "radians",
        },
        "column_names": {
            "DIHEDRAL_FORCE_CONSTANT": "K",
            "DIHEDRAL_PERIODICITY": "n",
            "DIHEDRAL_PHASE": "d",
        }
    },
    "nonbond": {
        "order": None,
        "form": {
            "name": "LJ",
            "form": "AB"
        },
        "units": {
            "A": "kcal * mol ** -1 * angstrom ** 12",
            "B": "kcal * mol ** -1 * angstrom ** 6",
        },
        "column_names": {
            "LENNARD_JONES_ACOEF": "A",
            "LENNARD_JONES_BCOEF": "B",
            "NONBONDED_PARM_INDEX": ""
        },
    },
}

# Can be either keyword - equivalent
mixing_rule = ["lorentz_berthelot", "arithmetic"]

## Exclusions and scalings

# NUMBER_EXCLUDED_ATOMS has NATOMS elements and contains the number of atoms that need to be excluded from the non-bonded
# calculation loop for atom i because i is involved in a bond, angle, or torsion with those atoms.
# Points to EXCLUDED_ATOMS_LIST, which contains the atom indices for atoms which should be excluded from ith atom
# calcuation.
exclusion_sections = ["NUMBER_EXCLUDED_ATOMS", "EXCLUDED_ATOMS_LIST"]

# "canonical" amber simulations will have values of 2 (SCNB) or 1.2 (SCEE), which should allow us to use the
# dl.set_nb_scaling_factor function. If not,
scaling_sections = ["SCEE_SCALE_FACTOR", "SCNB_SCALE_FACTOR"]

#----------------------------------------------------------------------------

store_other = []
for k, v in forcefield_parameters.items():
    store_other.extend(list(v["column_names"]))
store_other.extend(residue_store_names)
store_other.extend(molecule_store_names)
store_other.extend(exclusion_sections)


def parse_format(string):
    """
    Parses an AMBER style format string.

    Example:

    >>> string = "%FORMAT(10I8)"
    >>> _parse_format(string)
    [10, int, 8]
    """
    if "FORMAT" not in string:
        raise ValueError(
            "AMBER: Did not understand format line '%s'." % string)

    pstring = string.replace("%FORMAT(", "").replace(")", "").strip()
    ret = [x for x in re.split('(\d+)', pstring) if x]

    if ret[1] == "I":
        if len(ret) != 3:
            raise ValueError(
                "AMBER: Did not understand format line '%s'." % string)
        ret[1] = int
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
    elif ret[1] == "E":
        if len(ret) != 5:
            raise ValueError(
                "AMBER: Did not understand format line '%s'." % string)
        ret[1] = float
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
        # The .8 is not interesting to us
    elif ret[1] == "a":
        if len(ret) != 3:
            raise ValueError(
                "AMBER: Did not understand format line '%s'." % string)
        ret[1] = str
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
    else:
        raise ValueError(
            "AMBER: Type symbol '%s' not understood from line '%s'." %
            (ret[1], string))

    return ret


def build_format(fmt):
    if fmt[1] == str:
        fmt = "%-" + str(fmt[2]) + "s"
    elif fmt[1] == float:
        fmt = " % " + str(fmt[2] - 1) + "." + str(fmt[4]) + "E"
    elif fmt[1] == int:
        fmt = " % " + str(fmt[2] - 1) + "d"
    else:
        raise TypeError("Type (%s) not recognized" % type(fmt[1]))
    return fmt
