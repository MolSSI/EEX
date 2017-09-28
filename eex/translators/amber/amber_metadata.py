
# Possible size keys to look for in the header
size_keys = [
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH", "MPHIA", "NHPARM", "NPARM", "NNB", "NRES",
    "NBONA", "NTHETA", "NPHIA", "NUMBND", "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT", "NBPER", "NGPER", "NDPER",
    "MBPER", "MGPER", "MDPER", "IFBOX", "NMXRS", "IFCAP", "NUMEXTRA"
]

data_labels = {
    "ATOM_NAME": ["NATOM"],
    "CHARGE": ["NATOM"],
    "ATOMIC_NUMBER": ["NATOM"],
    "MASS": ["NATOM"],
    "ATOM_TYPE_INDEX": ["NATOM"],
    "NUMBER_EXCLUDED_ATOMS": ["NATOM"],
    "NONBONDED_PARM_INDEX": ["NTYPES ** 2"],
    "RESIDUE_LABEL": ["NRES"],
    "RESIDUE_POINTER": ["NRES"],
    "BOND_FORCE_CONSTANT": ["NUMBND"],
    "BOND_EQUIL_VALUE": ["NUMBND"],
    "ANGLE_FORCE_CONSTANT": ["NUMANG"],
    "ANGLE_EQUIL_VALUE": ["NUMANG"],
    "DIHEDRAL_FORCE_CONSTANT": ["NPTRA"],
    "DIHEDRAL_PERIODICITY": ["NPTRA"],
    "DIHEDRAL_PHASE": ["NPTRA"],
    "SCEE_SCALE_FACTOR": ["NPTRA"],
    "SCNB_SCALE_FACTOR": ["NPTRA"],
    "SOLTY": ["NATYP"],
    "LENNARD_JONES_ACOEF": ["(NTYPES * (NTYPES + 1)) / 2"],
    "LENNARD_JONES_BCOEF": ["(NTYPES * (NTYPES + 1)) / 2"],
    "BONDS_INC_HYDROGEN": ["3 * NBONH"],
    "BONDS_WITHOUT_HYDROGEN": ["3 * NBONA"],
    "ANGLES_INC_HYDROGEN": ["4 * NTHETH"],
    "ANGLES_WITHOUT_HYDROGEN": ["4 * NTHETA"],
    "DIHEDRALS_INC_HYDROGEN": ["5 * NPHIH"],
    "DIHEDRALS_WITHOUT_HYDROGEN": ["5 * NPHIA"],
    "EXCLUDED_ATOMS_LIST": ["NNB"],
    "HBOND_ACOEF": ["NPHB"],
    "HBOND_BCOEF": ["NPHB"],
    "HBCUT": ["NPHB"],
    "AMBER_ATOM_TYPE": ["NATOM"],
    "TREE_CHAIN_CLASSIFICATION": ["NATOM"],
    "JOIN_ARRAY": ["NATOM"],  # This section is filled with zeros, but is unused. We should not store it.
    "IROTAT": ["NATOM"],
    "SOLVENT_POINTERS": ["3 if IFBOX else 0"],
    "ATOMS_PER_MOLECULE": ["NATOM"],
    # "ATOMS_PER_MOLECULE": ["SOLVENT_POINTERS[1] if IFBOX else 0"], # SOLVENT_POINTERS[1] == NPSM
    "BOX_DIMENSIONS": [4],
    "CAP_INFO": ["1 if IFCAP else 0"],
    "CAP_INFO2": ["4 if IFCAP else 0"],
    "RADIUS_SET": [1],
    "RADII": ["NATOM"],
    "SCREEN": ["NATOM"],
    "IPOL": [1],
    "POLARIZABILITY": [0]
    # "POLARIZABILITY": ["NATOM if IPOL else 0"]
}

atom_data_units = {
    "charge": "18.2223 * e",  # Internal units
    "mass": "g * mol ** -1",
}

term_data_units = {
    # "BOND_FORCE_CONSTANT": ["kcal * mol ** -1 angstrom ** -2"],
    # "BOND_EQUIL_VALUE": ["angstrom"],
    # "ANGLE_FORCE_CONSTANT": ["kcal * mol ** -2 radian ** 2"],
    # "ANGLE_EQUIL_VALUE": ["radian"],
    "DIHEDRAL_FORCE_CONSTANT": ["kcal * mol ** -1"],
    "DIHEDRAL_PHASE": ["radian"],
    "LENNARD_JONES_ACOEFF": ["kcal * mol ** -12"],
    "LENNARD_JONES_BCOEFF": ["kcal * mol ** -6"],
}
register_forms = [
    (2, "harmonic", {"K": "kcal * mol ** -1 angstrom ** -2", "R0": "angstrom"}),
    (3, "harmonic", {"K": "kcal * mol ** -2 radian ** 2", "theta0": "radian"}),

]

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

topology_store_names = [
    "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN", "ANGLES_WITHOUT_HYDROGEN",
    "DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"
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
            "K": "kcal * mol ** -1 angstrom ** -2",
            "theta0": "angstrom"
        },
        "column_names": {
            "BOND_FORCE_CONSTANT": "K",
            "BOND_EQUIL_VALUE": "theta0"
        }
    }
    # "dihedral": ["DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PERIODICITY", "DIHEDRAL_PHASE"],
    # "non-bonded": ["LENNARD_JONES_ACOEFF", "LENNARD_JONES_BCOEFF"]
}

store_other = []
for k, v in forcefield_parameters.items():
    store_other.extend(list(v["column_names"]))
store_other.extend(residue_store_names)
