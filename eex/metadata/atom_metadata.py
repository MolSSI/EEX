"""
A file containing the valid categories and fields found in energy expressions.
"""

### A list of all valid atom properties that EEX can consume
atom_metadata = {
    "molecule_index": {
        "required_columns": ["molecule_index"],
        "description": "The associated molecular index for a given atom.",
        "unique": True,
        "dtype": int,
        "units": None,
    },
    "atom_name": {
        "required_columns": ["atom_name"],
        "description": "The unique name of a given atom.",
        "unique": True,
        "dtype": str,
        "units": None,
    },
    "atom_type": {
        "required_columns": ["atom_type"],
        "description": "The integer index of the atom type.",
        "unique": True,
        "dtype": int,
        "units": None,
    },
    "atomic_number": {
        "required_columns": ["atomic_number"],
        "description": "The atomic number of a given atom.",
        "unique": True,
        "dtype": int,
        "units": None,
    },
    "charge": {
        "required_columns": ["charge"],
        "description": "The charge of a given atom.",
        "unique": False,
        "dtype": float,
        "units": {
            "charge": "[charge]"
        },
        "tol": 8,
    },
    "xyz": {
        "required_columns": ["X", "Y", "Z"],
        "description": "The XYZ coordinates of a given atom.",
        "unique": True,
        "dtype": float,
        "units": {
            "X": "[length]",
            "Y": "[length]",
            "Z": "[length]"
        },
        "tol": 8,
    },
    "mass": {
        "required_columns": ["mass"],
        "description": "The mass of the atom.",
        "unique": False,
        "dtype": float,
        "units": {
            "mass": "[mass] / [substance]"
        },
        "tol": 8,
    },
    "residue_index": {
        "required_columns": ["residue_index"],
        "description": "The associated residue of a given atom.",
        "unique": True,
        "dtype": int,
        "units": None,
    },
    "residue_name": {
        "required_columns": ["residue_name"],
        "description": "The name of the associated residue of a given atom.",
        "unique": False,
        "dtype": str,
        "units": None,
    },
    "lj_coeffs": {
        "required_columns": ["LJ A", "LJ B"],
        "description": "The LJ parameters in their A/B form",
        "unique": False,
        "dtype": float,
        "units": {
            "LJ A": "[energy] * [length] ** -12",
            "LJ B": "[energy] * [length] ** -6"
        },
        "tol": 8,
    }
    # "implicit_solvent_radius":
    # ["implicit_solvent_radius"],  # Amber radii section - used for implicit solvent calculations
}

### Build required dictionaries
atom_property_to_column = {k: v["required_columns"] for k, v in atom_metadata.items()}

# Lookup dictionaries
_temp_symbol = [
    "X", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y",
    "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE",
    "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR",
    "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM",
    "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "UUB", "UUT", "UUQ",
    "UUP", "UUH", "UUS", "UUO"
]

atom_number_to_symbol = {k: v for k, v in zip(range(108), _temp_symbol)}
atom_symbol_to_number = {k: v for v, k in atom_number_to_symbol.items()}
