"""
A file containing the valid categories and fields found in energy expressions.
"""

### A list of all valid atom properties that EEX can consume
atom_metadata = {
    "molecule_index": {
        "required_columns": ["molecule_index"],
        "description": "The associated molecular index for a given atom.",
        "unique": False,
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
        "unique": False,
        "dtype": int,
        "units": None,
    },
    "atomic_number": {
        "required_columns": ["atomic_number"],
        "description": "The atomic number of a given atom.",
        "unique": False,
        "dtype": int,
        "units": None,
    },
    "charge": {
        "required_columns": ["charge"],
        "description": "The charge of a given atom.",
        "unique": False,
        "dtype": int,
        "units": "[charge]",
        "tol": 8,
    },
    "xyz": {
        "required_columns": ["X", "Y", "Z"],
        "description": "The XYZ coordinates of a given atom.",
        "unique": True,
        "dtype": float,
        "units": "[length]",
        "tol": 8,
    },
    "mass": {
        "required_columns": ["mass"],
        "description": "The mass of the atom.",
        "unique": False,
        "dtype": float,
        "units": "[mass] / [substance]",
        "tol": 8,
    },
    "residue_index": {
        "required_columns": ["residue_index"],
        "description": "The associated residue of a given atom.",
        "unique": False,
        "dtype": int,
        "units": None,
    },
    "residue_name": {
        "required_columns": ["residue_name"],
        "description": "The name of the associated residue of a given atom.",
        "unique": False,
        "dtype": int,
        "units": None,
    },

    # "implicit_solvent_radius":
    # ["implicit_solvent_radius"],  # Amber radii section - used for implicit solvent calculations
}

### Build required dictionaries
atom_property_to_column = {k: v["required_columns"] for k, v in atom_metadata.items()}
