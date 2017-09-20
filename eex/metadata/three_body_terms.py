"""
Metadata for three-body terms listed in alphabetical order

Each term_style has the following values:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the three-body term style

This data will be validated during testing.
"""

_three_body_term_styles = {
    "none": {},
    "zero": {},
    "hybrid": {},  # Special case - allows for more than one angle type in a simulation
    "table": {},
    "charmm": {
        "parameters": ["K", "theta0", "K_ub", "r_ub"],
        "units": {
            "K": "[energy] [arcunits] ** 2",
            "theta0": "[arcunits]",
            "K_ub": "[energy] [distance] ** -2",
            "r_ub": "[distance]",
        }
    },
    "class2": {
        "parameters": ["theta0", "K2", "K3", "K4"],
        "units": {
            "theta0": "[arcunits]",
            "K2": "[energy] [arcunits] ** -2",
            "K3": "[energy] [arcunits] ** -3",
            "K4": "[energy] [arcunits] ** -4",
        }
    },
    "cosine": {
        "parameters": ["K"],
        "units": {
            "K": "[energy]",
        }
    },
    "cosine/delta": {
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy]",
            "theta0": "[arcunits]"
        }
    },
    "cosine/periodic": {
        "parameters": ["C", "B", "n"],
        "units": {
            "C": "[energy]",
            "B": "N/A",
            "n": "N/A"
        }
    },
    "cosine/squared": {
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy]",
            "theta0": "[arcunits]",
        }
    },
    "harmonic": {
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy] [arcunits] ** -2",
            "theta0": "[arcunits]"
        }
    },
}

### Do NOT edit below this line

three_body_metadata = {}
# Internal store name
three_body_metadata["store_name"] = "3body"

# Valid variables used in all three-body terms
three_body_metadata["variables"] = {
    "theta": {
        "units": "[arcunits]",
        "description": "The angle between the indexed atom1, atom2, and atom3.",
    },
    "r12": {
        "units": "[distance]",
        "description": "The distance between the first and second indexed atoms."
    },
    "r23": {
        "units": "[distance]",
        "description": "The distance between the second and third indexed atoms."
    },
    "r13": {
        "units": "[distance]",
        "description": "The distance between the first and third indexed atoms."
    },
}

# Valid columns of the store
three_body_metadata["store_indices"] = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "atom_index3": "Index of the third atom.",
    "three_body_type": "Index of three_body_type stored in the DataLayer.",
}

three_body_metadata["styles"] = _three_body_term_styles