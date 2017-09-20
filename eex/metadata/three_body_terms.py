"""
Metadata for two-body terms
"""

# Internal store name
_three_body_store_name = "2body"

# Valid variables used in all three-body terms
_three_body_variables = {
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
_three_body_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "atom_index3": "Index of the third atom.",
    "three_body_type": "Index of three_body_type stored in the DataLayer.",
}

# Three-body styles listed in alphabetical order
# Each style has the following values:
#   - form: The overall mathematical expression of the term
#   - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
#   - description: A short word description of the two-body term style

_angle_styles = {
    "none": {},
    "class2": {
        "parameters": ["theta0", "K2", "K3", "K4"],
        "theta0": "[arcunits]",
        "K2": "[energy] [arcunits] ** -2",
        "K3": "[energy] [arcunits] ** -3",
        "K4": "[energy] [arcunits] ** -4",
    },
    "cosine/squared": {
        "parameters": ["K", "theta0"],
        "K": "[energy]",
        "thetea0": "[arcunits]s",
    },
    "zero": {},
    "cosine": {
        "parameters": ["K"],
        "K": "[energy]",
    },
    "harmonic": {
        "parameters": ["K", "theta0"],
        "K": "[energy] [arcunits] ** -2",
        "theta0": "[arcunits]"
    },
    "hybrid": {},  # Special case - allows for more than one angle type in a simulation
    "cosine/delta": {
        "parameters": ["K", "theta0"],
        "K": "[energy]",
        "theta0": "[arcunits]"
    },
    "table": {},
    "charmm": {
        "parameters": ["K", "theta0", "K_ub", "r_ub"],
        "K": "[energy] [arcunits] ** 2",
        "theta0": "[arcunits]",
        "K_ub": "[energy] [distance] ** -2",
        "r_ub": "[distance]",
    },
    "cosine/periodic": {
        "parameters": ["C", "B", "n"],
        "C": "[energy]",
        "B": "N/A",
        "n": "N/A"
    },
}
