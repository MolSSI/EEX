"""
Metadata for three-body terms listed in alphabetical order

Each functional_form has the following values:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the three-body term style

This data will be validated during testing.
"""

_three_body_functional_forms = {
    # "none": {},
    # "class2": {},
    # "zero": {},
    # "hybrid": {},  # Special case - allows for more than one angle type in a simulation
    # "table": {},
    # "dipole": {},
    # "fourier": {},
    # "sdk": {},
    # "cosine/shift": {},
    "cosine/squared": {
        "form": "K*(cos(theta)-cos(theta0))**2",
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy]",
            "theta0": "[arcunit]"
        },
        "description": "This is a cosine/squared angle"
    },
    "cosine": {
        "form": "K*(1+cos(theta))",
        "parameters": ["K"],
        "units": {
            "K": "[energy]"
        },
        "description": "This is a cosine potential"
    },
    "harmonic": {
        "form": "K * (theta - theta0) ** 2",
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy] * [arcunit] ** -2",
            "theta0": "[arcunit]"
        },
        "description": "A harmonic angle"
    },
    "cosine/delta": {
        "form": "K*(1+cos(theta-theta0))",
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy]",
            "theta0": "[arcunit]"
        },
        "description": "This is a cosine/delta potential"
    },
    "charmm": {
        "form": "K*(theta-theta0)**2 + K_ub*(r13-R_ub)**2",
        "parameters": ["K", "theta0", "K_ub", "R_ub"],
        "units": {
            "K": "[energy] * [arcunit] ** 2",
            "theta0": "[arcunit]",
            "K_ub": "[energy] * [length] ** -2",
            "R_ub": "[length]",
        },
        "description": "A CHARMM angle term?"
    },
    "cosine/periodic": {
        "form": "C * (1-B*((-1)**n) * cos(n*theta))",
        "parameters": ["C", "B", "n"],
        "units": {
            "C": "[energy]",
            "B": "phase",  #1 or -1
            "n": "count"  # 1 2 3 4 5 or 6
        },
        "description": "This is a cosine/periodic potential"
    },
    "fourier": {
        "form": "K*(c0+c1*cos(theta)+c2*cos(2*theta))",
        "parameters": ["K", "c0", "c1", "c2"],
        "units": {
            "K": "[energy]",
            "c0": "dimensionless",
            "c1": "dimensionless",
            "c2": "dimensionless"
        },
        "description": "This is a fourier potential"
    },
    "quartic": {
        "form": "K2*(theta-theta0)**2+K3*(theta-theta0)**3+K4*(theta-theta0)**4",
        "parameters": ["K2", "K3", "K4", "theta0"],
        "units": {
            "K2": "[energy] * [arcunit]**-2",  #Lammps uses radians
            "K3": "[energy] * [arcunit]**-3",
            "K4": "[energy] * [arcunit]**-4",
            "theta0": "[arcunit]"  #Lammps converts this to radians
        },
        "description": "This is a quartic angle"
    },
    "cosine/shift": {
        "form": "-U_min / 2 * (1 + cos(theta - theta0))",
        "parameters": ["U_min", "theta0"],
        "units": {
            "U_min": "[energy]",
            "theta0": "[arcunit]"
        },
        "description": "This is a cosine/shift angle"
    },
    "fourier/simple": {
        "form": "K*(1 + c*cos(n*theta))",
        "parameters": ["K", "c", "n"],
        "units": {
            "K": "[energy]",
            "c": "dimensionless",
            "n": "phase",
        },
        "description": "This is a fourier/simple angle"
    },

    #"class2": {
    #    "form": "NYI",
    #    "parameters": ["theta0", "K2", "K3", "K4"],
    #    "units": {
    #        "theta0": "[arcunit]",
    #        "K2": "[energy] * [arcunit] ** -2",
    #        "K3": "[energy] * [arcunit] ** -3",
    #        "K4": "[energy] * [arcunit] ** -4",
    #    },
    #    "description": "A generalized class 2 force field."
    #},
    # "cosine": {
    #     "parameters": ["K"],
    #     "units": {
    #         "K": "[energy]",
    #     }
    # },
    # "cosine/delta": {
    #     "parameters": ["K", "theta0"],
    #     "units": {
    #         "K": "[energy]",
    #         "theta0": "[arcunit]"
    #     }
    # },
    # "cosine/periodic": {
    #     "parameters": ["C", "B", "n"],
    #     "units": {
    #         "C": "[energy]",
    #         "B": "N/A",
    #         "n": "N/A"
    #     }
    # },
    # "cosine/squared": {
    #     "parameters": ["K", "theta0"],
    #     "units": {
    #         "K": "[energy]",
    #         "theta0": "[arcunit]",
    #     }
    # },
}

### Do NOT edit below this line

three_body_metadata = {}

# Valid variables used in all three-body terms
three_body_metadata["variables"] = {
    "theta": {
        "units": "[arcunit]",
        "description": "The angle between the indexed atom1, atom2, and atom3.",
    },
    "r12": {
        "units": "[length]",
        "description": "The length between the first and second indexed atoms."
    },
    "r23": {
        "units": "[length]",
        "description": "The length between the second and third indexed atoms."
    },
    "r13": {
        "units": "[length]",
        "description": "The length between the first and third indexed atoms."
    },
}

# Add store data
three_body_metadata["store_name"] = "3body"
three_body_metadata["store_indices"] = {
    "atom1": "Index of the first atom.",
    "atom2": "Index of the second atom.",
    "atom3": "Index of the third atom.",
    "term_index": "Index of three_body_type stored in the DataLayer.",
}

three_body_metadata["index_columns"] = ["atom1", "atom2", "atom3"]
three_body_metadata["forms"] = _three_body_functional_forms
