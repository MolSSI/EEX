"""
Metadata for two-body terms listed in alphabetical order

Each style has the following values:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the two-body term style

"""

_two_body_functional_forms = {
    # "none": {},
    # "zero": {},
    # "nonlinear": {},
    # "hybrid": {},  # Special case - allows for more than one bond type in a simulation
    # "table": {},  # Special case - creation of interpolation tables.
    "class2": {
        "form": "K2 * (r-r0) ** 2 + K3 * (r-r0) ** 3 + K4 * (r-r0) ** 4",
        "parameters": ["r0", "K2", "K3", "K4"],
        "units": {
            "r0": "[length]",
            "K2": "[energy] * [length] ** -2",
            "K3": "[energy] * [length] ** -3",
            "K4": "[energy] * [length] ** -4"
        },
        "description": "This is a class2 bond",
        "canonical_form": "class2",
    },
    "fene": {
        "form": "-0.5*K*R0 ** 2 * log(1-(r/R0) ** 2) + 4*epsilon*((sigma/r) ** 12 - (sigma/r) ** 6) + epsilon",
        "parameters": ["K", "R0", "epsilon", "sigma"],
        "units": {
            "K": "[energy] * [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]"
        },
        "description": "This is a fene bond!",
        "canonical_form": "fene",
    },
    "fene/expand": {
        "form":
        "-0.5*K*R0 ** 2*log(1-((r-delta)/R0) ** 2 + 4*epsilon*((sigma/(r-delta)) ** 12 - (sigma/(r-delta)) ** 6)) + epsilon",
        "parameters": ["K", "R0", "epsilon", "sigma", "delta"],
        "units": {
            "K": "[energy] * [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]",
            "delta": "[length]"
        },
        "description": "This is fene/expand bond",
        "canonical_form": "fene/expand",
    },
    "harmonic": {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "[energy] * [length] ** -2",
            "R0": "[length]"
        },
        "description": "This is a harmonic bond",
        "canonical_form": "harmonic",
    },
    "morse": {
        "form": "D * (1 - exp(-alpha * (r-R0))) ** 2",
        "parameters": ["D", "alpha", "R0"],
        "units": {
            "D": "[energy]",
            "alpha": "[length] ** -1",
            "R0": "[length]"
        },
        "description": "This is a class2 bond",
        "canonical_form": "morse",
    },
    "nonlinear": {
        "form": "(epsilon * (r - R0) ** 2) / ((lam ** 2) - ((r - R0) ** 2))",
        "parameters": ["epsilon", "R0", "lam"],
        "units": {
            "epsilon": "[energy]",
            "R0": "[length]",
            "lam": "[length]"
        },
        "description": "This is a nonlinear bond",
        "canonical_form": "nonlinear",
    },
    "quartic": {
        "form": "K*(r-Rc) ** 2 * (r-Rc-B1)*(r-Rc-B2) + U0 + 4*epsilon*((sigma/r) ** 12 - (sigma/r) ** 6) + epsilon",
        "parameters": ["K", "B1", "B2", "Rc", "U0", "epsilon", "sigma"],
        "units": {
            "K": "[energy] * [length] ** -4",
            "B1": "[length]",
            "B2": "[length]",
            "Rc": "[length]",
            "U0": "[energy]",
            "epsilon": "[energy]",
            "sigma": "[energy]"
        },
        "description": "This is a quartic bond",
        "canonical_form": "quartic",
    },
    "harmonic/shift": {
        "form": "U_min/(R0 - R_c)**2 * ((r - R0) ** 2 - (R_c - R0) ** 2)",
        "parameters": ["U_min", "R0", "R_c"],
        "units": {
            "U_min": "[energy]",
            "R0": "[length]",
            "R_c": "[length]"
        },
        "description": "This is a harmonic/shift bond",
        "canonical_form": "harmonic/shift",
    },
    "oxdna/fene": {
        "form": "-epsilon / 2 * log(1 - ((r - R0) / delta)**2 ) ",
        "parameters": ["epsilon", "delta", "R0"],
        "units": {
            "epsilon": "[energy]",
            "delta": "[length]",
            "R0": "[length]"
        },
        "description": "This is a oxdna/fene bond",
        "canonical_form": "oxdna/fene",
    },
    "fourth_power": {  # Only found in gromacs
        "form": "0.25*K*(r**2 - R0**2)**2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "[energy] [length] ** -4",
            "R0": "[length]"
        },
        "description": "This is a fourth_power bond. Used in GROMOS force field",
        "canonical_form": "fourth_power",
    },
    "cubic_bond": {  # Only found in gromacs
        "form": "K*(r - R0)**2 + K * K_cub * (r-R0)**3",
        "parameters": ["K", "K_cub", "R0"],
        "units": {
            "K": "[energy] [length] ** -2",
            "K_cub": "[length] ** -1",
            "R0": "[length]"
        },
        "description": "This is a cubic_bond potential. Found in GROMACS",
        "canonical_form": "cubic_bond",
    },
}

# Please do not edit below this line

# Internal store name
two_body_metadata = {}

# Valid variables used in all two-body terms
two_body_metadata["variables"] = {"r": {"units": "[length]", "description": "Distance between the two indexed atoms."}}

_inverted_two_body_conversion = dict()
for k, v in _two_body_functional_forms.items():
    if v['canonical_form'] in _inverted_two_body_conversion.keys():
        _inverted_two_body_conversion[v['canonical_form']].append(k)
    else:
        _inverted_two_body_conversion[v['canonical_form']] = [k]

# Add store data
two_body_metadata["store_name"] = "2body"
two_body_metadata["store_indices"] = {
    "atom1": "Index of the first atom.",
    "atom2": "Index of the second atom.",
    "term_index": "Index of two_body_type stored in the DataLayer.",
}
two_body_metadata["index_columns"] = ["atom1", "atom2"]
two_body_metadata["forms"] = _two_body_functional_forms
two_body_metadata["group"] = _inverted_two_body_conversion
