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
    "class2": {
        "form": "K2 * (r-r0) ** 2 + K3 * (r-r0) ** 3 + K4 * (r-r0) ** 4",
        "parameters": ["R0", "K2", "K3", "K4"],
        "units": {
            "R0": "[length]",
            "K2": "[energy] [length] ** -2",
            "K3": "[energy] [length] ** -3",
            "K4": "[energy] [length] ** -4"
        },
        "description": "This is a class2 bond"
    },
    "fene": {
        "form": "-0.5*K*R0 ** 2*ln(1-(r/R0) ** 2) + 4*epsilon((sigma/r) ** 12 - (sigma/r) ** 6)) + epsilon",
        "parameters": ["K", "R0", "epsilon", "sigma"],
        "units": {
            "K": "[energy] [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]"
        },
        "description": "This is a fene bond!"
    },
    "fene/expand": {
        "form":
        "-0.5*K*R0 ** 2*ln(1-((r-delta)/R0) ** 2 + 4*epsilon((sigma/(r-delta)) ** 12 - (sigma/(r-delta)) ** 6)) + epsilon",
        "parameters": ["K", "R0", "epsilon", "sigma", "delta"],
        "units": {
            "K": "[energy] [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]",
            "delta": "[length]"
        },
        "description":
        "This is fene/expand bond"
    },
    # "hybrid": {},  # Special case - allows for more than one bond type in a simulation
    "harmonic": {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "[energy] [length] ** -2",
            "R0": "[length]"
        },
        "description": "This is a harmonic bond"
    },
    "morse": {
        "form": "D * (1 - e ** (-alpha * (r-R0))) ** 2",
        "parameters": ["D", "alpha", "R0"],
        "units": {
            "D": "[energy]",
            "alpha": "[length] ** -1",
            "R0": "[length]"
        },
        "description": "This is a class2 bond"
    },
    "nonlinear": {
        "form": "(epsilon*(r-R0) ** 2) / (lambda ** 2-(r-R0) ** 2)",
        "parameters": ["epsilon", "R0", "lambda"],
        "units": {
            "epsilon": "[energy]",
            "R0": "[length]",
            "lambda": "[length]"
        },
        "description": "This is a nonlinear bond"
    },
    # "table": {},  # Special case - creation of interpolation tables.
    "quartic": {
        "form": "K(r-Rc) ** 2 * (r-Rc-B1)*(r-Rc-B2) + U0 + 4*epsilon*((sigma/r) ** 12 - (sigma/r) ** 6) + epsilon",
        "parameters": ["K", "B1", "B2", "Rc", "U0"],
        "units": {
            "K": "[energy] [length] ** -4",
            "B1": "[length]",
            "B2": "[length]",
            "Rc": "[length]",
            "U0": "[energy]"
        },
        "description": "This is a quartic bond"
    },
}

### Please do not edit below this line

# Internal store name
two_body_metadata = {}

# Valid variables used in all two-body terms
two_body_metadata["variables"] = {"r": {"units": "[length]", "description": "Distance between the two indexed atoms."}}

# Add store data
two_body_metadata["store_name"] = "2body"
two_body_metadata["store_indices"] = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "two_body_type": "Index of two_body_type stored in the DataLayer.",
}
two_body_metadata["forms"] = _two_body_functional_forms
