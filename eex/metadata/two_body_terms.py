"""
Metadata for two-body terms
"""

from collections import OrderedDict

# Internal store name
_two_body_store_name = "2body"

# Valid variables used in all two-body terms
_two_body_variables = {"r": {"units": "[distance]", "description": "Distance between the two indexed atoms."}}

# Valid columns of the store
_two_body_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "two_body_type": "Index of two_body_type stored in the DataLayer.",
}

# Two-body styles listed in alphabetical order
# Each style has the following values:
#   - form: The overall mathematical expression of the term
#   - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
#   - description: A short word description of the two-body term style
_two_body_styles = {
    "none": {},
    "zero": {},
    "class2": {
        "form":
        "K2 * (r-r0) ** 2 + K3 * (r-r0) ** 3 + K4 * (r-r0) ** 4",
        "parameters":
        OrderedDict({
            "R0": "[distance]",
            "K2": "[energy] [distance] ** -2",
            "K3": "[energy] [distance] ** -3",
            "K4": "[energy] [distance] ** -4"
        }),
        "description":
        "This is a class2 bond"
    },
    "fene": {
        "form":
        "-0.5*K*R0 ** 2*ln(1-(r/R0) ** 2) + 4*epsilon((sigma/r) ** 12 - (sigma/r) ** 6)) + epsilon",
        "parameters":
        OrderedDict({
            "K": "[energy] [distance] ** -2",
            "R0": "[distance]",
            "epsilon": "[energy]",
            "sigma": "[distance]"
        }),
        "description":
        "This is a fene bond!"
    },
    "fene/expand": {
        "form":
        "-0.5*K*R0 ** 2*ln(1-((r-delta)/R0) ** 2 + 4*epsilon((sigma/(r-delta)) ** 12 - (sigma/(r-delta)) ** 6)) + epsilon",
        "parameters":
        OrderedDict({
            "K": "[energy] [distance] ** -2",
            "r0": "[distance]",
            "epsilon": "[energy]",
            "sigma": "[distance]",
            "delta": "[distance]"
        }),
        "description":
        "This is fene/expand bond"
    },
    "hybrid": {},  # Special case - allows for more than one bond type in a simulation
    "harmonic": {
        "form": "K*(r-r0) ** 2",
        "parameters": OrderedDict({
            "K": "[energy] [distance] ** -2",
            "r0": "[distance]"
        }),
        "description": "This is a harmonic bond"
    },
    "morse": {
        "form": "D * (1 - e ** (-alpha * (r-r0))) ** 2",
        "parameters": OrderedDict({
            "D": "[energy]",
            "alpha": "[distance] ** -1",
            "r0": "[distance]"
        }),
        "description": "This is a class2 bond"
    },
    "nonlinear": {
        "form": "(epsilon*(r-r0) ** 2) / (lambda ** 2-(r-r0) ** 2)",
        "parameters": OrderedDict({
            "epsilon": "[energy]",
            "r0": "[distance]",
            "lambda": "[distance]"
        }),
        "description": "This is a nonlinear bond"
    },
    "table": {},  # Special case - creation of interpolation tables.
    "quartic": {
        "form":
        "K(r-Rc) ** 2 * (r-Rc-B1)*(r-Rc-B2) + U0 + 4*epsilon*((sigma/r) ** 12 - (sigma/r) ** 6) + epsilon",
        "parameters":
        OrderedDict({
            "K": "[energy] [distance] ** -4",
            "B1": "[distance]",
            "B2": "[distance]",
            "Rc": "[distance]",
            "U0": "[energy]"
        }),
        "description":
        "This is a quartic bond"
    },
}
