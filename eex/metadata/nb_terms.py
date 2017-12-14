"""
Contains all of the metadata for non-bonded terms.
"""

_nb_functional_forms = {
    "LJ": {
        "AB": {
            "form": "A/(r ** 12) - B/(r ** 6)",
            "parameters": ["A", "B"],
            "units": {
                "A": "[energy] * [length] ** 12",
                "B": "[energy] * [length] ** 6",
            },
            "description": "This is the AB LJ form",
        },
        "epsilon/sigma": {
            "form": "-4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)",
            "parameters": ["epsilon", "sigma"],
            "units": {
                "epsilon": "[energy]",
                "sigma": "[length]",
            },
            "description": "This is the classic LJ non-bonded",
        },
        "epsilon/Rmin": {
            "form": "-epsilon * ((Rmin / r) ** 12 - 2 * (Rmin / r) ** 6)",
            "parameters": ["epsilon", "Rmin"],
            "units": {
                "epsilon": "[energy]",
                "Rmin": "[length]",
            },
            "description": "This is the Rmin LJ form",
        },
    },
    "Buckingham": {
        "Buckingham": {
            "form": "A * exp(-r / rho) - C / r ** 6",
            "parameters": ["A", "rho", "C"],
            "units": {
                "A": "[energy]",
                "rho": "[length]",
                "C": "[energy] * [length] ** 6",
            },
            "description": "This is Buckingham form",
        },
    },
}

_defaults = {
    "LJ" : "AB",
    "Buckingham" : "Buckingham",
}

nb_metadata = {}

# Valid variables used in all two-body terms
nb_metadata["variables"] = {"r": {"units": "[length]", "description": "Distance between the two indexed atoms."}}

# Add store data
nb_metadata["store_name"] = "nb"
nb_metadata["store_indices"] = {
    "atom_type1": "Index of the first atom.",
    "atom_type2": "Index of the second atom.",
    # "term_index": "Index of nb_type stored in the DataLayer.",
}
# nb_metadata["index_columns"] = ["atom1", "atom2"]

nb_metadata["forms"] = _nb_functional_forms
nb_metadata["defaults"] = _defaults
