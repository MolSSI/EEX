"""
Metadata for four-body terms listed in alphabetical order

Each functional_form has the following values:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the four-body term style

This data will be validated during testing.
"""

_four_body_functional_forms = {
    "harmonic": {
        "form": "K * (theta - theta0) ** 2",
        "parameters": ["K", "theta0"],
        "units": {
            "K": "[energy] [arcunits] ** -2",
            "theta0": "[arcunits]"
        },
        "description": "A harmonic dihedral"
    },
}

### Do NOT edit below this line

four_body_metadata = {}

# Valid variables used in all four-body terms
four_body_metadata["variables"] = {
    "theta": {
        "units": "[arcunits]",
        "description": "The dihedral between the indexed atom1, atom2, atom3, and atom4.",
    },
}

# Add store data
four_body_metadata["store_name"] = "4body"
four_body_metadata["store_indices"] = {
    "atom1": "Index of the first atom.",
    "atom2": "Index of the second atom.",
    "atom3": "Index of the third atom.",
    "atom4": "Index of the fourth atom.",
    "term_index": "Index of four_body_type stored in the DataLayer.",
}

four_body_metadata["index_columns"] = ["atom1", "atom2", "atom3", "atom4"]
four_body_metadata["forms"] = _four_body_functional_forms
