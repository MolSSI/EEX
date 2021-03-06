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
    # "none": {},
    # "zero": {},
    # "class2": {},
    # "cosine/shift/exp": {},
    # "hybrid": {},
    # "spherical": {},
    # "table": {},
    "charmmfsw": {
        "form": "K*(1 + cos(n*phi-d))",
        "parameters": ["K", "n", "d"],
        "units": {
            "K": "[energy]",
            "n": "count",
            "d": "[arcunit]",
        },
        "description": "This is a charmm dihedral",
        "canonical_form": "RB",
    },
    "multiharmonic": {
        "form":
        "A_1 + A_2 * (cos(phi)) + A_3 * (cos(phi)) ** 2 + A_4 * (cos(phi)) ** 3 + A_5 * (cos(phi)) ** (4)",
        "parameters": ["A_1", "A_2", "A_3", "A_4", "A_5"],
        "units": {
            "A_1": "[energy]",
            "A_2": "[energy]",
            "A_3": "[energy]",
            "A_4": "[energy]",
            "A_5": "[energy]",
        },
        "description":
        "This is a multi/harmonic dihedral",
        "canonical_form":
        "RB",
    },
    "RB": {
        "form":
        "A_0 + A_1 * (cos(phi)) + A_2 * (cos(phi)) ** 2 + A_3 * (cos(phi)) ** 3 + A_4 * (cos(phi)) ** (4) + A_5 * (cos(phi)) ** (5)",
        "parameters": ["A_0", "A_1", "A_2", "A_3", "A_4", "A_5"],
        "units": {
            "A_0": "[energy]",
            "A_1": "[energy]",
            "A_2": "[energy]",
            "A_3": "[energy]",
            "A_4": "[energy]",
            "A_5": "[energy]",
        },
        "description":
        "This is a ryckaert_bellemans",
        "canonical_form":
        "RB",
    },
    # "fourier": {
    #     "form": "Sum(k_i * (1.0 + cos(n_i * phi - d_i)))",
    #     "parameters": ["k_i", "n_i", "d_i"],
    #     "units": {
    #         "k_i": "[energy]",
    #         "n_i": "phase",
    #         "d_i": "[arcunit]"
    #     },
    #     "description": "This is a fourier dihedral"
    # },
    "harmonic": {
        "form": "K * (1 + d * cos(n * phi))",
        "parameters": ["K", "n", "d"],
        "units": {
            "K": "[energy]",
            "n": "phase",  # + / 1
            "d": "count"
        },
        "description": "A harmonic dihedral",
        "canonical_form": "RB",
    },
    "helix": {
        "form": "A*(1-cos(phi)) + B*(1+cos(3*phi)) + C*(1+cos(phi+PI/4))",
        "parameters": ["A", "B", "C"],
        "units": {
            "A": "[energy]",
            "B": "[energy]",
            "C": "[energy]"
        },
        "description": "This is a helix dihedral",
        "canonical_form": "helix",
    },
    # "nharmonic": {
    #     "form": "Sum( A_n*(cos(phi))**(n-1))",
    #     "parameters": ["A_n", "n"],
    #     "units": {
    #         "A_n": "[energy]",
    #         "n": "count"
    #     },
    #     "description": "This is a nharmonic dihedral"
    # },
    "opls": {
        "form":
        "0.5*K_1*(1+cos(phi))+0.5*K_2*(1-cos(2*phi))+0.5*K_3*(1+cos(3*phi))+0.5*K_4*(1-cos(4*phi))",
        "parameters": ["K_1", "K_2", "K_3", "K_4"],
        "units": {
            "K_1": "[energy]",
            "K_2": "[energy]",
            "K_3": "[energy]",
            "K_4": "[energy]",
        },
        "description":
        "This is a opls dihedral",
        "canonical_form":
        "RB",
    },
    "quadratic": {
        "form": "K*(phi-phi0)**2",
        "parameters": ["K", "phi0"],
        "units": {
            "K": "[energy]",
            "phi0": "[arcunit]"
        },
        "description": "This is a quadratic dihedral",
        "canonical_form": "quadratic",
    },
    "restricted": {
        "form": "0.5*K*(cos(phi)-cos(phi0))**2/(sin(phi)**2)",
        "parameters": ["K", "phi0"],
        "units": {
            "K": "[energy]",
            "phi0": "[arcunit]",
        },
        "description": "This is a restricted bending angle found in Gromacs",
        "canonical_form": "restricted",
    },

    # IMPROPERS START HERE
    #"class2_improper": {},
    #"hybrid_improper": {},
    #"none_improper": {},
    #"ring_improper": {},
    #"umbrella_improper": {},
    #"zero_improper": {},
    "cossq_improper": {
        "form": "0.5*K*(cos(chi-chi0))**2",
        "parameters": ["K", "chi0"],
        "units": {
            "K": "[energy]",
            "chi0": "[arcunit]"
        },
        "description": "This is a cossq improper",
        "canonical_form": "cossq_improper",
    },
    "cvff_improper": {
        "form": "K*(1+d*cos(n*chi))",
        "parameters": ["K", "d", "n"],
        "units": {
            "K": "[energy]",
            "d": "phase",
            "n": "count"
        },
        "description": "This is a cvff improper",
        "canonical_form": "RB",
    },
    "distance_improper": {
        "form": "K_2*r**2+K_4*r**4",
        "parameters": ["K_2", "K_4"],
        "units": {
            "K_2": "[energy] * [length]**2",
            "K_4": "[energy] * [length]**4",
        },
        "description": "This is a distance improper",
        "canonical_form": "distance_improper",
    },
    "fourier_improper": {
        "form": "K*(C0+C1*cos(omega)+C2*cos(2*omega))",
        "parameters": ["K", "C0", "C1", "C2"],
        "units": {
            "K": "[energy]",
            "C0": "dimensionless",
            "C1": "dimensionless",
            "C2": "dimensionless",
        },
        "description": "This is a fourier improper",
        "canonical_form": "RB",
    },
    "harmonic_improper": {
        "form": "K*(chi-chi0)**2",
        "parameters": ["K", "chi0"],
        "units": {
            "K": "[energy]",
            "chi0": "[arcunit]",
        },
        "description": "This is a harmonic improper",
        "canonical_form": "harmonic_improper",
    },
}

# Do NOT edit below this line

four_body_metadata = {}

# Valid variables used in all four-body terms
four_body_metadata["variables"] = {
    "phi": {
        "units":
        "[arcunit]",
        "description":
        "The dihedral between the indexed atom1, atom2, atom3, and atom4.",
    },
    "theta": {
        "units": "[arcunit]",
        "description": "Angle between three given atoms"
    },
    "r": {
        "units": "[distance]",
        "description": "Distance between two given atoms"
    },
    "chi": {
        "units": "[arcunit]",
        "description": "Improper angle"
    },
    "omega": {
        "units":
        "[arcunit]",
        "description":
        "Angle between the vector formed by a non-central atom and the plane formed by the other three atoms"
    }
}

_inverted_four_body_conversion = dict()
for k, v in _four_body_functional_forms.items():
    if v['canonical_form'] in _inverted_four_body_conversion.keys():
        _inverted_four_body_conversion[v['canonical_form']].append(k)
    else:
        _inverted_four_body_conversion[v['canonical_form']] = [k]

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
four_body_metadata["group"] = _inverted_four_body_conversion
