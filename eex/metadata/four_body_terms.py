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
            "n": "[]",  # must be type int, no units
            "d":
            "[arcunits]",  # must be type int. Differs because units must be [degrees] regardless of units command ?
        },
        "description": "This is a charmm dihedral"
    },
    "multi/harmonic": {
        "form": "Sum(A_n * (cos(phi))**(n-1))",  #n goes from 1 to 5
        "parameters": ["A_n"],
        "units": {
            "A_n": "[energy]",
        },
        "description": "This is a multi/harmonic dihedral"
    },
    "fourier": {
        "form": "Sum(k_i * (1.0 + cos(n_i * phi - d_i)))",
        "parameters": ["k_i", "n_i", "d_i"],
        "units": {
            "k_i": "[energy]",
            "n_i": "[]",
            "d_i": "[arcunits]"
        },
        "description": "This is a fourier dihedral"
    },
    "harmonic": {
        "form": "K * (1 + d * cos( n * phi))",
        "parameters": ["K", "n", "d"],
        "units": {
            "K": "[energy]",
            "n": "[]", # + / 1
            "d": "[]"
        },
        "description": "A harmonic dihedral"
    },
    "helix": {
        "form": "A*(1-cos(phi)) + B*(1+cos(3*phi)) + C*(1+cos(phi+pi/4))",
        "parameters": ["A", "B", "C"],
        "units": {
            "A": "[energy]",
            "B": "[energy]",
            "C": "[energy]"
        },
        "description": "This is a helix dihedral"
    },
    "nharmonic": {
        "form": "Sum( A_n*(cos(phi))**(n-1))",
        "parameters": ["A_n", "n"],
        "units": {
            "A_n": "[energy]",
            "n": "[]"
        },
        "description": "This is a nharmonic dihedral"
    },
    "opls": {
        "form": "0.5*K_1*(1+cos(phi))+0.5*K_2*(1-cos(2*phi))+0.5*K_3+(1+cos(3*phi))+0.5*K_4*(1-cos(4*phi))",
        "parameters": ["K_1", "K_2", "K_3", "K_4"],
        "units": {
            "K_1": "[energy]",
            "K_2": "[energy]",
            "K_3": "[energy]",
            "K_4": "[energy]",
        },
        "description": "This is a opls dihedral"
    },
    "quadratic": {
        "form": "K*(phi-phi0)**2",
        "parameters": ["K", "phi0"],
        "units": {
            "K": "[energy]",
            "phi0": "[arcunits]"
        },
        "description": "This is a quadratic dihedral"
    },
    #######IMPROPERS START HERE
    #"class2_improper": {},
    #"hybrid_improper": {},
    #"none_improper": {},
    #"ring_improper": {},
    #"umbrella_improper": {},
    #"zero_improper": {},
    "cossq": {
        "form": "0.5*K*(cos(chi-chi0))**2",
        "parameters": ["K", "chi0"],
        "units": {
            "K": "[energy]",
            "chi0": "[arcunits]"
        },
        "description": "This is a cossq improper"
    },
    "cvff": {
        "form": "K*(1+d*cos(n*chi))",
        "parameters": ["K", "d", "n"],
        "units": {
            "K": "[energy]",
            "d": "[]",
            "n": "[]"
        },
        "description": "This is a cvff improper"
    },
    "distance": {
        "form": "K_2*r**2+K_4*r**4",
        "parameters": ["K_2", "K_4"],
        "units": {
            "K_2": "[energy] [length]**2",
            "K_4": "[energy] [length]**4",
        },
        "description": "This is a distance improper"
    },
    "fourier_improper": {
        "form": "K*(C0+C1*cos(omega)+C2*cos(2*omega))",
        "parameters": ["K", "C0", "C1", "C2"],
        "units": {
            "K": "[energy]",
            "C0": "[]",
            "C1": "[]",
            "C2": "[]",
        },
        "description": "This is a fourier improper"
    },
    "harmonic": {
        "form": "K*(chi-chi0)**2",
        "parameters": ["K", "chi0"],
        "units": {
            "K": "[energy]",
            "chi0": "[arcunits]",
        },
        "description": "This is a harmonic improper"
    },
}

### Do NOT edit below this line

four_body_metadata = {}

# Valid variables used in all four-body terms
four_body_metadata["variables"] = {
    "phi": {
        "units": "[arcunits]",
        "description": "The dihedral between the indexed atom1, atom2, atom3, and atom4.",
    },
    "theta": {
        "units": "[arcunits]",
        "description": "Angle between three given atoms"
    },
    "r": {
        "units": "[distance]",
        "description": "Distance between two given atoms"
    },
    "chi": {
        "units": "[arcunits]",
        "description": "Improper angle"
    },
    "omega": {
        "units":
        "[arcunits]",
        "description":
        "Angle between the vector formed by a non-central atom and the plane formed by the other three atoms"
    }
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
