"""
Metadata the amber forcefield.

Each style follows conventions outlined for EEX internal metadata.

Bonded terms (typically bonds, angles, and dihedrals) should be a dictionary named using the convention `_n_body_functional_forms`

Each dictionary has the has the following keys:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the two-body term style

Nonbonded terms are similar, except that an extra key sub-key (again, following the convention set by the internal EEX metadata.

"""
_two_body_functional_forms = {
    "harmonic": {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "kcal * mol ** -1 angstrom ** -2",
            "R0": "angstrom"
        },
        "description": "This is a harmonic bond"
    },
}

_three_body_functional_forms = {
    "harmonic": {
        "form": "K*(theta-theta0)**2",
        "parameters": ["K", "theta0"],
        "units": {
            "K": "kcal radian**-2",
            "theta0": "degree"
        },
        "description": "This is a harmonic angle"
    },
}

_four_body_functional_forms = {
    "charmmfsw": {
        "form": "K*(1 + cos(n*phi-d))",
        "parameters": ["K", "n", "d"],
        "units": {
            "K": "kcal * mol ** -1",
            "n": "phase",
            "d": "radians",
        },
        "description": "This is a charmm dihedral"
    },

}

_nonbond_functional_forms = {
    "LJ": {
        "AB": {
            "form": "A/(r ** 12.0) - B/(r ** 6.0)",
            "parameters": ["A", "B"],
            "units": {
                "A": "kcal * mol ** -1 * angstrom ** 12",
                "B": "kcal * mol ** -1 * angstrom ** 6",
            },
        },
    }
}


term_data = {}
term_data[2] = _two_body_functional_forms
term_data[3] = _three_body_functional_forms
term_data[4] = _four_body_functional_forms

