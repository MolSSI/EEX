# Steps for writing a plugin for EEX

This guide is currently *under construction*, but when completed, should explain EEX's reader/writer strategy
and contain instructions for writing reader/writer plugins.

## 1. Build forcefield dictionary.

filename: `program_ff.py`

where `program` refers to the software pacakge for which you are developing a reader/writer plugin.

The general format for this file is given below. This file contains dictionaries for two body,
three body, four body, and nonbonded terms which are valid for the software package. The bottom of the file builds another
dictionary, `term_data`, where the keys are the order of the term. This metadata will be used first to perform compatibility checks,
to see if a particular translation can be performed, then in the reader and writer plugins.

The keywords for each term dictionary (`_two_body_functional_forms`, etc), should match (excluding constant coefficients) the keyword used in EEX's
internal representation. Valid term forms can be found in `metadata`.

```
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
```

In the example above (the amber force field metadata), the units are always the same for a given functional
form. However, in some instances (such as LAMMPS), it may be necessary to specify the general units and build the units used
in

## 2. Build list of required metadata

