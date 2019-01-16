# Steps for writing a plugin for EEX

This guide is currently *under construction*, but when completed, should explain EEX's reader/writer strategy
and contain instructions for writing reader/writer plugins.

## 1. Build forcefield files

**filename**: `program_ff.py`

**examples**:
- [amber_ff.py](./amber/amber_ff.py)
- [lammps_ff.py](./lammps/lammps_ff.py)

where `program` refers to the software pacakge for which you are developing a reader/writer plugin.

This file contains dictionaries for two body (bonds), three body (angles), four body (dihedrals), and nonbonded terms which are valid for the software package.
Bonded terms (i.e. bonds, angles, and dihedrals) should be a dictionary named using the convention `_n_body_functional_forms`.

Each the keys of each dictionary are keywords (names) for functional forms which are valid in the program. The keywords for each
term dictionary (`_two_body_functional_forms`, etc), should match (excluding constant coefficients) the keyword used in EEX's
internal representation. Valid term forms can be found in `metadata`.

For bonded terms, each key should have the following sub-keys:
  - form: The overall mathematical expression of the term
  - parameters: The ordered name of the terms as they will be stored with their expected unit contexts
  - units: A dictionary that contains the unit context for each parameter
  - description: A short word description of the two-body term style

For example, the entry for two body (bonded) terms in Amber is shown here:

```
_two_body_functional_forms = {
    "harmonic": {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "kcal * mol ** -1 angstrom ** -2",
            "R0": "angstrom"
        },
        "description": "This is an amber harmonic bond"
    },
}
```

The bottom of the file builds the dictionary `term_data`, where the keys are the order of the term. This metadata will be used first to perform compatibility checks,
to see if a particular translation can be performed, then in the reader and writer plugins.

In the example above (the amber force field metadata), the units are always the same for a given functional
form. However, in some instances (such as LAMMPS), it may be necessary to specify the general units and build the units used based on keywords.

For example, the LAMMPS entry for a harmonic bond does not contain specific units.

```
    "harmonic": {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "[energy] * [length] ** -2",
            "R0": "[length]"
        },
        "description": "This is a harmonic bond"
    },
```

There is an additional level where the preferred style is set for nonbonded terms.

```
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
    },
}
```

TODO - program metadata validator

## 2. Build list of required metadata

