# Dictionaries for coeff format and units - http://lammps.sandia.gov/doc/Section_commands.html
# NYI
_bond_styles = {
    "none" : {},

    "fene": {
        "terms": ["K", "r0", "epsilon", "sigma"],
        "K" : "energy distance^-2",
        "r0": "distance",
        "epsilon": "energy",
        "sigma": "distance",
    },

    "nonlinear": {
        "terms": ["epsilon", "r0", "lambda"]
        "epsilon" : "energy",
        "r0" : "distance",
        "lambda": "distance"
    },

    "zero": {},

    "fene/expand": {
        "terms": ["K", "r0", "epsilon", "sigma", "delta"],
        "K": "energy distance^-2",
        "r0": "distance",
        "epsilon": "energy",
        "sigma": "distance",
        "delta": "distance,"
    },

    "quartic": {
        "terms": ["K", "B1", "B2", "Rc", "U0"],
        "K": "energy distance^-4",
        "B1": "distance",
        "B2": "distance",
        "Rc": "distance",
        "U0": "energy"
    },

    "hybrid": {}, # Special case - allows for more than one bond type in a simulation

    "harmonic": {
        "terms": ["K","r0"],
        "K" : "energy distance^-2",
        "r0" : "distance"
    },

    "table": {}, # Special case - creation of interpolation tables.

    "class2": {
        "terms": ["R0", "K2", "K3", "K4"]
        "R0": "distance",
        "K2": "energy distance^-2",
        "K3": "energy distance^-3",
        "K4": "energy distnace^-4",
    },

    "morse":{
        "terms": ["D", "alpha", "r0"],
        "D": "energy",
        "alpha": "distance^-1",
        "r0": "distance"
    },
}

_angle_styles = {
    "none": {},
    "class2": {
        "terms": ["theta0", "K2", "K3", "K4"],
        "theta0": "degree",
        "K2": "energy radian^-2",
        "K3": "energy radian^-3",
        "K4": "energy radian^-4,
    },

    "cosine/squared": {
        "terms": ["K", "theta0"],
        "K": "energy",
        "thetea0": "degrees",
    },

    "zero": {},

    "cosine": {
        "terms": ["K"],
        "K": "energy",
    },

    "harmonic": {
        "terms": ["K", "theta0"],
        "K": "energy radian^-2",
        "theta0": "degree"
    },

    "hybrid": {}, # Special case - allows for more than one angle type in a simulation

    "cosine/delta": {
        "terms": ["K", "theta0"],
        "K": "energy",
        "theta0": "degree"
    },

    "table": {}

    "charmm": {
        "terms": ["K", "theta0", "K_ub", "r_ub"],
        "K": "energy radian^2",
        "theta0": "degree",
        "K_ub": "energy distance^-2",
        "r_ub": "distance",
    },

    "cosine/periodic": {
        "terms": ["C", "B", "n"]
        "C": "energy",
        "B": "N/A",
        "n": "N/A"
    },
}

_dihedral_styles = {
    "none": {},
    "charmmfsw": {},
    "multi/harmonic": {},
    "zero": {},
    "class2": {},
    "opls": {},
    "hybrid": {},
    "harmonic": {},
    "charmm": {},
    "helix": {}
}

_improper_styles = {
    "none": {},
    "cvff": {},
    "zero": {},
    "harmonic": {},
    "hybrid": {},
    "umbrella": {},
    "class2": {},
}