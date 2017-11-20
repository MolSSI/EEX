"""
A helper file to produce GROMACS metadata based on the type of computation
"""

data_keys = [
    "defaults", "atomtypes", "nonbond_params", "moleculetype", "atoms", 
    "bonds", "pairs", "angles", "dihedrals", "system", "molecules"
]

# Units based GROMACS unit styles (http://lammps.sandia.gov/doc/units.html)

#units = {
#    "[length]": "nanometer",
#    "[mass]": "u",
#    "[time]": "picosecond",
#    "[charge]": "e",
#    "[temperature]": "kelvin",
#    "[energy]": "(kilojoules * mol ** -1)",
#    "[force]": "(kilojoules * mol ** -1 nanometer ** -1)",
#    "[pressure]": "bar",
#    "[velocity]": "(nanometer picosecond** -1)",
#    "[dipole]": "e * nanometer",
#    "[electric field]": "(kilojoules * mol ** -1 * nanometer ** -1 * e ** -1)",
#    "[electric potential]": "(kilojoules * mol ** -1 * e ** -1)",
#    }

bond_styles = { 
    "harmonic": {
        "order": 2,
        "units": {
            "K": "0.5 * kilojoules * mol ** -1 * nanometer ** -2",
            "R0": "nanometer"
        },  
        "gmx_type": 6,
    },  
    "fourth_power": {
        "order": 2,
        "form": "fourth_power",
        "units": {
            "K": "kilojoules * mol ** -1 * nanometer ** -4",
            "R0": "nanometer"
        },  
        "gmx_type": 2,
    },  
    "morse": {
        "order": 2,
        "form": "morse",
        "units": {
            "D": "kilojoules * mol ** -1",
            "alpha": "nanometer ** -1",
            "R0": "nanometer",
        },  
        "gmx_type": 3,
    },  
    "cubic_bond": {
        "order": 2,
        "form": "cubic_bond",
        "units": {
            "K": "kilojoules * mol ** -1 * nanometer ** -2",
            "K_cub": "nanometer ** -1",
            "R0": "nanometer",
        },  
        "gmx_type": 4, 
    },  
    "fene": {
        "order": 2,
        "form": "fene",
        "units": {
            "K": "kilojoules * mol ** -1 * nanometer ** -2",
            "R0": "nanometer",
	    "epsilon": "kilojoules * mol ** -1", # Note that this is always zero in GMX!
	    "sigma": "nanometer", # Note that this is always zero in GMX! 
        },  
        "gmx_type": 7,
    },  
}


#angles not included: 
# bond bond cross term
# bond angle cross term
angle_styles = { 
    "harmonic": {
        "order": 2,
        "units": {
            "K": "0.5 * kilojoules * mol ** -1 * radian ** -2",
            "theta0": "degree"
        },  
    },  
    "cosine": {
        "order": 2,
        "units": {
            "K": "0.5 * kilojoules * mol ** -1 * radian ** -2",
            "theta0": "degree"
        },  
    },  
    "restricted": {
        "order": 2,
        "form": "restricted",
        "units": {
            "K": "kilojoules * mol ** -1",
            "theta0": "degree",
        },  
    },  
    "urey-bradley": {
        "order": 2,
        "units": {
            "K": "0.5 * kilojoules * mol ** -1 ** radian ** -2",
            "theta0": "degree",
            "K_ub": "0.5 * kilojoules * mol ** -1 ** nanometer ** -2",
            "R_ub": "nanometer",
        },  
    },  
    "quartic": {
        "order": 2,
        "form": "quartic_gmx",
        "units": {
            "K0": "kilojoules * mol",
            "K1": "kilojoules * mol ** -1",
            "K2": "kilojoules * mol ** -2",
            "K3": "kilojoules * mol ** -3",
            "K4": "kilojoules * mol ** -4",
            "theta0": "degree",
        },  
    },  
}


dihedral_styles = { 
    "charmm": {
        "order": 4,
        "form": "charmmfsw",
        "units": {
            "K": "kilojoules * mol ** -1",
            "n": "count",
            "d": "degree"
        },  
    },  
    "ryckaert-bellemans": {
        "order": 4,
        "form": "opls",
        "units": {
            "K_1": "kilojoules * mol ** -1",
            "K_2": "kilojoules * mol ** -1",
            "K_3": "kilojoules * mol ** -1",
            "K_4": "kilojoules * mol ** -1",
        },  
    },  
    "opls": {
        "order": 4,
        "form": "opls",
        "units": {
            "K_1": "kilojoules * mol ** -1",
            "K_2": "kilojoules * mol ** -1",
            "K_3": "kilojoules * mol ** -1",
            "K_4": "kilojoules * mol ** -1",
        },  
    },  
    "restricted": {
        "order": 4,
        "form": "restricted",
        "units": {
            "K": "kilojoules * mol ** -1",
            "phi0": "degree",
        },  
    },  
}

improper_styles = { 
    "harmonic": {
        "order": 4,
        "units": {
            "K": "0.5 * kilojoules * mol ** -1",
            "chi0": "degree",
        },  
    },  
    "periodic": {
        "order": 4,
        "form": "charmmfsw", 
        "units": {
            "K": "kilojoules * mol ** -1",
            "n": "count",
            "d": "degree"
        },  
    },  
}

term_data = {}
term_data[2] = bond_styles
term_data[3] = angle_styles
term_data[4] = dihedral_styles

