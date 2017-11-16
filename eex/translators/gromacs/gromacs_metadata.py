"""
A helper file to produce GROMACS metadata based on the type of computation
"""

import copy

from . import gromacs_ff 
import eex

data_keys = [
    "defaults", "atomtypes", "nonbond_params", "moleculetype", "atoms", 
    "bonds", "pairs", "angles", "dihedrals", "system", "molecules"
]

# Units based GROMACS unit styles (http://lammps.sandia.gov/doc/units.html)

bond_styles = { 
    "harmonic": {
        "order": 2,
        "form": "0.5 * K * (r - R0) ** 2",
        "units": {
            "K": "kjoule * mol ** -1 * nanometer ** -2",
            "R0": "nanometer"
        },  
    },  
    "fourth_power": {
        "order": 2,
        "form": "fourth_power",
        "units": {
            "K": "kjoule * mol ** -1 * nanometer ** -4",
            "R0": "nanometer"
        },  
    },  
    "morse": {
        "order": 2,
        "form": "morse",
        "units": {
            "D": "kjoule * mol ** -1",
            "alpha": "nanometer ** -1",
            "R0": "nanometer",
        },  
    },  
    "cubic_bond": {
        "order": 2,
        "form": "cubic_bond",
        "units": {
            "K": "kjoule * mol ** -1 * nanometer ** -2",
            "K_cub": "nanometer ** -1",
            "R0": "nanometer",
        },  
    },  
    "fene": {
        "order": 2,
        "form": "fene",
        "units": {
            "K": "kjoule * mol ** -1 * nanometer ** -2",
            "R0": "nanometer",
	    "epsilon": "kjoule * mol ** -1", # Note that this is always zero in GMX!
	    "sigma": "nanometer", # Note that this is always zero in GMX! 
        },  
    },  
}


#angles not included: 
# bond bond cross term
# bond angle cross term
angle_styles = { 
    "harmonic": {
        "order": 2,
        "form": "0.5 * K * (theta - theta0) ** 2",
        "units": {
            "K": "kjoule * mol ** -1 * radian ** -2",
            "R0": "degree"
        },  
    },  
    "cosine": {
        "order": 2,
        "form": "0.5*K*(cos(theta)-cos(theta0))**2",
        "units": {
            "K": "kjoule * mol ** -1 * radian ** -2",
            "theta0": "degree"
        },  
    },  
    "restricted": {
        "order": 2,
        "form": "restricted",
        "units": {
            "K": "kjoule * mol ** -1",
            "R0": "degree",
        },  
    },  
    "urey-bradley": {
        "order": 2,
        "form": "0.5*K*(theta-theta0)**2 + 0.5*K_ub*(r13-R_ub)**2",
        "units": {
            "K": "kjoule * mol ** -1 ** radian ** -2",
            "theta0": "degree",
            "K_ub": "kjoule * mol ** -1 ** nanometer ** -2",
            "R_ub": "nanometer",
        },  
    },  
    "quartic": {
        "order": 2,
        "form": "quartic_gmx",
        "units": {
            "K0": "kjoule * mol",
            "K1": "kjoule * mol ** -1",
            "K2": "kjoule * mol ** -2",
            "K3": "kjoule * mol ** -3",
            "K4": "kjoule * mol ** -4",
            "R0": "degree",
        },  
    },  
}


dihedral_types= { 
    "charmm": {
        "order": 4,
        "form": "charmmfsw",
        "units": {
            "K": "kjoule * mol ** -1",
            "n": "count",
            "d": "degree"
        },  
    },  
    "ryckaert-bellemans": {
        "order": 4,
        "form": "opls",
        "units": {
            "K_1": "kjoule * mol ** -1",
            "K_2": "kjoule * mol ** -1",
            "K_3": "kjoule * mol ** -1",
            "K_4": "kjoule * mol ** -1",
        },  
    },  
    "opls": {
        "order": 4,
        "form": "opls",
        "units": {
            "K_1": "kjoule * mol ** -1",
            "K_2": "kjoule * mol ** -1",
            "K_3": "kjoule * mol ** -1",
            "K_4": "kjoule * mol ** -1",
        },  
    },  
    "restricted": {
        "order": 4,
        "form": "restricted",
        "units": {
            "K": "kjoule * mol ** -1",
            "phi0": "degree",
        },  
    },  
}

improper_types= { 
    "harmonic": {
        "order": 4,
        "form": "0.5*K*(chi-chi0)**2",
        "units": {
            "K": "kjoule * mol ** -1",
            "chi0": "degree",
        },  
    },  
    "periodic": {
        "order": 4,
        "form": "charmmfsw", 
        "units": {
            "K": "kjoule * mol ** -1",
            "n": "count",
            "d": "degree"
        },  
    },  
}
