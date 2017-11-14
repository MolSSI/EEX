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
