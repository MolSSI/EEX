"""
Dictionaries for coeff format and units - http://lammps.sandia.gov/doc/Section_commands.html
"""

from collections import OrderedDict

_valid_bond_variables = {
    "r": {
        "units": "[length]", 
        "description": "Distance between the two index atoms"
    },
}

_valid_angle_variables = {
    "theta": {
        "units": "[radian]",
        "description": "Angle between three consecutive atoms"
    },
    "r": {
        "units": "[length]",
        "description": "Distance between two given atoms"
    }
}

_valid_dihedral_variables = {
    "phi": {
        "units": "[radian]",
        "description": "Dihedral angle arising from four consecutive atoms"
    },
    "theta": {
        "units": "[radian]",
        "description": "Angle between three given atoms"
    },
    "r": {
        "units": "[distance]",
        "description": "Distance between two given atoms"
    }
}

_valid_improper_variables = {
    "chi": {
        "units": "[radian]",
        "description": "Improper angle"
    },
    "r": {
        "units": "[distance]",
        "description": "Distance between the central atom and the plane formed by the other three atoms"
    },
    "omega": {
        "units": "[radian]",
        "description": "Angle between the vector formed by a non-central atom and the plane formed by the other three atoms"
    }
}

# Valid columns of dataframe
_valid_bond_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",

    # DataLayers knows the bondtype
    "bond_type": "Index of bond_type stored in the DataLayer",

    # DataLayer needs to currate input and form unique bond_types
    "bond_style": "Bond style name from the _bond_styles dictionary",
    "coeffs": "..."
}

_valid_angle_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "atom_index3": "Index of the third atom.",

    # DataLayers knows the bondtype
    "angle_type": "Index of angle_type stored in the DataLayer",

    # DataLayer needs to currate input and form unique bond_types
    "angle_style": "Bond style name from the _bond_styles dictionary",
    "coeffs": "..."
}

_valid_dihedral_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "atom_index3": "Index of the third atom.",
    "atom_index4": "Index of the fourth atom.",

    # DataLayers knows the bondtype
    "dihedral_type": "Index of dihedral_type stored in the DataLayer",

    # DataLayer needs to currate input and form unique bond_types
    "dihedral_style": "improper style name from the _bond_styles dictionary",
    "coeffs": "..."
}

_valid_improper_indices = {
    "atom_index1": "Index of the first atom.",
    "atom_index2": "Index of the second atom.",
    "atom_index3": "Index of the third atom.",
    "atom_index4": "Index of the fourth atom.",

    # DataLayers knows the bondtype
    "improper_type": "Index of improper_type stored in the DataLayer",

    # DataLayer needs to currate input and form unique bond_types
    "improper_style": "Dihedral style name from the _bond_styles dictionary",
    "coeffs": "..."
}

_bond_styles = {
    "none": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "fene": {
        "form": "-0.5*K*R0 ** 2*ln(1-(r/R0) ** 2) + 4*epsilon((sigma/r) ** 12 - (sigma/r) ** 6)) + epsilon",
        "parameters": ["K","R0","epsilon","sigma"],
        "units": {
            "K": "[energy] [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]"
        },
        "description": "This is a fene bond!"
    },
    "nonlinear": {
        "form": "(epsilon*(r-r0) ** 2) / (lambda ** 2-(r-r0) ** 2)",
        "parameters": ["epsilon","r0","lambda"],
        "units": {
            "epsilon": "[energy]",
            "r0": "[length]",
            "lambda": "[length]"
        },
        "description": "This is a nonlinear bond"
    },
    "zero": {
        "form": "NYI",
        "parameters": "NYI":
        "units": "NYI",
        "description": "NYI"
    },
    "fene/expand": {
        "form": "-0.5*K*R0 ** 2*ln(1-((r-delta)/R0) ** 2 + 4*epsilon((sigma/(r-delta)) ** 12 - (sigma/(r-delta)) ** 6)) + epsilon",
        "parameters": ["K", "R0", "epsilon", "sigma", "delta"],
        "units": {
            "K": "[energy] [length] ** -2",
            "R0": "[length]",
            "epsilon": "[energy]",
            "sigma": "[length]",
            "delta": "[length]"
        },
        "description": "This is fene/expand bond"
    },
    "quartic": {
        "form": "K(r-Rc) ** 2 * (r-Rc-B1)*(r-Rc-B2) + U0 + 4*epsilon*((sigma/r) ** 12 - (sigma/r) ** 6) + epsilon",
        "parameters": ["K", "B1", "B2", "Rc", "U0"],
        "units": {
            "K": "[energy] [length] ** -4",
            "B1": "[length]",
            "B2": "[length]",
            "Rc": "[length]",
            "U0": "[energy]"
        },
        "description": "This is a quartic bond"
    },
    "hybrid": {  # Special case - allows for more than one bond type in a simulation
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "harmonic": {
        "form": "K*(r-r0) ** 2",
        "parameters": ["K", "r0"],
        "units": {
            "K": "[energy] [length] ** -2",
            "r0": "[length]"
        },
        "description": "This is a harmonic bond"
    },
    "table": {  # Special case - creation of interpolation tables.
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "class2": {
        "form": "K2 * (r-R0) ** 2 + K3 * (r-R0) ** 3 + K4 * (r-R0) ** 4",
        "parameters": ["R0","K2","K3","K4"],
        "units": {
            "R0": "[length]",
            "K2": "[energy] [length] ** -2",
            "K3": "[energy] [length] ** -3",
            "K4": "[energy] [length] ** -4"
        },
        "description": "This is a class2 bond"
    },
    "morse": {
        "form": "D * (1 - e ** (-alpha * (r-r0))) ** 2",
        "parameters": ["D","alpha","r0"],
        "units": {
            "D": "[energy]",
            "alpha": "[length] ** -1",
            "r0": "[length]"
        },
        "description": "This is a class2 bond"

    },
    "oxdna/fene": {
        "form": "-epsilon / 2 * ln (1 - ( (r - r0) / delta)**2 ) ",
        "parameters": ["epsilon","delta","r0"],
        "units": { 
            "epsilon": "[energy]",
            "delta": "[length]",
            "r0": "[length]"
        },
        "description": "This is a oxdna/fene bond"
    },
}

_angle_styles = {
    "none": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "class2": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "cosine/squared": {
        "form": "K*(cos(theta)-cos(theta0))**2",
        "parameters": ["K","theta0"],
        "units": {
            "K": "[energy]", 
            "theta0": "[degrees]"
         },
         "description": "This is a cosine/squared angle"
    },
    "zero": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "cosine": {
        "form": "K*(1+cos(theta))",
        "parameters": ["K"],
        "units": {
            "K": "[energy]"
         },
        "description": "This is a cosine potential"
    },
    "harmonic": {
        "form": "K*(theta-theta0)**2",
        "parameters": ["K","theta0"],
        "units": {
            "K": "[energy] [radian]**-2",
            "theta0": "[degree]"
        },
        "description": "This is a harmonic"
    },
    "hybrid": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },  # Special case - allows for more than one angle type in a simulation
    "cosine/delta": {
        "form": "K*(1+cos(theta-theta0))",
        "parameters": ["K","theta0"],
        "units": { 
            "K": "[energy]",
            "theta0": "[degree]"
        }, 
        "description": "This is a cosine/delta potential"
    },
    "table": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "charmm": {
        "form": "k*(theta-theta0)**2 + k_ub*(r-r_ub)**2",
        "parameters": ["k","theta0","k_ub","r_ub"],
        "units": {
            "k": "[energy] [radian]**-2",
            "theta0": "[degree]",
            "k_ub": "[energy] length**-2",
            "r_ub": "length"
        },
        "description": "This is a charmm angle potential"
    },
    "cosine/periodic": {
        "form": "C * (1-B*(-1)**n*cos(n*theta))",
        "parameters": ["C","B","n"],
        "units": {
            "C": "[energy]",
            "B": "[dimensionless]",  #1 or -1
            "n": "[dimensionless]"  # 1 2 3 4 5 or 6
        },
        "description": "This is a cosine/periodic potential"
    },
    "dipole": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "fourier": {
        "form": "K*(c0+c1*cos(theta)+c2*cos(2*theta))",
        "parameters": ["K","c0","c1","c2"],
        "units": {
            "K": "[energy]",
            "c0": "[dimensionless]",
            "c1": "[dimensionless]", 
            "c2": "[dimensionless]" 
        },
        "description": "This is a fourier potential"
    },
    "quartic": {
        "form": "K2*(theta-theta0)**2+K3*(theta-theta0)**3+K4*(theta-theta0)**4",
        "parameters": ["K2","K3","K4","theta0"],
        "units": { 
            "K2": "[energy] radian**-2",
            "K3": "[energy] radian**-3",
            "K4": "[energy] radian**-4",
            "theta0": "[degrees]" #Lammps converts this to radians
        },
        "description": "This is a quartic bond"
    },
    "sdk": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    }
}

_dihedral_styles = {
    "none": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "charmmfsw": {
        "form": "K*(1 + cos(n*phi-d))",
        "parameters": ["K","n","d"],
        "units": { 
            "K": "[energy]",
            "n": "[dimensionless]",  # must be type int, no units
            "d": "[degrees]",  # must be type int. Differs because units must be [degrees] regardless of units command ?
        },
        "description": "This is a charmm dihedral"
    },

    "multi/harmonic": {
        "form": "Sum(A_n * (cos(phi))**(n-1))", #n goes from 1 to 5
        "parameters": ["A_n"], 
        "units": {,
        "A_n": "[energy]",
         },
         "description": "This is a multi/harmonic dihedral"
    }
    "zero": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    # Class2 is complicated special case - see http://lammps.sandia.gov/doc/dihedral_class2.html
    "class2": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "cosine/shift/exp": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "fourier": {
        "form": "Sum(k_i * (1.0 + cos(n_i * phi - d_i)))",
        "parameters": ["k_i","n_i","d_i"],
        "units": { 
            "k_i": "[energy]",
            "n_i": "[dimensionless]",
            "d_i": "[degrees]"
        },
        "description": "This is a fourier dihedral"
    },
    "harmonic": {
        "form": "K*(1+d*cos(n*phi))",
        "parameters": ["K","n","d"],
        "units": { 
            "K": "[energy]",
            "n": "[dimensionless]",
            "d": "[dimensionless]"
        },
        "description": "This is a harmonic dihedral"
    },
    "helix": {
        "form": "A*(1-cos(phi)) + B*(1+cos(3*phi)) + C*(1+cos(phi+pi/4))",
        "parameters": ["A","B","C"],
        "units": { 
            "A": "[energy]",
            "B": "[energy]",
            "C": "[energy]"
        },
        "description": "This is a helix dihedral"
    },
    "hybrid": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "nharmonic": {
        "form": "Sum( A_n*(cos(phi))**(n-1))",
        "parameters": ["A_n","n"],
        "units": { 
            "A_n": "[energy]",
            "n": "[dimensionless]"
        },
        "description": "This is a nharmonic dihedral"
    },
    "opls": {
        "form": "0.5*K_1*(1+cos(phi))+0.5*K_2*(1-cos(2*phi))+0.5*K_3+(1+cos(3*phi))+0.5*K_4*(1-cos(4*phi))",
        "parameters": ["K_1","K_2","K_3", "K_4"],
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
        "parameters": ["K","phi0"],
        "units": { 
            "K": "[energy]",
            "phi0": "[degrees]"
        },
        "description": "This is a quadratic dihedral"
    },
    "spherical":
    {  #This type includes dihedrals phi and angles theta
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "table": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
}

_improper_styles = {
    "class2": {
        "form": "NYI",
        "parameters": ["","",""],
        "units": "NYI",
        "description": "NYI"
    },
    "cossq": {
        "form": "0.5*K*(cos(chi-chi0))**2",
        "parameters": ["K","chi0"],
        "units": { 
            "K": "[energy]",
            "chi0": "[degrees]"
        },
        "description": "This is a cossq improper"
    },
    "cvff": {
        "form": "K*(1+d*cos(n*chi))",
        "parameters": ["K","d","n"],
        "units": { 
            "K": "[energy]",
            "d": "[dimensionless]",
            "n": "[dimensionless]"
        },
        "description": "This is a cvff improper"
    },
    "distance": {
        "form": "K_2*r**2+K_4*r**4",
        "parameters": ["K_2","K_4"],
        "units": { 
            "K_2": "[energy] [length]**2",
            "K_4": "[energy] [length]**4",
        },
        "description": "This is a distance improper"
    },
    "fourier": {
        "form": "K*(C0+C1*cos(omega)+C2*cos(2*omega))",
        "parameters": ["K","C0","C1","C2"],
        "units": { 
            "K": "[energy]",
            "C0": "[dimensionless]",
            "C1": "[dimensionless]",
            "C2": "[dimensionless]",
        },
        "description": "This is a fourier improper"
    },
    "harmonic": {
        "form": "K*(chi-chi0)**2",
        "parameters": ["K","chi0"],
        "units": { 
            "K": "[energy]",
            "chi0": "[degrees]",
        },
        "description": "This is a harmonic improper"
    },
    "hybrid": {        
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "none": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "ring": {
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "umbrella": {  #Used in the dreiding force field
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
    "zero": { 
        "form": "NYI",
        "parameters": "NYI",
        "units": "NYI",
        "description": "NYI"
    },
}
