"""
Contains metadata for additional simulation details 
"""

box_metadata = {
    "dimensions": {
        "a": "[length]",
        "b": "[length]",
        "c": "[length]",
        "alpha": "radian",
        "beta": "radian",
        "gamma": "radian",
        },

        "center": {
            "x": "[length]",
            "y": "[length]",
            "z": "[length]",
        },

        "boundary_conditions" : {
            "x": ["shrink-wrapped", "periodic", "fixed"],
            "y": ["shrink-wrapped", "periodic", "fixed"],
            "z": ["shrink-wrapped", "periodic", "fixed"]
        },

    "description": "This is information defines the simulation box",
}
# possible boundary conditions. Can be periodic (molecule moves through one side to the other),
# fixed - particle is lost when moves beyond box boundary, or "shrink-wrapped", meaning that that
# the box size adjusts


# These should be the keyword possibilities - function for computing based on keyword. May not go here
_box_center = ["origin", "center", "custom"]

_electrostatics = {
        "ewald": {
            "parameters": ["alpha", "accuracy", "kmax"],
            "units": {
                "alpha": "[length]**-1",
                "accuracy": "dimensionless",
                "kmaxx": "count",
                "kmaxy": "count",
                "kmaxz": "count",
            },
        "description": "Information for the Ewald method", 
        },
        "pme": {
            "parameters": ["g_ewald", "grid_size", "order", "accuracy"],
            "units": {
                "g_ewald": "[length]**-1",
                "grid_size_x": "dimensionless",
                "grid_size_y": "dimensionless",
                "grid_size_z": "dimensionless",
                "order": "dimensionless",
                "accuracy": "dimensionless",
            },
        "description": "Information for the PME method", 
        },
        "pppm": {
            "parameters": ["g_ewald", "grid_size", "order", "accuracy"],
            "units": {
                "g_ewald": "[length]**-1",
                "grid_size_x": "dimensionless",
                "grid_size_y": "dimensionless",
                "grid_size_z": "dimensionless",
                "order": "dimensionless",
                "accuracy": "dimensionless",
            },
        "description": "Information for the PPME method", 
        },
        "wolf": {
            "parameters": ["alpha", "accuracy"],
            "units": {
                "alpha": "[length] ** -1",
            },
        "description": "Information for the Wolf method", 
        },
        "dsf": {
            "parameters": ["alpha"],
            "units": {
                "alpha": "[length] ** -1",
            },
        "description": "Information for the DSF method", 
        },
        "reaction": {
            "parameters": ["A"],
            "units": {
                "A": "dimensionless", #This term depends on dielectric
            },
        "description": "Information for the reaction method", 
        },
        "cut": {
            "parameters": [],
            "units": {
            },
        "description": "Information for the cut method", 
        },
        "shift":{
            "parameters": [],
            "units": {
            },
        "description": "Information for the shift method", 
        },
        "debye": {
            "parameters": ["alpha"],
            "units": {
                "alpha": "[length]**-1", 
            },
        "description": "Information for the Debye method", 
        },
    "cutoff": "[length]",
}

_van_der_waals = {
    "cutoff": "[length]",
    "long": 
        {
            "standard":{
                "parameters": [],
                "units": {}
                },
            "switching": {
                "parameters":["switching_radius"],
                "units": {
                    "switching_radius": "[length]",
                    },
                },
            "shift": {
                "parameters": [],
                "units": {}
                },
            "cut": {
                "parameters": [],
                "units": {}
                },
    }
}

_mixing_rule = ["lorentz_berthelot", "geometric", "kong", "sixth_power", "custom"]

_neighbor = {
        "verlet": {
        "skin": "[length]",
        "frequency": "count",
    },
}

#_special_bonds = {
#}

exclusions = {
    "coul":{
        "scale12": "dimensionless",
        "scale13": "dimensionless",
        "scale14": "dimensionless",
    }, 
    "lj":{
        "scale12": "dimensionless",
        "scale13": "dimensionless",
        "scale14": "dimensionless",
    }
}

_torsion_convention = ["180_is_trans", "0_is_trans"]

#_groups = {
#}
