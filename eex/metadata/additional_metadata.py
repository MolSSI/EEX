"""
Contains metadata for additional simulation details 
"""

_box_info = {
    "dimensions": {
        "lx": "[length]",
        "ly": "[length]",
        "lz": "[length]",
        "alpha": "[length]",
        "beta": "[length]",
        "gamma": "[length]",
        },
    "description": "This is information defines the simulation box",
}

_boundary = {
    "x": ["p", "f"], 
    "y": ["p", "f"], 
    "z": ["p", "f"], 
}

_box_origin = ["origin", "center", "custom"]

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
            "shift"{
                "parameters": [],
                "units": {}
                },
            "cut"{
                "parameters": [],
                "units": {}
                },
}

_mixing_rule = ["lorentz_berthelot", "geometric", "kong", "sixth_power", "custom"]

_neighbor = {
    "verlet" = {
        {
        "skin": "[length]",
        "frequency": "count",
        },
}

#_special_bonds = {
#}

_exclusions = {
    "electrostatic":{
        "scale12": "dimensionless",
        "scale13": "dimensionless",
        "scale14": "dimensionless",
    }
    "van_der_waals":{
        "scale12": "dimensionless",
        "scale13": "dimensionless",
        "scale14": "dimensionless",
    }
}

_torsion_convention = ["180_is_trans", "0_is_trans"]

#_groups = {
#}
