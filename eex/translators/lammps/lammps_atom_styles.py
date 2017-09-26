"""
Dictionary for LAMMPS atom_style - determines format of LAMMPS atom section

"""

_atom_style ={

    "angle"	: ["atom_ID", "molecule_ID", "atom_type" "x" "y" "z"],
    "atomic": ["atom_ID", "atom_type", "x", "y", "z"],
    "body":	["atom_ID", "atom_type", "bodyflag", "mass", "x", "y", "z"],
    "bond":	["atom_ID", "molecule_ID", "atom_type", "x", "y",  "z"],
    "charge": ["atom_id", "atom_type", "charge", "x", "y","z"],
    "dipole": ["atom_id", "atom_type", "charge", "x", "y", "z", "mux", "muy", "muz"],
    "dpd": ["atom_id", "atom_type", "theta", "x", "y","z"],
    "edpd":	["atom_id", "atom_type", "edpd_temp", "edpd_cv", "x", "y","z"],
    "mdpd":	["atom_id", "atom_type", "x", "y","z"],
    # "tdpd"	["atom_id", "atom_type", "x", "y", "z", "cc1", "cc2", ccNspecies], #NYI
    "electron":	["atom_id", "atom_type", "charge", "spin eradius", "x", "y","z"],
    "ellipsoid": ["atom_id", "atom_type", "ellipsoidflag density", "x", "y","z"],
    "full":	["atom_id", "molecule-ID", "atom_type", "charge", "x", "y","z"],
    "line": ["atom_id", "molecule-ID", "atom_type", "lineflag", "density", "x", "y","z"],
    "meso":	["atom_id", "atom_type", "rho", "e", "cv", "x", "y","z"],
    "molecular": ["atom_id", "molecule-ID" "atom_type", "x", "y","z"],
    "peri":	["atom_id", "atom_type", "volume density", "x", "y","z"],
    "smd":	["atom_id", "atom_type", "molecule volume", "mass", "kernel_radius", "contact_radius" "x", "y","z"],
    "sphere": ["atom_id", "atom_type", "diameter", "density", "x", "y","z"],
    "template": ["atom_id", "molecule_ID" "template_index", "template_atom", "atom_type", "x", "y","z"],
    "tri":	["atom_id", "molecule_ID", "atom_type", "triangleflag", "density", "x", "y","z"],
    "wavepacket":["atom_id", "atom_type", "charge", "spin", "eradius", "etag", "cs_re", "cs_im", "x", "y","z"],
    # "hybrid":["atom_id", "atom_type", "x", "y", "z", "sub-style1 sub-style2 ..."]


}

