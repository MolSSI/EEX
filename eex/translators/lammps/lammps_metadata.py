"""
A helper file to produce LAMMPS metadata based on the type of computation
"""
import copy

from . import lammps_ff
import eex

# Keywords for different outputs in log file
output_keywords = {
            'Step':{
                'units': 'count', 
                'canonical': ['step'],
            },
            'Elapsed':{
                'units': '[time]',
                'canonical': ['elapsed_time'],
            },
            'Elaplong':{
                'units': '[time]',
                'canonical': ['elapsed_long'],
            },
            'Dt':{
                'units': '[time]',
                'canonical': ['timestep'],
            },
            'Time':{
                'units': '[time]',
                'canonical': ['simulation_time'],
            },
            'CPU':{
                'units': '[time]',
                'canonical': ['cpu_time'],
            },
            'T/CPU':{
                'units': '[time]',
                'canonical': ['time_per_cpu'],
            },
            'S/CPU':{
                'units': 'count', 
                'canonical': ['timestep_per_cpu'],
            },
            'CPULeft':{
                'units': '[time]',
                'canonical': ['cpu_left'],
            },
            'Part':{
                'units': 'count',
                'canonical': ['part'],
            },
            'TimeoutLeft':{
                'units': '[time]',
                'canonical': ['time_remain'],
            },
            'Atoms':{
                'units': 'count',
                'canonical': ['num_atoms'],
            },
            'Temp':{
                'units': '[temperature]',
                'canonical': ['temperature'],
            },
            'Press':{
                'units': '[pressure]',
                'canonical': ['pressure'],
            },
            'PotEng':{
                'units': '[energy]',
                'canonical': ['potential'],
            },
            'KinEng':{
                'units': '[energy]',
                'canonical': ['kinetic'],
            },
            'TotEng':{
                'units': '[energy]',
                'canonical': ['total'],
            },
            'Enthalpy':{
                'units': '[energy]',
                'canonical': ['enthalpy'],
            },
            'E_vdwl':{
                'units': '[energy]',
                'canonical': ['vdw_total'],
            },
            'E_coul':{
                'units': '[energy]',
                'canonical': ['coulomb_total'],
            },
            'E_pair':{
                'units': '[energy]',
                'canonical': ['pairwise'],
            },
            'E_bond':{
                'units': '[energy]',
                'canonical': ['bond'],
            },
            'E_angle':{
                'units': '[energy]',
                'canonical': ['angle'],
            },
            'E_dihed':{
                'units': '[energy]',
                'canonical': ['proper'],
            },
            'E_impro':{
                'units': '[energy]', 
                'canonical': ['improper'],
            },
            'E_mol':{
                'units': '[energy]',
                'canonical': ['molecular'],
            },
            'E_long':{
                'units': '[energy]',
                'canonical': ['kspace'],
            },
            'E_tail':{
                'units': '[energy]',
                'canonical': ['tail'],
            },
            'Volume':{
                'units': '[length] ** 3',
                'canonical': ['volume'],
            },
            'Density':{
                'units': '[mass] * [length] ** -3',
                'canonical': ['density'],
            },
            'Lx':{
                'units': '[length]',
                'canonical': ['lx'],
            },
            'Ly':{
                'units': '[length]',
                'canonical': ['ly'],
            },
            'Lz':{
                'units': '[length]',
                'canonical': ['lz'],
            },
            'Xlo':{
                'units': '[length]',
                'canonical': ['xlo'],
            },
            'Xhi':{
                'units': '[length]',
                'canonical': ['xhi'],
            },
            'Ylo':{
                'units': '[length]',
                'canonical': ['ylo'],
            },
            'Yhi':{
                'units': '[length]',
                'canonical': ['yhi'],
            },
            'Zlo':{
                'units': '[length]',
                'canonical': ['zlo'],
            },
            'Zhi':{
                'units': '[length]',
                'canonical': ['zhi'],
            },
            'Xy':{
                'units': '[length]',
                'canonical': ['xy'],
            },
            'Xz':{
                'units': '[length]',
                'canonical': ['xz'],
            },
            'Yz':{
                'units': '[length]',
                'canonical': ['yz'],
            },
            'Xlat':{
                'units': '[length]',
                'canonical': ['xlat'],
            },
            'Ylat':{
                'units': '[length]',
                'canonical': ['ylat'],
            },
            'Zlat':{
                'units': '[length]',
                'canonical': ['zlat'],
            },
            'Bonds':{
                'units': 'count',
                'canonical': ['num_bonds'],
            },
            'Angles':{
                'units': 'count',
                'canonical': ['num_angles'],
            },
            'Dihedrals':{
                'units': 'count',
                'canonical': ['num_proper'],
            },
            'Impros':{
                'units': 'count',
                'canonical': ['num_improper'],
            },
            'Pxx':{
                'units': '[pressure]',
                'canonical': ['num_pxx'],
            },
            'Pyy':{
                'units': '[pressure]',
                'canonical': ['num_pyy'],
            },
            'Pzz':{
                'units': '[pressure]',
                'canonical': ['num_pzz'],
            },
            'Pxy':{
                'units': '[pressure]',
                'canonical': ['num_pxy'],
            },
            'Pxz':{
                'units': '[pressure]',
                'canonical': ['num_pxz'],
            },
            'Pyz':{
                'units': '[pressure]',
                'canonical': ['num_pyz'],
            },
            'Fmax':{
                'units': '[force]',
                'canonical': ['fmax'],
            },
            'Fnorm':{
                'units': '[force]',
                'canonical': ['fnorm'],
            },
            'Nbuild':{
                'units': 'count',
                'canonical': ['nbuild'],
            },
            'Ndanger':{
                'units': 'count',
                'canonical': ['ndanger'],
            },
            'Cella':{
                'units': '[length]',
                'canonical': ['cella'],
            },
            'Cellb':{
                'units': '[length]',
                'canonical': ['cellb'],
            },
            'Cellc':{
                'units': '[length]',
                'canonical': ['cellc'],
            },
            'CellAlpha':{
                'units': 'radian',
                'canonical': ['cellalpha'],
            },
            'CellBeta':{
                'units': 'radian',
                'canonical': ['cellbeta'],
            },
            'CellGamma':{
                'units': 'radian',
                'canonical': ['cellgamma'],
            },
        }


# Possible types of the keyword 'variable' found in the input file
variable_types = ['delete', 'index', 'loop', 'world', 'universe', 'uloop', 'string', 'format', 'getenv', 'file', 'atomfile', 'python', 'internal', 'equal', 'vector', 'atom']

mixing_rules = ["geometric", "arithmetic", "sixthpower"]

exclusions = {
    "default": {
        "coul":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.0,
        },
        "lj":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.0,
        },
    },
    "amber": {
        "coul":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.8333,
        },
        "lj":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.5,
        },
    },
    "charmm": {
        "coul":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.0,
        },
        "lj":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 0.0,
        },
    },
    "dreiding": {
        "coul":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 1.0,
        },
        "lj":{
            "scale12": 0.0,
            "scale13": 0.0, 
            "scale14": 1.0,
        },
    },
    "fene": {
        "coul":{
            "scale12": 0.0,
            "scale13": 1.0, 
            "scale14": 1.0,
        },
        "lj":{
            "scale12": 0.0,
            "scale13": 1.0, 
            "scale14": 1.0,
        },
    },
    "lj/coul": {
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
    },
    "lj": {
        "lj":{
            "scale12": "dimensionless",
            "scale13": "dimensionless", 
            "scale14": "dimensionless",
        },
    },
    "coul": {
        "coul":{
            "scale12": "dimensionless",
            "scale13": "dimensionless", 
            "scale14": "dimensionless",
        },
    },
    "additional_keywords": ["angle", "dihedral"]
}


# Possible size keys to look for in the header
size_keys = [
    "atoms", "atom types", "bonds", "bond types", "angles", "angle types", "dihedrals", "dihedral types", "impropers",
    "improper types"
]

# Units based LAMMPS unit styles (http://lammps.sandia.gov/doc/units.html)
units_style = {
    "default": "real",
    "lj": {},
    "real": {
        "[mass]": "(gram mol ** -1)",
        "[length]": "angstrom",
        "[time]": "femtosecond",
        "[energy]": "(kcal mol ** -1)",
        "[velocity]": "(angstrom femtosecond ** -1)",
        "[force]": "(kcal * mol ** -1 angstrom ** -1)",
        "[torque]": "(kcal * mol ** -1)",
        "[temperature]": "kelvin",
        "[pressure]": "atmosphere",
        "[dynamic viscosity]": "poise",
        "[charge]": "e",
        "[dipole]": "e * angstrom",
        "[electric field]": "(volt angstrom ** -1)",
        "[density]": "(gram cm ** -dim)"
    },
    "metal": {
        "[mass]": "(gram mol ** -1)",
        "[length]": "(angstrom)",
        "[time]": "(picosecond)",
        "[energy]": "(eV)",
        "[velocity]": "(angstrom picosecond ** -1)",
        "[force]": "(eV angstrom ** -1)",
        "[torque]": "(eV)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(bar)",
        "[dynamic viscosity]": "(poise)",
        "[charge]": "(e)",
        "[dipole]": "(e * angstrom)",
        "[electric field]": "(volt angstrom ** -1)",
        "[density]": "(gram cm ** -dim)",
    },
    "si": {
        "[mass]": "(kilogram)",
        "[length]": "(meter)",
        "[time]": "(second)",
        "[energy]": "(joule)",
        "[velocity]": "(meter second ** -1)",
        "[force]": "(newtons)",
        "[torque]": "(newton meter)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(pascal)",
        "[dynamic viscosity]": "(pascal second)",
        "[charge]": "(coulomb)",
        "[dipole]": "(coulomb meter)",
        "[electric field]": "(volt meter ** -1)",
        "[density]": "(kilogram meter ** -dim)",
    },
    "cgs": {
        "[mass]": "(gram)",
        "[length]": "(centimeters)",
        "[time]": "(second)",
        "[energy]": "(ergs)",
        "[velocity]": "(centimeters second ** -1)",
        "[force]": "(dyne)",
        "[torque]": "(dyne centimeters)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(dyne cm ** -2)",
        "[dynamic viscosity]": "(poise)",
        "[charge]": "(statcoulombs)",
        "[dipole]": "(statcoulombs centimeter)",
        "[electric field]": "(statvolt cm ** -1)",
        "[density]": "(gram cm ** -dim)",
    },
    "electron": {
        "[mass]": "(amu)",
        "[length]": "(Bohr)",
        "[time]": "(femtosecond)",
        "[energy]": "(hartree)",
        "[velocity]": "(Bohr (atomic time unit) ** -1)",
        "[force]": "(hartree Bohr ** -1)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(pascal)",
        "[charge]": "(e)",
        "[dipole] moment": "(debye)",
        "[electric field]": "(volt cm ** -1)",
    },
    "micro": {
        "[mass]": "(picogram)",
        "[length]": "(micrometer)",
        "[time]": "(microsecond)",
        "[energy]": "(picogram micrometer ** 2 microsecond ** -2)",
        "[velocity]": "(micrometers microsecond ** -1)",
        "[force]": "(picogram micrometer microsecond ** -2)",
        "[torque]": "(picogram micrometer ** 2 microsecond ** -2)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(picogram micrometer ** -1 microsecond ** -2)",
        "[dynamic viscosity]": "(picogram micrometer ** -1 microsecond ** -1)",
        "[charge]": "(picocoulomb)",
        "[dipole]": "(picocoulomb micrometer)",
        "[electric field]": "(volt micrometer ** -1)",
        "[density]": "(picogram micrometer ** -dim)"
    },
    "nano": {
        "[mass]": "(attogram)",
        "[length]": "(nanometer)",
        "[time]": "(nanosecond)",
        "[energy]": "(attogram nanometer ** 2 nanosecond ** -2)",
        "[velocity]": "(nanometers nanosecond ** -1)",
        "[force]": "(attogram nanometer nanosecond ** -2)",
        "[torque]": "(attogram nanometer ** 2 nanosecond ** -2)",
        "[temperature]": "(kelvin)",
        "[pressure]": "(attogram nanometer ** -1 nanosecond ** -2)",
        "[dynamic viscosity]": "(attogram nanometer ** -1 nanosecond ** -1)",
        "[charge]": "(e)",
        "[dipole]": "(e nanometer)",
        "[electric field]": "(volt nanometer ** -1)",
    },
}

# Possible atom styles
atom_style = {
    "default": "full",
    "angle": ["atom_index", "molecule_index", "atom_type", "X", "Y", "Z"],
    "atomic": ["atom_index", "atom_type", "X", "Y", "Z"],
    "body": ["atom_index", "atom_type", "bodyflag", "mass", "X", "Y", "Z"],
    "bond": ["atom_index", "molecule_index", "atom_type", "X", "Y", "Z"],
    "charge": ["atom_index", "atom_type", "charge", "X", "Y", "Z"],
    "dipole": ["atom_index", "atom_type", "charge", "X", "Y", "Z", "mux", "muy", "muz"],
    "dpd": ["atom_index", "atom_type", "theta", "X", "Y", "Z"],
    "edpd": ["atom_index", "atom_type", "edpd_temp", "edpd_cv", "X", "Y", "Z"],
    "mdpd": ["atom_index", "atom_type", "X", "Y", "Z"],
    # "tdpd"    ["atom_index", "atom_type", "X", "Y", "Z", "cc1", "cc2", ccNspecies], #NYI
    "electron": ["atom_index", "atom_type", "charge", "spin eradius", "X", "Y", "Z"],
    "ellipsoid": ["atom_index", "atom_type", "ellipsoidflag density", "X", "Y", "Z"],
    "full": ["atom_index", "molecule_index", "atom_type", "charge", "X", "Y", "Z"],
    "line": ["atom_index", "molecule_index", "atom_type", "lineflag", "density", "X", "Y", "Z"],
    "meso": ["atom_index", "atom_type", "rho", "e", "cv", "X", "Y", "Z"],
    "molecular": ["atom_index", "molecule_index"
                  "atom_type", "X", "Y", "Z"],
    "peri": ["atom_index", "atom_type", "volume density", "X", "Y", "Z"],
    "smd": ["atom_index", "atom_type", "molecule volume", "mass", "kernel_radius", "contact_radius"
            "X", "Y", "Z"],
    "sphere": ["atom_index", "atom_type", "diameter", "density", "X", "Y", "Z"],
    "template": ["atom_index", "molecule_index"
                 "template_index", "template_atom", "atom_type", "X", "Y", "Z"],
    "tri": ["atom_index", "molecule_index", "atom_type", "triangleflag", "density", "X", "Y", "Z"],
    "wavepacket": ["atom_index", "atom_type", "charge", "spin", "eradius", "etag", "cs_re", "cs_im", "X", "Y", "Z"],
    # "hybrid":["atom_index", "atom_type", "X", "Y", "Z", "sub-style1 sub-style2 ..."]
}

# Units for data labels
_atom_utypes = {"mass": "[mass]", "charge": "[charge]", "xyz": "[length]"}

_operation_table = {

    # Add atom data
    "Atoms": {
        "size": "atoms",
        "dl_func": "add_atoms",
        "df_cols": "atom_style",
        "kwargs": {
            "utype": None,
            "by_value": True,
        },
        "call_type": "single",
    },

    # Add any other atom parameters
    "Masses": {
        "size": "atom types",
        "dl_func": "add_atom_parameters",
        "call_type": "add_atom_parameters",
        "atom_property": "mass",
        "kwargs": {
            "utype": None
        },
    },

    # Add term data by index
    "Bonds": {
        "size": "bonds",
        "dl_func": "add_bonds",
        "df_cols": ["bond_index", "term_index", "atom1", "atom2"],
        "call_type": "single",
    },
    "Angles": {
        "size": "angles",
        "dl_func": "add_angles",
        "df_cols": ["angle_index", "term_index", "atom1", "atom2", "atom3"],
        "call_type": "single",
    },
    "Dihedrals": {
        "size": "dihedrals",
        "dl_func": "add_dihedrals",
        "df_cols": ["dihedral_index", "term_index", "atom1", "atom2", "atom3", "atom4"],
        "call_type": "single",
    },
    "Impropers": {
        "size": "impropers",
        "dl_func": "NYI",
        "call_type": "single",
    },

    # Add long range pair data
    "Pair Coeffs": {
        "size": "atom types",
        "n_uids": 1,
        "dl_func": "add_nb_parameter",
        "call_type": "nb_parameter",
        "kwargs": {
            "atom_type": None,
            "nb_parameters": None,
            "nb_name": "LJ",
            "nb_model": "epsilon/sigma",
            "utype": None
            }
    },

    "PairIJ Coeffs": {
        "size": "pair types",
        "n_uids": 2,
        "dl_func": "add_nb_parameter",
        "call_type": "nb_parameter",
        "kwargs": {
            "atom_type": None,
            "atom_type2": None,
            "nb_parameters": None,
            "nb_name": "LJ",
            "nb_model": "epsilon/sigma",
            "utype": None
            }
    },

    # Add term parameters
    "Bond Coeffs": {
        "size": "bond types",
        "dl_func": "add_parameters",
        "call_type": "parameter",
        "args": {
            "order": 2,
            "style_keyword": "bond_style"
        },
    },
    "Angle Coeffs": {
        "size": "angle types",
        "dl_func": "add_parameters",
        "call_type": "parameter",
        "args": {
            "order": 3,
            "style_keyword": "angle_style"
        },
    },
    "Dihedral Coeffs": {
        "size": "dihedral types",
        "dl_func": "add_parameters",
        "call_type": "parameter",
        "args": {
            "order": 4,
            "style_keyword": "dihedral_style"
        },
    },
    "Improper Coeffs": {
        "size": "improper types",
        "dl_func": "NYI",
        "call_type": "parameter",
        # "dl_func": "add_parameters"
    },
}

# Define a few temporaries

def get_context(utype, context):
    return units_style[utype][context]


def build_atom_units(utype):
    ret = {}

    unit_values = units_style[utype]
    for k, v in _atom_utypes.items():
        ret[k] = unit_values[v]
    return ret


def build_valid_category_list():
    return list(_operation_table)


def build_operation_table(extra_simulation_data, size_dict):

    # Get the first operations table
    ret = copy.deepcopy(_operation_table)

    aunits = build_atom_units(extra_simulation_data["units"])
    for k, v in ret.items():

        # Supply the Atom unit types
        if v["dl_func"] in ["add_atom_parameters", "add_atoms"]:
            v["kwargs"]["utype"] = aunits
 
        # Supply bond types, angle types, dihedral types
        if ("args" in v) and ("style_keyword" in v["args"]) and (v["args"]["style_keyword"] in ["bond_style", "angle_style", "dihedral_style"]):
            if (v["args"]["style_keyword"] in extra_simulation_data.keys()):
                v["args"]["style_keyword"] = extra_simulation_data[v["args"]["style_keyword"]]

        # Supply atom style
        if "df_cols" in v and "atom_style" in v["df_cols"]:
            v["df_cols"] = atom_style[extra_simulation_data["atom_style"]]


        # Build sizes dict
        try:
            v["size"] = size_dict[v["size"]]
        except KeyError:
            v["size"] = 0

        # Blank kwargs if needed
        if "kwargs" not in v:
            v["kwargs"] = {}

    return ret


def build_term_table(utype):

    ustyle = units_style[utype]
    ret = copy.deepcopy(lammps_ff.term_data)

    # Loop over order
    for ok, ov in ret.items():
        # Loop over terms
        for k, v in ov.items():
            # Loop over parameters
            utype = {}
            for pk, pv in v["units"].items():
                # print(ok, k, pk)
                utype[pk] = eex.units.convert_contexts(pv, ustyle)
            v["utype"] = utype
    return ret


def build_nb_table(utype):

    ustyle = units_style[utype]
    ret = copy.deepcopy(lammps_ff.nb_styles)

    for k, v in ret.items():
        utype = {}
        for pk, pv in v["units"].items():
            utype[pk] = eex.units.convert_contexts(pv, ustyle)
        v["utype"] = utype
    return ret


def build_prop_table(utype):

    ustyle = units_style[utype]
    ret = copy.deepcopy(output_keywords)

    for k, v in ret.items():
        utype = eex.units.convert_contexts(v['units'], ustyle)
        v['utype'] = utype
    return ret

# if __name__ == "__main__":
#     from eex.units import ureg

#     # Test unit md
#     for k, v in units_style.items():
#         for k, v in units_style[k].items():
#             try:
#                 u = ureg.parse_expression(v)
#                 # print(u)
#             except:
#                 print("Pint could not parse %s" % v)

#     build_atom_units("real")
#     build_operation_table("real")
#     build_term_table("real")
