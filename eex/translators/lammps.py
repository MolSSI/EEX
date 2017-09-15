"""
LAMMPS EEX I/O
"""

import pandas as pd
import math

from .. import datalayer
from .. import utility

import logging
logger = logging.getLogger(__name__)

# Possible size keys to look for in the header
_size_keys = [
    "atoms", "atom types", "bonds", "bond types", "angles", "angle types", "dihedrals", "dihedral types", "impropers",
    "improper types"
]

# Coeff labels to look for and their size
_coeff_labels = {
    "Pair Coeffs": ("atom types", "NYI"),
    "Bond Coeffs": ("bond types", "NYI"),
    "Angle Coeffs": ("angle types", "NYI"),
    "Dihedral Coeffs": ("dihedral types", "NYI"),
    "Improper Coeffs": ("improper types", "NYI"),
}

# Data labels to look for their size
_data_labels = {
    "Masses": ("atom types", "NYI"),
    "Atoms": ("atoms", "add_atoms", ["atom_index", "molecule_index", "atom_type", "charge", "X", "Y", "Z"]),
    "Bonds": ("bonds", "add_bonds", ["bond_index", "bond_type", "atom1_index", "atom2_index"]),
    "Angles": ("angles", "add_angles", ["bond_index", "angle_type", "atom1_index", "atom2_index", "atom3_index"]),
    "Dihedrals": ("dihedrals", "NYI"),
    "Impropers": ("impropers", "NYI"),
}

# Units for data labels
# Not sure how we'll implement yet - just getting information in file
_data_labels_units = {
    "Masses": "mass",
    "Atoms": ["N/A", "N/A", "atom_type", "charge", "distance", "distance", "distance"]),

}


# Dictionaries for coeff format and units - http://lammps.sandia.gov/doc/Section_commands.html
# NYI
_bond_styles = {
    "none" : {},
    "fene": {
        "bond_coeff" : {}
    },
    "nonlinear": {},
    "zero": {},
    "quartic": {},
    "hybrid": {},
    "harmonic": {},
    "table": {},
    "class2": {},
    "morse":{},
}

_angle_styles = {
    "none": {},
    "class2": {},
    "cosine/squared": {},
    "zero": {},
    "cosine": {},
    "harmonic": {},
    "hybrid": {},
    "cosine/delta": {},
    "charmm": {},
    "cosine/periodic": {},
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

_full_labels = _coeff_labels.copy()
_full_labels.update(_data_labels)

_full_labels_list = list(_full_labels)


def _get_size(label, sizes_dict):
    if label not in list(_full_labels):
        raise KeyError("LAMMPS: Label '%s' not recognized" % label)

    skey = _full_labels[label][0]
    return sizes_dict[skey]


def _get_dl_function(label):
    if label not in list(_full_labels):
        raise KeyError("LAMMPS: Label '%s' not recognized" % label)

    return _full_labels[label][1]

def _get_df_columns(label):
    if label not in list(_full_labels):
        raise KeyError("LAMMPS: Label '%s' not recognized" % label)

    if len(_full_labels[label]) >= 3:
        return _full_labels[label][2]
    else:
        return False


def read_lammps_file(dl, filename, blocksize=110):

    ### First we need to figure out system dimensions
    max_rows = 100  # How many lines do we attempt to search?
    with open(filename, "r") as infile:
        header_data = [next(infile).strip() for x in range(max_rows)]

    dim_dict = {
        "xlo": None,
        "xhi": None,
        "ylo": None,
        "yhi": None,
        "zlo": None,
        "zhi": None,
    }

    sizes_dict = {}

    startline = None
    current_data_category = None

    header = header_data[0]
    for num, line in enumerate(header_data[1:]):

        # Skip blanklines
        if line == "":
            continue

        # Skip comment line
        elif line[0] == "#":
            continue

        # We are
        elif utility.line_fuzzy_list(line, _full_labels_list)[0]:
            startline = num + 3  # Skips first row and two blank lines
            current_data_category = utility.line_fuzzy_list(line, _full_labels_list)[1]
            break

        # Figure out the dims
        elif ("lo" in line) and ("hi" in line):
            dline = line.split()
            if dline[-1] == "xhi":
                dim_dict["xlo"] = float(dline[0])
                dim_dict["xhi"] = float(dline[1])

            elif dline[-1] == "yhi":
                dim_dict["ylo"] = float(dline[0])
                dim_dict["yhi"] = float(dline[1])
            elif dline[-1] == "zhi":
                dim_dict["zlo"] = float(dline[0])
                dim_dict["zhi"] = float(dline[1])
            else:
                raise KeyError(
                    "LAMMPS Read: The following line looks like a dimension line, but does not match:\n%s" % line)

        # Are we a size line?
        elif utility.line_fuzzy_list(line, _size_keys)[0]:
            dline = line.split()
            size = int(dline[0])
            size_name = " ".join(dline[1:])

            if size_name in list(sizes_dict):
                raise KeyError("LAMMPS Read: KeyError size key %s already found." % size_name)
            elif size_name not in _size_keys:
                raise KeyError("LAMMPS Read: KeyError size key %s not recognized." % size_name)
            else:
                sizes_dict[size_name] = size

        else:
            raise IOError("LAMMPS Read: Line not understood!\n%s" % line)

    # Make sure we have what we need
    if startline is None:
        raise IOError("LAMMPS Read: Did not find data start in %d header lines." % max_rows)

    if sum((v != None) for k, v in dim_dict.items()) != 6:
        raise IOError("LAMMPS Read: Did not find dimension data in %d header lines." % max_rows)

    if ("atoms" not in list(sizes_dict)) or ("atom types" not in list(sizes_dict)):
        raise IOError("LAMMPS Read: Did not find size data on 'atoms' or 'atom types' in %d header lines." % max_rows)

    ### Iterate over the primary data portion of the object

    reader = pd.read_table(
        filename,
        header=None,
        iterator=True,
        names=range(10),
        engine="c",
        comment="#",
        delim_whitespace=True,
        skiprows=startline)

    while True:

        # Figure out the size of the chunk to read
        size = _get_size(current_data_category, sizes_dict)
        dl_func = _get_dl_function(current_data_category)
        df_cols = _get_df_columns(current_data_category)

        # Read in the data, in chunks
        remaining = size
        num_blocks = int(math.ceil(size / float(blocksize)))
        for block in range(num_blocks):

            # Figure out the size of the read
            read_size = blocksize
            if remaining < blocksize:
                read_size = remaining

            # Read and update DL
            data = reader.get_chunk(read_size).dropna(axis=1, how="all")
            if dl_func != "NYI":
                # print(data)
                data.columns = df_cols
            dl.call_by_string(dl_func, data)

            # Update remaining
            remaining -= blocksize

        # Figure out the next category to read
        try:
            tmp = reader.get_chunk(1).dropna(axis=1, how="any")
        except StopIteration:
            break

        current_data_category = " ".join(str(x) for x in list(tmp.iloc[0]))

    # raise Exception("")
    data = {}
    data["sizes"] = sizes_dict
    data["dimensions"] = dim_dict

    return data
