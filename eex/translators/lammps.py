"""
LAMMPS EEX I/O
"""

import pandas as pd
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
    "Masses": "atom types",
    "Pair Coeffs": "atom types",
    "Bond Coeffs": "bond types",
    "Angle Coeffs": "angle types",
    "Dihedral Coeffs": "dihedral types",
    "Improper Coeffs": "improper types"
}

# Data labels to look for their size
_data_labels = {
    "Atoms": "atoms",
    "Bonds": "bonds",
    "Angles": "angles",
    "Dihedrals": "dihedrals",
    "Impropers": "impropers"
}

_full_labels = _coeff_labels.copy()
_full_labels.update(_data_labels)

_full_labels_list = list(_full_labels)


def read_lammps_file(dl, filename):

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

    breakline = None
    first_data_type = None

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
            breakline = num
            first_data_type = utility.line_fuzzy_list(line, _full_labels_list)[1]
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
    if breakline is None:
        raise IOError("LAMMPS Read: Did not find data start in %d header lines." % max_rows)

    if sum((v != None) for k, v in dim_dict.items()) != 6:
        raise IOError("LAMMPS Read: Did not find dimension data in %d header lines." % max_rows)

    if ("atoms" not in list(sizes_dict)) or ("atom types" not in list(sizes_dict)):
        raise IOError("LAMMPS Read: Did not find size data on 'atoms' or 'atom types' in %d header lines." % max_rows)

    print(sizes_dict)
    print(dim_dict)
    ### Iterate over the primary data portion of the object
    # print(breakline)
    # print(first_data_type)
    # while True:

    #     if first_data_type not in list(sizes_dict):
    #         raise KeyError("")

    # raise Exception("")
    data = {}
    data["sizes"] = sizes_dict
    data["dimensions"] = dim_dict

    return data
