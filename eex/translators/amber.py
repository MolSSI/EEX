"""
LAMMPS EEX I/O
"""

import pandas as pd
import math
import re

from .. import datalayer
from .. import utility

import logging
logger = logging.getLogger(__name__)

# Possible size keys to look for in the header
_size_keys = [
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH", "MPHIA", "NHPARM", "NPARM", "NNB", "NRES",
    "NBONA", "NTHETA", "NPHIA", "NUMBND", "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT", "NBPER", "NGPER", "NDPER",
    "MBPER", "MGPER", "MDPER", "IFBOX", "NMXRS", "IFCAP", "NUMEXTRA"
]

_data_labels = {
    "ATOM_NAME": "NATOM",
    "CHARGE": "NATOM",
    "ATOMIC_NUMBER": "NATOM",
    "MASS": "NATOM",
    "ATOM_TYPE_INDEX": "NATOM",
    "NUMBER_EXCLUDED_ATOMS": "NATOM",
    "NONBONDED_PARM_INDEX": "NATOM ** 2",
    "RESIDUE_LABEL": "NRES",
    "RESIDUE_POINTER": "NRES",
    "BOND_FORCE_CONSTANT": "NUMND",
    "BOND_EQUIL_VALUE": "NUMND",
    "ANGLE_FORCE_CONSTANT": "NUMANG",
    "ANGLE_EQUIL_VALUE": "NUMND",
    "DIHEDRAL_FORCE_CONSTANT": "NPTRA",
    "DIHEDRAL_PERIODICITY": "NPTRA",
    "DIHEDRAL_PHASE": "NPTRA",
    "SCEE_SCALE_FACTOR": "NPTRA",
    "SCNB_SCALE_FACTOR": "NPTRA",
    "SOLTY": "NATYP",
    "LENNARD_JONES_ACOEF": "(NTYPES * (NTYPES + 1)) / 2",
    "LENNARD_JONES_BCOEF": "(NTYPES * (NTYPES + 1)) / 2",
    "BONDS_INC_HYDROGEN": "3 * NBONH",
    "BONDS_WITHOUT_HYDROGEN": "3 * NBONA",
    "ANGLES_INC_HYDROGEN": "4 * NTHETH",
    "ANGLES_WITHOUT_HYDROGEN": "4 * NTHETA",
    "DIHEDRALS_INC_HYDROGEN": "5 * NPHIH",
    "DIHEDRALS_WITHOUT_HYDROGEN": "5 * NPHIA",
    "EXCLUDED_ATOMS_LIST": "NNB",
    "HBOND_ACOEF": "NPHB",
    "HBOND_BCOEF": "NPHB",
    "HBCUT": "NPHB",
    "AMBER_ATOM_TYPE": "NATOM",
    "TREE_CHAIN_CLASSIFICATION": "NATOM",
    "JOIN_ARRAY": "NATOM",
    "IROTAT": "NATOM",
    "SOLVENT_POINTERS": 3,
    "ATOMS_PER_MOLECULE": "NPSM",
    "BOX_DIMENSIONS": 4,
    "CAP_INFO": "1 if IFCAP else 0",
    "CAP_INFO2": "4 if IFCAP else 0",
    "RADIUS_SET": 1,
    "RADII": "NATOM",
    "IPOL": 1,
    "POLARIZABILITY": "NATOM if IPOL else 0"
}


def _parse_format(string):
    """
    Parses an AMBER style format string.

    Example:

    >>> string = "%FORMAT(10I8)"
    >>> _parse_format(string)
    [10, int, 8]
    """
    if "FORMAT" not in string:
        raise ValueError("AMBER: Did not undstand format line '%s'." % string)

    string = string.replace("%FORMAT(", "").replace(")", "").strip()
    ret = [x for x in re.split('(\d+)', string) if x]
    if len(ret) != 3:
        raise ValueError("AMBER: Did not undstand format line '%s'." % string)

    ret[0] = int(ret[0])
    ret[2] = int(ret[2])

    if ret[1] == "I":
        ret[1] = int
    elif ret[1] == "F":
        ret[1] = float
    else:
        raise ValueError("AMBER: Type symbol '%s' not understood." % ret[1])

    return ret


def read_amber_file(dl, filename, blocksize=110):

    ### First we need to figure out system dimensions
    max_rows = 100  # How many lines do we attempt to search?
    with open(filename, "r") as infile:
        header_data = [next(infile).strip() for x in range(max_rows)]

    sizes = {}
    data = {}
    found_sizes = False
    found_version = False
    for num, line in enumerate(header_data):

        # First try to figure out the version
        if "VERSION" in line:
            sline = line.strip().split()
            if "VERSION_STAMP" == sline[1]:
                data["VERSION"] = sline[3]
                if data["VERSION"] != "V0001.000":
                    raise ValueError("AMBER Read: Did not recognize version '%s'." % data["VERSION"])
                found_version = True
            else:
                raise ValueError("AMBER Read: Could not understand version line.")

        # Get size_information
        elif "FLAG POINTERS" in line:
            ncols, dtype, width = _parse_format(header_data[num + 1])

            parsed_sizes = []
            for shift in range(4):
                dline = header_data[num + 2 + shift]
                dline = [dtype(dline[i:i + width]) for i in range(0, len(dline), width)]
                parsed_sizes.extend(dline)

            if len(parsed_sizes) != len(_size_keys):
                raise IndexError("AMBER Read: The size of the FLAG POINTERS does not match the FLAG NAMES.")

            for k, v in zip(_size_keys, parsed_sizes):
                sizes[k] = v

            found_sizes = True

    if not found_sizes:
        raise KeyError("AMBER Read: Did not find FLAG POINTERS data.")

    if not found_version:
        raise KeyError("AMBER Read: Did not find VERSION_STAMP data.")

    ### Iterate over the primary data portion of the object

    current_data_category = None
    current_data_type = None
    file_handle = open(filename, "r")
    for x in file_handle:
        if "FLAG ATOM_NAME" in x:

            break

    return  data
    print('--')
    print(pd.read_fwf(file_handle, nrows=3, widths=[80], header=None))
    print('--')
    # print(pd.read_fwf(file_handle, nrows=3, widths=[80], header=None))
    # print('--')
    # print(pd.read_fwf(file_handle, nrows=3, widths=[80], header=None))
    # print('--')
    # print(pd.read_fwf(file_handle, nrows=3, widths=[80], header=None))
    # print('--')
    # print(pd.read_fwf(file_handle, nrows=3, widths=([4] * 20), header=None))
    # print(pd.read_fwf(file_handle, nrows=3, widths=([4] * 20), header=None))
    # print(next(file_handle).strip())

    raise Exception()
    while True:

        # Figure out the size of the chunk to read
        size = _get_size(current_data_category, sizes_dict)
        dl_func = _get_dl_function(current_data_category)
        df_cols = _get_df_columns(current_data_category)

        # Read in the data, in chunks
        remaining = size
        for block in range(int(math.ceil(size / blocksize))):

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


    file_handle.close()
    return data
