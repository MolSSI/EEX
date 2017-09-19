"""
LAMMPS EEX I/O
"""

import pandas as pd
import math
import re
import numpy as np

# Python 2/3 compat
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

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
    "ATOM_NAME": ["NATOM"],
    "CHARGE": ["NATOM"],
    "ATOMIC_NUMBER": ["NATOM"],
    "MASS": ["NATOM"],
    "ATOM_TYPE_INDEX": ["NATOM"],
    "NUMBER_EXCLUDED_ATOMS": ["NATOM"],
    "NONBONDED_PARM_INDEX": ["NTYPES ** 2"],
    "RESIDUE_LABEL": ["NRES"],
    "RESIDUE_POINTER": ["NRES"],
    "BOND_FORCE_CONSTANT": ["NUMBND"],
    "BOND_EQUIL_VALUE": ["NUMBND"],
    "ANGLE_FORCE_CONSTANT": ["NUMANG"],
    "ANGLE_EQUIL_VALUE": ["NUMANG"],
    "DIHEDRAL_FORCE_CONSTANT": ["NPTRA"],
    "DIHEDRAL_PERIODICITY": ["NPTRA"],
    "DIHEDRAL_PHASE": ["NPTRA"],
    "SCEE_SCALE_FACTOR": ["NPTRA"],
    "SCNB_SCALE_FACTOR": ["NPTRA"],
    "SOLTY": ["NATYP"],
    "LENNARD_JONES_ACOEF": ["(NTYPES * (NTYPES + 1)) / 2"],
    "LENNARD_JONES_BCOEF": ["(NTYPES * (NTYPES + 1)) / 2"],
    "BONDS_INC_HYDROGEN": ["3 * NBONH"],
    "BONDS_WITHOUT_HYDROGEN": ["3 * NBONA"],
    "ANGLES_INC_HYDROGEN": ["4 * NTHETH"],
    "ANGLES_WITHOUT_HYDROGEN": ["4 * NTHETA"],
    "DIHEDRALS_INC_HYDROGEN": ["5 * NPHIH"],
    "DIHEDRALS_WITHOUT_HYDROGEN": ["5 * NPHIA"],
    "EXCLUDED_ATOMS_LIST": ["NNB"],
    "HBOND_ACOEF": ["NPHB"],
    "HBOND_BCOEF": ["NPHB"],
    "HBCUT": ["NPHB"],
    "AMBER_ATOM_TYPE": ["NATOM"],
    "TREE_CHAIN_CLASSIFICATION": ["NATOM"],
    "JOIN_ARRAY": ["NATOM"],  # This section is filled with zeros, but is unused. We should not store it.
    "IROTAT": ["NATOM"],
    "SOLVENT_POINTERS": ["3 if IFBOX else 0"],
    "ATOMS_PER_MOLECULE": ["NATOM"],
    # "ATOMS_PER_MOLECULE": ["SOLVENT_POINTERS[1] if IFBOX else 0"], # SOLVENT_POINTERS[1] == NPSM
    "BOX_DIMENSIONS": [4],
    "CAP_INFO": ["1 if IFCAP else 0"],
    "CAP_INFO2": ["4 if IFCAP else 0"],
    "RADIUS_SET": [1],
    "RADII": ["NATOM"],
    "SCREEN": ["NATOM"],
    "IPOL": [1],
    "POLARIZABILITY": [0]
    # "POLARIZABILITY": ["NATOM if IPOL else 0"]
}

_data_units = {
    """
    Gives units of sections with units - for conversions NYI

    """

    "CHARGE": ["18.2223*e"],  # Internal units
    "MASS": ["g mol^-1"],
    "BOND_FORCE_CONSTANT": ["kcal mol^-1 Angstrom^-2"],
    "BOND_EQUIL_VALUE": ["Angstrom"],
    "ANGLE_FORCE_CONSTANT": ["kcal mol^-2 radian^2"],
    "ANGLE_EQUIL_VALUE": ["radian"],
    "DIHEDRAL_FORCE_CONSTANT": ["kcal mol^-1"],
    "DIHEDRAL_PHASE": ["radian"],
    "LENNARD_JONES_ACOEFF": ["kcal mol^-12"],
    "LENNARD_JONES_BCOEFF": ["kcal mol^-6"],
}

_atom_property_names = {"ATOM_NAME": "atom_name",
                        "CHARGE": "charge",
                        "MASS": "mass",
                        "ATOM_TYPE_INDEX": "atom_type",
                        "ATOMIC_NUMBER": "atomic_number",
                        # "AMBER_ATOM_TYPE": "atom_type_name",
                        # "RADII" : "implicit_solvent_radius",
                        }

_residue_store_names = ["RESIDUE_LABEL", "RESIDUE_POINTER"]

_interaction_store_names = ["BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"]

_forcefield_parameters = ["BOND_FORCE_CONSTANT", "BOND_EQUIL_VALUE", "ANGLE_FORCE_CONSTANT", "ANGLE_EQUIL_VALUE",
                          "DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PERIODICITY", "DIHEDRAL_PHASE", "LENNARD_JONES_ACOEFF",
                          "LENNARD_JONES_BCOEFF", ]


def _parse_format(string):
    """
    Parses an AMBER style format string.

    Example:

    >>> string = "%FORMAT(10I8)"
    >>> _parse_format(string)
    [10, int, 8]
    """
    if "FORMAT" not in string:
        raise ValueError("AMBER: Did not understand format line '%s'." % string)

    pstring = string.replace("%FORMAT(", "").replace(")", "").strip()
    ret = [x for x in re.split('(\d+)', pstring) if x]

    if ret[1] == "I":
        if len(ret) != 3:
            raise ValueError("AMBER: Did not understand format line '%s'." % string)
        ret[1] = int
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
    elif ret[1] == "E":
        if len(ret) != 5:
            raise ValueError("AMBER: Did not understand format line '%s'." % string)
        ret[1] = float
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
        # The .8 is not interesting to us
    elif ret[1] == "a":
        if len(ret) != 3:
            raise ValueError("AMBER: Did not understand format line '%s'." % string)
        ret[1] = str
        ret[0] = int(ret[0])
        ret[2] = int(ret[2])
    else:
        raise ValueError("AMBER: Type symbol '%s' not understood from line '%s'." % (ret[1], string))

    return ret


def _data_flatten(data, column_name, category_index, df_index_name):
    # Reorganize the data 2D -> 1D packing
    flat_data = data.values.flatten()
    index = np.arange(category_index, flat_data.shape[0] + category_index)
    dl_col_name = column_name

    # Build and curate the data
    df = pd.DataFrame({df_index_name: index, dl_col_name: flat_data})
    df.dropna(axis=0, how="any", inplace=True)
    return df

def _data_reshape(data, num_columns):
    data_values = data.values
    # Remove nans
    data_process = data_values[~np.isnan(data_values)]
    # Reshape data
    data_shape = data_process.reshape(-1,num_columns)
    df = pd.DataFrame(data=data_shape)
    return df



def read_amber_file(dl, filename, blocksize=5000):
    ### First we need to figure out system dimensions
    max_rows = 100  # How many lines do we attempt to search?
    with open(filename, "r") as infile:
        header_data = [next(infile).strip() for x in range(max_rows)]

    sizes_dict = {}
    ret_data = {}
    found_sizes = False
    found_version = False
    for num, line in enumerate(header_data):

        # First try to figure out the version
        if "VERSION" in line:
            sline = line.strip().split()
            if "VERSION_STAMP" == sline[1]:
                ret_data["VERSION"] = sline[3]
                if ret_data["VERSION"] != "V0001.000":
                    raise ValueError("AMBER Read: Did not recognize version '%s'." % ret_data["VERSION"])
                found_version = True
            else:
                raise ValueError("AMBER Read: Could not understand version line.")

        # Get size_information
        elif "FLAG POINTERS" in line:
            ncols, dtype, width = _parse_format(header_data[num + 1])

            parsed_sizes = []
            for shift in range(4):
                dline = header_data[num + 2 + shift]
                # Make sure the leading white space is not being cut off each line
                if shift < 3:
                    dline = dline.rjust(80)
                dline = [dtype(dline[i:(i + width)]) for i in range(0, len(dline), width)]
                parsed_sizes.extend(dline)

            if len(parsed_sizes) != len(_size_keys):
                raise IndexError("AMBER Read: The size of the FLAG POINTERS does not match the FLAG NAMES.")

            for k, v in zip(_size_keys, parsed_sizes):
                sizes_dict[k] = v

            found_sizes = True

    if not found_sizes:
        raise KeyError("AMBER Read: Did not find FLAG POINTERS data.")

    if not found_version:
        raise KeyError("AMBER Read: Did not find VERSION_STAMP data.")

    ### Iterate over the primary data portion of the object

    # Figure out the size of each label
    label_sizes = {}
    for k, v in _data_labels.items():
        if isinstance(v[0], int):
            label_sizes[k] = v[0]
        elif v[0] in list(sizes_dict):
            label_sizes[k] = sizes_dict[v[0]]
        else:
            # print("%30s %40s %d" % (k, v[0], int(eval(v[0], sizes_dict))))
            label_sizes[k] = int(eval(v[0], sizes_dict))

    # Find the start
    current_data_category = None
    current_data_type = None
    file_handle = open(filename, "r")
    counter = 0
    for line in file_handle:
        if counter > 100:
            raise KeyError("AMBER Read: Could not find the line a data category in the first 100 lines.")
        matched, name = utility.line_fuzzy_list(line, list(_data_labels))
        if matched:
            current_data_category = name
            current_data_type = _parse_format(next(file_handle))

            break
        counter += 1

    while True:

        # Type out the sizes and types
        nsize = label_sizes[current_data_category]
        nrows = int(math.ceil(nsize / float(current_data_type[0])))
        dtypes = [current_data_type[1]] * current_data_type[0]
        widths = [current_data_type[2]] * current_data_type[0]

        # Read in the data, in chunks
        remaining = nrows
        # print(remaining)
        # print(current_data_category, nsize, nrows, current_data_type)
        num_blocks = int(math.ceil(nrows / float(blocksize)))
        category_index = 0
        for block in range(num_blocks):

            # Figure out the size of the read
            read_size = blocksize
            if remaining < blocksize:
                read_size = remaining

            # read_fwf will push the file pointer *at least 4* so lets just read it in
            if read_size < 4:
                data = []
                for x in range(read_size):
                    data.append(next(file_handle))
                    tmp_handle = StringIO("".join(data))
            else:
                tmp_handle = file_handle

            # Read in the data
            data = pd.read_fwf(tmp_handle, nrows=read_size, widths=widths, dtypes=dtypes, header=None)

            # 1D atom properties
            if current_data_category in list(_atom_property_names):

                df = _data_flatten(data, _atom_property_names[current_data_category], category_index, "atom_index")
                category_index += df.shape[0]
                # Add the data to DL
                dl.add_atoms(df, by_value=True)

            elif current_data_category in list(_residue_store_names):
                df = _data_flatten(data, current_data_category, category_index, "res_index")

                # Force residue pointer to be type int
                if current_data_category == "RESIDUE_POINTER":
                    df = df.astype(int)
                category_index += df.shape[0]

                # Add the data to DL
                dl.add_other(current_data_category, df)

            elif current_data_category in list(_interaction_store_names):
                df = _data_flatten(data, "Bonds", category_index, "index")
                dl.add_other("bonds", df.astype(int, inplace=True))

            else:
                # logger.debug("Did not understand data category.. passing")
                pass
            # elif current_data_category == "ATOM_NAME":

            # print(data.head())
            # print(data.tail())
            # if dl_func != "NYI":
            #     # print(data)
            #     data.columns = df_cols
            # dl.call_by_string(dl_func, data)

            # Update remaining
            remaining -= blocksize

        # If we are doing nothing, we still have a blank line
        if nrows == 0:
            next(file_handle)

        # Figure out the next category to read
        try:
            category_line = next(file_handle)
            format_line = next(file_handle)
        except StopIteration:
            break

        # Bad hack for solvent pointers
        if "SOLVENT_POINTERS" in category_line:
            dline = next(file_handle).rstrip()
            data = [int(dline[i:i + 8]) for i in range(0, len(dline), 8)]
            label_sizes["ATOMS_PER_MOLECULE"] = data[1]
            try:
                category_line = next(file_handle)
                format_line = next(file_handle)
            except StopIteration:
                break

        matched, current_data_category = utility.line_fuzzy_list(category_line, list(_data_labels))
        if not matched:
            raise KeyError("AMBER Read: Data category '%s' not understood" % category_line)

        current_data_type = _parse_format(format_line)
        # break
        # raise Exception("")

    ### Handle any data we added to the other columns

    # Expand residue values

    res_df = dl.get_other(_residue_store_names)

    sizes = np.diff(res_df["RESIDUE_POINTER"])
    last_size = sizes_dict["NATOM"] - res_df["RESIDUE_POINTER"].iloc[-1] + 1
    sizes = np.concatenate((sizes, [last_size])).astype(np.int)

    res_df["residue_index"] = np.arange(0, res_df.shape[0])
    res_df = pd.DataFrame({"residue_index": np.repeat(res_df["residue_index"].values, sizes, axis=0),
                           "residue_name": np.repeat(res_df["RESIDUE_LABEL"].values.astype('str'), sizes, axis=0)})

    res_df.index.name = "atom_index"
    dl.add_atoms(res_df, by_value=True)


    # Handle bonds

    # Get stored bond data from data layer
    bond_df = dl.get_other("bonds")

    # Reshape data
    bond_reshape = _data_reshape(bond_df, 3)

    # Calculate atom indices for bonds based on internal amber format (see ambermd.org/prmtop.pdf)
    bond_reshape.loc[:,0:1] = (bond_reshape.loc[:,0:1]/3+1).astype(int)

    # Add names to columns
    bond_reshape.columns = ["atom1_index", "atom2_index", "bond_type"]
    bond_reshape['bond_index'] = bond_reshape.index

    # Add bonds to data layer
    dl.add_bonds(bond_reshape)

    file_handle.close()
    return ret_data
