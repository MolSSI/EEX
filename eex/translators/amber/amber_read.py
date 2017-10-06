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

import eex
import logging

# AMBER local imports
from . import amber_metadata as amd

logger = logging.getLogger(__name__)


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
    df.set_index(df_index_name, inplace=True)
    return df


def read_amber_file(dl, filename, inpcrd=None, blocksize=5000):
    """

    Parameters
    ----------
    dl : eex.DataLayer
        The datalayer to add data to
    filename : str
        The name of the prmtop file
    inpcrd : str, optional
        If None, attempts to read the file filename.replace("prmtop", "inpcrd") otherwise passes.


    """

    ### First we need to figure out system dimensions
    max_rows = 100  # How many lines do we attempt to search?
    header_data = eex.utility.read_lines(filename, max_rows)

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

            if len(parsed_sizes) != len(amd.size_keys):
                raise IndexError("AMBER Read: The size of the FLAG POINTERS does not match the FLAG NAMES.")

            for k, v in zip(amd.size_keys, parsed_sizes):
                sizes_dict[k] = v

            found_sizes = True

    if not found_sizes:
        raise KeyError("AMBER Read: Did not find FLAG POINTERS data.")

    if not found_version:
        raise KeyError("AMBER Read: Did not find VERSION_STAMP data.")

    ### Iterate over the primary data portion of the object

    # Figure out the size of each label
    label_sizes = {}
    for k, v in amd.data_labels.items():
        if isinstance(v[0], int):
            label_sizes[k] = v[0]
        elif v[0] in list(sizes_dict):
            label_sizes[k] = sizes_dict[v[0]]
        else:
            # print("%30s %40s %d" % (k, v[0], int(eval(v[0], sizes_dict))))
            label_sizes[k] = int(eval(v[0], sizes_dict))

    # print(label_sizes)
    # Find the start
    current_data_category = None
    current_data_type = None
    file_handle = open(filename, "r")
    counter = 0
    for line in file_handle:
        if counter > 100:
            raise KeyError("AMBER Read: Could not find the line a data category in the first 100 lines.")
        matched, name = eex.utility.line_fuzzy_list(line, list(amd.data_labels))
        if matched:
            current_data_category = name
            current_data_type = _parse_format(next(file_handle))

            break
        counter += 1

    # Build any required temporaries
    _current_topology_indices = {
        "bonds": [3, 0, np.array([])],
        "angles": [4, 0, np.array([])],
        "dihedrals": [5, 0, np.array([])],
    }

    # Iterate over the file
    while True:

        # Type out the sizes and types
        nsize = label_sizes[current_data_category]
        nrows = int(math.ceil(nsize / float(current_data_type[0])))
        dtypes = [current_data_type[1]] * current_data_type[0]
        widths = [current_data_type[2]] * current_data_type[0]

        # Read in the data, in chunks
        remaining_read = nrows
        num_blocks = int(math.ceil(nrows / float(blocksize)))
        category_index = 1
        for block in range(num_blocks):

            # Figure out the size of the read
            read_size = blocksize
            if remaining_read < blocksize:
                read_size = remaining_read

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
            if current_data_category in list(amd.atom_property_names):

                df = _data_flatten(data, amd.atom_property_names[current_data_category], category_index, "atom_index")
                category_index += df.shape[0]

                # Add the data to DL
                dl.add_atoms(df, by_value=True, utype=amd.atom_data_units)

            elif current_data_category in amd.store_other:
                df = _data_flatten(data, current_data_category, category_index, "res_index")

                # Force residue pointer to be type int
                if current_data_category == "RESIDUE_POINTER":
                    df = df.astype(int)
                category_index += df.shape[0]

                # Add the data to DL
                dl.add_other(current_data_category, df)

            # Store bond, angle, dihedrals
            elif current_data_category in list(amd.topology_store_names):
                category = current_data_category.split("_")[0].lower()

                mod_size, current_size, remaining_data = _current_topology_indices[category]

                data = data.values.ravel()
                data = data[np.isfinite(data)]

                # Prepend remaining data
                if remaining_data.size:
                    data = np.hstack((remaining_data, data))

                # Get current remaining and
                remaining = data.size % mod_size
                if remaining:
                    data = data[:-remaining]
                    _current_topology_indices[category][-1] = data[-remaining:].copy()
                else:
                    _current_topology_indices[category][-1] = np.array([])

                data = data.reshape(-1, mod_size).astype(int)

                # Weird AMBER indexing, we have: atom1, atom2, ..., term_index
                # Atom indices are (index / 3 + 1)
                data[:, :-1] = data[:, :-1] / 3 + 1

                # Build column names and atom sizes
                col_name = ["atom" + str(x) for x in range(1, mod_size)]
                col_name.append("term_index")

                index = np.arange(current_size, data.shape[0] + current_size, dtype=int)

                # Build and curate the data
                df_dict = {}
                for num, name in enumerate(col_name):
                    df_dict[name] = data[:, num]

                # Form the DF and add!
                df = pd.DataFrame(df_dict, index=index)
                dl.add_terms(mod_size - 1, df)

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
            remaining_read -= blocksize

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

        matched, current_data_category = eex.utility.line_fuzzy_list(category_line, list(amd.data_labels))
        if not matched:
            raise KeyError("AMBER Read: Data category '%s' not understood" % category_line)

        current_data_type = _parse_format(format_line)

    # Close out the file handle
    file_handle.close()

    ### Handle any data we added to the other columns

    # Expand residue values
    res_df = dl.get_other(amd.residue_store_names)

    sizes = np.diff(res_df["RESIDUE_POINTER"])
    last_size = sizes_dict["NATOM"] - res_df["RESIDUE_POINTER"].iloc[-1] + 1
    sizes = np.concatenate((sizes, [last_size])).astype(np.int)

    res_df["residue_index"] = np.arange(res_df.shape[0])
    res_df = pd.DataFrame({
        "residue_index": np.repeat(res_df["residue_index"].values, sizes, axis=0),
        "residue_name": np.repeat(res_df["RESIDUE_LABEL"].values.astype('str'), sizes, axis=0)
    })

    res_df.index = np.arange(1, res_df.shape[0] + 1)
    res_df.index.name = "atom_index"
    dl.add_atoms(res_df, by_value=True)

    # Handle term parameters
    other_tables = dl.list_other_tables()
    for key, param_data in amd.forcefield_parameters.items():
        param_col_names = list(param_data["column_names"])

        # No data to store
        if len(set(param_col_names) - set(other_tables)):
            continue

        cnt = 1  # Start counting from one
        for ind, row in dl.get_other(param_col_names).iterrows():
            params = {}
            for k, v in param_data["column_names"].items():
                params[v] = row[k]
            uid = dl.add_parameters(
                param_data["order"], param_data["form"], params, uid=cnt, utype=param_data["units"])
            cnt += 1

    ### Try to pull in an inpcrd file for XYZ coordinates
    inpcrd_file = filename.replace('.prmtop', '.inpcrd')
    try:
        header_data = eex.utility.read_lines(inpcrd_file, 2)
    except OSError:
        header_data = []

    if header_data:
        inpcrd_size = header_data[1].split()
        if len(inpcrd_size) > 1:
            raise Exception("Cannot handle velocities or pressure in INPCRD file yet.")

        read_size = math.ceil(float(inpcrd_size[0]) / 2)

        file_handle = open(inpcrd_file, "r")
        data = pd.read_fwf(
            file_handle, nrows=read_size, widths=([12] * 6), dtypes=([float] * 6), header=None, skiprows=2)
        file_handle.close()

        df = pd.DataFrame(data.values.reshape(-1, 3), columns=["X", "Y", "Z"])
        df.dropna(axis=0, how="any", inplace=True)
        df.index = np.arange(1, df.shape[0] + 1)

        df.index.name = "atom_index"
        dl.add_atoms(df, utype={"XYZ": "angstrom"})

    return ret_data
