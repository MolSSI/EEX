"""
AMBER EEX I/O
"""

import time
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

        # Get size_information - read in "FLAG POINTERS"
        elif "FLAG POINTERS" in line:
            ncols, dtype, width = amd.parse_format(header_data[num + 1])

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
        # If section always has same size
        if isinstance(v[0], int):
            label_sizes[k] = v[0]
        # Section is variable size (eg NATOMS)
        elif v[0] in list(sizes_dict):
            label_sizes[k] = sizes_dict[v[0]]
        # Section is variable size and has to be evaluated, (eg "(NTYPES * (NTYPES + 1)) / 2")
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
        matched, name = eex.utility.fuzzy_list_match(line, list(amd.data_labels))
        if matched:
            current_data_category = name
            current_data_type = amd.parse_format(next(file_handle))

            break
        counter += 1

    # Build any required temporaries
    _current_topology_indices = {
        "bonds": [3, 0, np.array([])],
        "angles": [4, 0, np.array([])],
        "dihedrals": [5, 0, np.array([])],
    }
    _nonbonded_params = {}

    # Iterate over the file
    while True:
        # Type out the sizes and types
        nsize = label_sizes[current_data_category]
        nrows = int(math.ceil(nsize / float(current_data_type[0])))
        dtypes = [current_data_type[1]] * current_data_type[0]
        widths = [current_data_type[2]] * current_data_type[0]
        # print(current_data_category)

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
            tmp_handle = None
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

            # Store forcefield parameters as "other" for later processing
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
                # A negative indexed value in position 3 indicates 1-4 NB interactions for this dihedral
                # should not be counted (multi-term dihedral or cyclic system). Store as
                # dihedral for now, neglecting negative sign

                data = np.absolute(data)

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

            # Get box information from prmtop if here. Will be overwritten by inpcrd if information is provided.
            elif current_data_category == "BOX_DIMENSIONS":
                box_size = {}
                box_center = []
                a = data[1].values[0]
                b = data[2].values[0]
                c = data[3].values[0]

                box_size["alpha"] = data[0].values[0]
                box_size["beta"] = data[0].values[0]
                box_size["gamma"] = data[0].values[0]



                for v in amd.box_units["center"]:
                    box_center.append(eval(v))

                box_center = dict(zip(['x', 'y', 'z'], box_center))

                box_size["a"] = a
                box_size["b"] = b
                box_size["c"] = c

                dl.set_box_size(box_size, utype={"a": amd.box_units["length"], "b": amd.box_units["length"],
                                                 "c" : amd.box_units["length"], "alpha": amd.box_units["angle"],
                                                 "beta": amd.box_units["angle"], "gamma": amd.box_units["angle"],})

                dl.set_box_center(box_center, utype={"x": amd.box_units["length"], "y": amd.box_units["length"],
                                                 "z" : amd.box_units["length"]})

            else:
                # logger.debug("Did not understand data category.. passing")
                pass
            # elif current_data_category == "ATOM_NAME":

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

        matched, current_data_category = eex.utility.fuzzy_list_match(category_line, list(amd.data_labels))
        if not matched:
            raise KeyError("AMBER Read: Data category '%s' not understood" % category_line)

        current_data_type = amd.parse_format(format_line)

    # Close out the file handle
    file_handle.close()

    ### Handle any data we added to the other columns

    # Expand residue values
    res_df = dl.get_other(amd.residue_store_names)

    sizes = np.diff(res_df["RESIDUE_POINTER"])
    last_size = sizes_dict["NATOM"] - res_df["RESIDUE_POINTER"].iloc[-1] + 1
    sizes = np.concatenate((sizes, [last_size])).astype(np.int)

    res_df["residue_index"] = np.arange(1, res_df.shape[0] + 1)
    res_df = pd.DataFrame({
        "residue_index": np.repeat(res_df["residue_index"].values, sizes, axis=0),
        "residue_name": np.repeat(res_df["RESIDUE_LABEL"].values.astype('str'), sizes, axis=0)
    })

    res_df.index = np.arange(1, res_df.shape[0] + 1)
    res_df.index.name = "atom_index"
    dl.add_atoms(res_df, by_value=True)

    # Expand molecule values - if periodic simulation information will be given. If not, we assign all atoms the same
    # molecule number.
    if sizes_dict["IFBOX"] > 0:
        molecule_df = dl.get_other(amd.molecule_store_names)

        # This is only set if SOLVENT_POINTERS is present in the prmtop (ie - a periodic simulation)
        number_of_molecules = label_sizes["ATOMS_PER_MOLECULE"]

        molecule_range = np.arange(1, number_of_molecules+1)

        # Next, we need to create a dataframe with the column header "molecule_index". The variable
        # molecule_df contains a list of the number of atoms in each molecule, while the variable number_of_molecules
        # gives the total number of molecules in the system.

        molecule_df = pd.DataFrame({
            "molecule_index": np.repeat(molecule_range, [int(x) for x in molecule_df["ATOMS_PER_MOLECULE"].values], axis=0)
        })

        molecule_df.index = np.arange(1, molecule_df.shape[0] +1)
        molecule_df.index.name = "atom_index"

    else:
        molecule_df = pd.DataFrame({
            "molecule_index": np.ones(sizes_dict["NATOM"])
        })

        molecule_df.index = np.arange(1, molecule_df.shape[0] +1)
        molecule_df.index.name = "atom_index"

    dl.add_atoms(molecule_df, by_value=True)



    # Handle forcefield parameters
    other_tables = set(dl.list_other_tables())
    for key, param_data in amd.forcefield_parameters.items():
        param_col_names = list(param_data["column_names"])
        # No data to store
        if len(set(param_col_names) - other_tables):
            continue

        # Bond parameters (bond, angle, dihedral) will have an "order", the order for nonbond parameters is None
        if param_data["order"] is not None:
            cnt = 1  # Start counting from one
            for ind, row in dl.get_other(param_col_names).iterrows():
                params = {}
                for k, v in param_data["column_names"].items():
                    params[v] = row[k]
                uid = dl.add_term_parameter(
                    param_data["order"], param_data["form"], params, uid=cnt, utype=param_data["units"])
                cnt += 1
        else:
            # Get info for grabbing LJ parameters
            nb_parm_index = dl.get_other("NONBONDED_PARM_INDEX")
            A_coeff_list = dl.get_other("LENNARD_JONES_ACOEF")
            B_coeff_list = dl.get_other("LENNARD_JONES_BCOEF")
            stored_atom_types = np.unique(dl.get_atoms('atom_type'))
            ntypes = len(stored_atom_types)

            # Need relevant atom types. Should go 1...n where n is number of atom_types. Loop through stored to
            # get all combinations

            for x in range(len(stored_atom_types)):
                for y in range(x + 1):
                    # For amber, section NONBOND_PARM_INDEX gives pointer to LENNARD_JONES_ACOEF and _BCOEF sections.
                    # The atom types are used to compute NB_PARM_INDEX index, which is used to get ACOEF and BCOEF

                    # Get atom types - nb_key is (atom_type1, atom_type2)
                    nb_key = (stored_atom_types[y], stored_atom_types[x])

                    # Calcuate index into NB_PARM_INDEX
                    ind_nb_parm_index = ntypes * (nb_key[0] - 1) + nb_key[1]

                    # Get index into LENNARDJONES_ACOEF and LENNARDJONES_BCOEF
                    # Subtract 1 because Amber indexes from 1, but python indexes from 0
                    nb_index = nb_parm_index.iloc[ind_nb_parm_index - 1]

                    # Grab values
                    A_coeff = A_coeff_list.iloc[nb_index - 1]['LENNARD_JONES_ACOEF'].values[0]
                    B_coeff = B_coeff_list.iloc[nb_index - 1]['LENNARD_JONES_BCOEF'].values[0]

                    # Store in datalayer
                    dl.add_nb_parameter(
                        atom_type=nb_key[0],
                        atom_type2=nb_key[1],
                        nb_parameters={"A": A_coeff,
                                       "B": B_coeff},
                        nb_name=amd.forcefield_parameters["nonbond"]["form"]["name"],
                        nb_model=amd.forcefield_parameters["nonbond"]["form"]["form"],
                        utype=amd.forcefield_parameters["nonbond"]["units"])

    ### Try to pull in an inpcrd file for XYZ coordinates and box information
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

        ## Read in extra line if simulation is periodic - NOTE - this will not work for files with velocity information
        ## (mentioned above in exception)
        if sizes_dict["IFBOX"] > 0:
            read_size += 1

        file_handle = open(inpcrd_file, "r")
        data = pd.read_fwf(
            file_handle, nrows=read_size, widths=([12] * 6), dtypes=([float] * 6), header=None, skiprows=2)

        file_handle.close()

        if sizes_dict["IFBOX"] > 0:
            box_information = data.tail(1).values[0]

            box_sizes = {"a": box_information[0], "b": box_information[1], "c": box_information[2],
                         "alpha": box_information[3], "beta": box_information[3], "gamma": box_information[3],
                         }

            dl.set_box_size(box_sizes, utype={"a": amd.box_units["length"], "b": amd.box_units["length"],
                                                 "c" : amd.box_units["length"], "alpha": amd.box_units["angle"],
                                                 "beta": amd.box_units["angle"], "gamma": amd.box_units["angle"],})

            # Drop box info from atom coordinates
            data.drop(data.index[-1], inplace=True)


        df = pd.DataFrame(data.values.reshape(-1, 3), columns=["X", "Y", "Z"])
        df.dropna(axis=0, how="any", inplace=True)
        df.index = np.arange(1, df.shape[0] + 1)

        df.index.name = "atom_index"
        dl.add_atoms(df, utype={"XYZ": "angstrom"})


    return ret_data
