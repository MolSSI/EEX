"""
LAMMPS EEX I/O
"""
import os
import pandas as pd
import math

import eex

from . import lammps_metadata as lmd

import logging
logger = logging.getLogger(__name__)

def read_lammps_data_file(dl, filename, blocksize=110):
    
    ### Figure out system dimensions and general header data
    max_rows = 100  # How many lines do we attempt to search?
    header_data = eex.utility.read_lines(filename, max_rows)

    box_size = {}
    sizes_dict = {}

    startline = None
    current_data_category = None
    category_list = lmd.build_valid_category_list()

    header = header_data[0]
    for num, line in enumerate(header_data[1:]):

        # Skip blanklines
        if line == "":
            continue

        # Skip comment line
        elif line[0] == "#":
            continue

        # We are
        elif eex.utility.fuzzy_list_match(line, category_list)[0]:
            startline = num + 3  # Skips first row and two blank lines
            current_data_category = eex.utility.fuzzy_list_match(line, category_list)[1]
            break

        # Figure out the dims
        elif ("lo" in line) and ("hi" in line):
            dline = line.split()
            if dline[-1] == "xhi":
                box_size["x"] = (float(dline[0]), float(dline[1]))

            elif dline[-1] == "yhi":
                box_size["y"] = (float(dline[0]), float(dline[1]))
            elif dline[-1] == "zhi":
                box_size["z"] = (float(dline[0]), float(dline[1]))
            else:
                raise KeyError(
                    "LAMMPS Read: The following line looks like a dimension line, but does not match:\n%s" % line)

        # Are we a size line?
        elif eex.utility.fuzzy_list_match(line, lmd.size_keys)[0]:
            dline = line.split()
            size = int(dline[0])
            size_name = " ".join(dline[1:])

            if size_name in list(sizes_dict):
                raise KeyError("LAMMPS Read: KeyError size key %s already found." % size_name)
            elif size_name not in lmd.size_keys:
                raise KeyError("LAMMPS Read: KeyError size key %s not recognized." % size_name)
            else:
                sizes_dict[size_name] = size

        else:
            raise IOError("LAMMPS Read: Line not understood!\n%s" % line)

    # Set the box size
    dl.set_box_size(box_size, utype=lmd.get_context("real", "[length]"))

    # Make sure we have what we need
    if startline is None:
        raise IOError("LAMMPS Read: Did not find data start in %d header lines." % max_rows)

    if ("atoms" not in list(sizes_dict)) or ("atom types" not in list(sizes_dict)):
        raise IOError("LAMMPS Read: Did not find size data on 'atoms' or 'atom types' in %d header lines." % max_rows)

    ### Create temporaries specific to the current unit specification
    op_table = lmd.build_operation_table("real", sizes_dict)
    term_table = lmd.build_term_table("real")

    # term_table = {"Bond Coeffs": {"order": 2, "name":"harmonic", "utype":""}, "Angle Coeffs": {}}

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
        op = op_table[current_data_category]

        # Read in the data, in chunks
        remaining = op["size"]
        num_blocks = int(math.ceil(op["size"] / float(blocksize)))
        for block in range(num_blocks):

            # Figure out the size of the read
            read_size = blocksize
            if remaining < blocksize:
                read_size = remaining

            # Read and update DL
            data = reader.get_chunk(read_size).dropna(axis=1, how="all")

            # Nothing defined
            if op["dl_func"] == "NYI":
                pass

            # Single call
            elif op["call_type"] == "single":
                # print(op["dl_func"], op["df_cols"])
                if "df_cols" in op:
                    data.columns = op["df_cols"]
                dl.call_by_string(op["dl_func"], data, **op["kwargs"])

            elif op["call_type"] == "add_atom_parameters":
                atom_prop = op["atom_property"]
                utype = op["kwargs"]["utype"][atom_prop]
                for idx, row in data.iterrows():
                    dl.add_atom_parameter(atom_prop, row.iloc[1], uid=row.iloc[0], utype=utype)

            # Adding parameters
            elif op["call_type"] == "parameter":
                order = op["args"]["order"]
                fname = op["args"]["form_name"]
                cols = term_table[order][fname]["parameters"]
                data.columns = ["uid"] + cols

                for idx, row in data.iterrows():
                    params = list(row[cols])
                    utype = term_table[order][fname]["utype"]
                    dl.add_term_parameter(order, fname, params, uid=int(row["uid"]), utype=utype)

            else:
                raise KeyError("Operation table call '%s' not understoop" % op["call_type"])

            # Update remaining
            remaining -= blocksize

        # Figure out the next category to read
        try:
            tmp = reader.get_chunk(1).dropna(axis=1, how="any")
        except StopIteration:
            break

        current_data_category = " ".join(str(x) for x in list(tmp.iloc[0]))

    # Mass is missing its index, we can copy the data over
    dl.store.copy_table("atom_type", "mass", {"atom_type": "mass"})

    # raise Exception("")
    data = {}
    data["sizes"] = sizes_dict
    data["header"] = header

    return data

def get_bond_coeff():
    pass
def get_angle_coeff():
    pass
def get_diherdal_coeff():
    pass
def get_include():
    pass
def get_variable():
    pass
def get_units(opts):
    pass
def get_atom_style():
    pass
def get_pair_style():
    pass
def kspace_style():
    pass
def pair_modify():
    pass
def special_bonds():
    pass
def bond_style():
    pass
def angle_style():
    pass
def dihedral_style():
    pass

keyword_dispatcher = {
    "read_data": read_lammps_data_file,
    #"bond_coeff": get_bond_coeff,
    #"angle_coeff": get_angle_coeff, 
    #"dihedral_coeff": get_dihedral_coeff, 
    #"include": get_include
    #"variable": get_variable, 
    "units": get_units, 
    #"atom_style": get_atom_style, 
    #"pair_style": get_pair_style, 
    #"kspace_style": get_kspace_style, 
    #"pair_modify": get_pair_modify, 
    #"special_bonds": get_special_bonds, 
    #"bond_style": get_bond_style, 
    #"angle_style": get_angle_style, 
    #"dihedral_style": get_dihedral_style, 
    }

def read_lammps_file(dl, fname):
    """
        Reads a LAMMPS input file
    """
    input_dir = os.path.dirname(fname)
    with open(fname, 'r') as input_file:
        for line in input_file:
            if line.strip():
                line = line.split()
                keyword = line[0]
                keyword_opts = line[1:]
                if keyword in keyword_dispatcher:
                    f = keyword_dispatcher[keyword]
                    if keyword == "read_data":
                        data = f(dl, os.path.join(input_dir, keyword_opts[0]))
    return data
