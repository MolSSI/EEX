"""
LAMMPS EEX I/O
"""
import os
import pandas as pd
import math
import re
import eex

from . import lammps_metadata as lmd

import logging
logger = logging.getLogger(__name__)


def read_lammps_data_file(dl, filename, extra_simulation_data, blocksize=110):

    if not isinstance(extra_simulation_data, dict):
        raise TypeError("Validate term dict: Extra simulation data type '%s' not understood" % str(type(utype)))

    # This is a list of keywords needed from the input file. 
    # This list depends on the molecule topology (for instance, ethane does
    # not require 'dihedral_style' information, but butane does. 
    # More keywords will be appended as required when reader the header
    # section
    needed_keywords = ["units"]

    # Figure out system dimensions and general header data
    max_rows = 100  # How many lines do we attempt to search?
    header_data = eex.utility.read_lines(filename, max_rows)

    box_size = {}
    box_center = {}
    tilt_factors = {'xy': 0, 'xz': 0, 'yz': 0}
    sizes_dict = {}

    startline = None
    current_data_category = None


    # Category_list contains keywords of data file e.g. Atoms, Masses, 
    # Bond Coeffs, etc.
    category_list = lmd.build_valid_category_list()

    header = header_data[0]
    for num, line in enumerate(header_data[1:]):

        # Skip blanklines
        if line == "":
            continue

        # Skip comment line
        elif line[0] == "#":
            continue

        # Evaluate if line contains any keyword contained in category list. If
        # True, it means we have read all the required information in the 
        # header section. The variable start line is where the header ends
        # and data begins. At this point, the variable current_data_category 
        # holds the first data section in the data file
        elif eex.utility.fuzzy_list_match(line, category_list)[0]:
            startline = num + 3  # Skips first row and two blank lines
            current_data_category = eex.utility.fuzzy_list_match(line, category_list)[1]
            break

        # Read tilt factors
        elif ("xy" in line) and ("xz" in line) and ("yz" in line):
            dline = line.split()
            tilt_factors["xy"] = float(dline[0])
            tilt_factors["xz"] = float(dline[1])
            tilt_factors["yz"] = float(dline[2])

        # Figure out the dims
        elif ("lo" in line) and ("hi" in line):
            dline = line.split()
            if dline[-1] == "xhi":
                box_size["x"] = float(dline[1]) - float(dline[0])
                box_center["x"] = box_size["x"]/2. + float(dline[0])
            elif dline[-1] == "yhi":
                box_size["y"] = float(dline[1]) - float(dline[0])
                box_center["y"] = box_size["y"] / 2. + float(dline[0])
            elif dline[-1] == "zhi":
                box_size["z"] = float(dline[1]) - float(dline[0])
                box_center["z"] = box_size["z"] / 2. + float(dline[0])

            else:
                raise KeyError(
                    "LAMMPS Read: The following line looks like a dimension line, but does not match:\n%s" % line)

        # Are we a size line?
        # Lmd.size_keys contains ['atoms', 'atom types', 'bonds', 'bond types']
        # located at the top of the file
        # This chunk of code populates sizes_dict
        elif eex.utility.fuzzy_list_match(line, lmd.size_keys)[0]:
            dline = line.split()
            size = int(dline[0])
            size_name = " ".join(dline[1:])

            if size_name in list(sizes_dict):
                raise KeyError("LAMMPS Read: KeyError size key %s already found." % size_name)
            if size_name not in lmd.size_keys:
                raise KeyError("LAMMPS Read: KeyError size key %s not recognized." % size_name)

            sizes_dict[size_name] = size

            # If we found keywords such as 'Angle types' that means we need
            # to provide 'Angle style' in the input file (and in the extra
            # simulation data dictionary). This is important for smaller 
            # molecules such as ethane or propane that do not have Angles or
            # Dihedrals. 

            if 'types' in size_name and size > 0:
                needed_keywords.append('%s_style' %(split_size[0]))

        else:
            raise IOError("LAMMPS Read: Line not understood!\n%s" % line)

    # We have now a list of the needed keywords based on the topology. 
    # Needed_keywords contains at least ["units"], but can also contain
    # ["units", "bond_style", "angle_style", ...]

    for keyword in needed_keywords:
        if keyword not in extra_simulation_data.keys():
            raise KeyError("The keyword '%s' is missing" % keyword)
        else:
            if keyword == 'units':
                unit_style = extra_simulation_data["units"]
                if unit_style not in lmd.units_style:
                    raise KeyError("Could not find unit style '%s'." % unit_style)
            elif keyword == 'bond_style':
                bond_style = extra_simulation_data['bond_style']
                if bond_style not in lmd.lammps_ff.term_data[2]:
                    raise KeyError("Could not find bond style '%s'." % bond_style)
            elif keyword == 'angle_style':
                angle_style = extra_simulation_data['angle_style']
                if angle_style not in lmd.lammps_ff.term_data[3]:
                    raise KeyError("Could not find angle style '%s'." % angle_style)
            elif keyword == 'dihedral_style':
                dihedral_style = extra_simulation_data['dihedral_style']
                if dihedral_style not in lmd.lammps_ff.term_data[4]:
                    raise KeyError("Could not find dihedral style '%s'." % dihedral_style)
            else:
                raise KeyError("Key %s not understood. " % k)

    # Set the box size
    lattice_constants = eex.utility.compute_lattice_constants(box_size, tilt_factors)

    dl.set_box_size(lattice_constants)
    dl.set_box_center(box_center)

    # Make sure we have what we need
    if startline is None:
        raise IOError("LAMMPS Read: Did not find data start in %d header lines." % max_rows)

    if ("atoms" not in list(sizes_dict)) or ("atom types" not in list(sizes_dict)):
        raise IOError("LAMMPS Read: Did not find size data on 'atoms' or 'atom types' in %d header lines." % max_rows)

    # sizes_dict["Pair types"] is relevant if we need to specify 
    # cross interactions explicitly
    sizes_dict["pair types"] = sizes_dict["atom types"] * (sizes_dict["atom types"] + 1) / 2

    # Create temporaries (op_table, term_table), specific to the current unit 
    # specification
    # op_table is a dictionary of dictionaries. Its keys are the keywords
    # of the Lammps data file i.e. 'Atoms', 'Bonds', 'Dihedral coeffs', etc.
    # Each key, i.e. Atoms, has a value that is a dictionary. The keys for
    # these nested dictionaries are 
    # ['size', 'dl_func', 'df_cols', 'kwargs', 'call_type']
    op_table = lmd.build_operation_table(unit_style, sizes_dict)

    # Term table is a dictionary with all the metadata for each two
    # three and four body potential functional forms using the lammps
    # metadata file. The units of these functional forms are consistent with
    # the information in the input file.
    # E.g. term_table[2]["fene"]["form"]
    term_table = lmd.build_term_table(unit_style)

    # nb_term_table is a dictionary of all the different pair_styles
    # with the desired units
    nb_term_table = lmd.build_nb_table(unit_style)

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

        # Read in the current section, in chunks
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
                fname = extra_simulation_data[op["args"]["style_keyword"]]
                cols = term_table[order][fname]["parameters"]
                data.columns = ["uid"] + cols
                for idx, row in data.iterrows():
                    params = list(row[cols])
                    utype = term_table[order][fname]["utype"]
                    dl.add_term_parameter(order, fname, params, uid=int(row["uid"]), utype=utype)

            elif op["call_type"] == "nb_parameter":
                fname = op["kwargs"]["nb_name"]
                fform = op["kwargs"]["nb_model"]
                utype = nb_term_table[fname]["utype"]
                op["kwargs"]["utype"] = utype
                cols = nb_term_table[fname]["parameters"]
                n_required_params = len(cols)
                n_given_params = len(data.columns)
                difference = n_given_params - n_required_params
                if difference != op["n_uids"]:
                    raise Exception("Incorrect number of arguments for pair_coeff")
                uid_list = ["uid" + str(x) for x in range(0, difference)]
                data.columns = uid_list + cols
                for idx, row in data.iterrows():

                    op["kwargs"]["nb_parameters"] = row[['epsilon', 'sigma']].to_dict()
                    if len(uid_list)==1:
                        op["kwargs"]["atom_type"] = int(row['uid0'])
                    elif len(uid_list)==2:
                        op["kwargs"]["atom_type"] = int(row['uid0'])
                        op["kwargs"]["atom_type2"] = int(row['uid1'])
                    dl.call_by_string(op["dl_func"], **op["kwargs"])

                print(dl.list_nb_parameters('LJ'))


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

    ## raise Exception("")
    #data = {}
    #data["sizes"] = sizes_dict
    #data["header"] = header

    #return data

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
def get_atom_style():
    pass
def get_pair_style():
    pass
def kspace_style():
    pass
def pair_modify():
    pass

def read_lammps_input_file(dl, fname, blocksize=110):
    """
        Reads a LAMMPS input file
    """
    extra_simulation_data = {}

    variable_list = {}

    input_dir = os.path.dirname(fname)

    input_file = eex.utility.read_lines(fname)

    for lnum, line in enumerate(input_file):
        if len(line) == 0: continue

        # Handle variables
        if "$" in line:
            m = re.search('\$\{(.*?)\}', line)
            variable_name = m.group(1)
            try:
                variable_values = variable_list[variable_name]["values"]
                # Need to fix this
                line = re.sub('\$\{(.*?)\}', variable_values[0], line)
                
            except KeyError:
                raise KeyError("The variable %s has not been defined" % variable_name)

        line_split = line.split()
        keyword = line_split[0]
        keyword_opts = line_split[1:]
        # Handle keywords
        if keyword == "read_data":
            # This needs to be generalized in case a path is given. (?)
            file_path = keyword_opts[0]
            if not os.path.isabs(file_path):
                file_path = os.path.join(input_dir, file_path)
            read_lammps_data_file(dl, file_path, extra_simulation_data, blocksize)
        elif keyword == "include_data":
            include_data = eex.utility.read_lines(keyword_opts[0])
            for inum, line_data in enumerate(include_data):
                input_file.insert(lnum + inum + 1, line_data)
        elif keyword == "variable":
            tmp = {}

            variable_name = keyword_opts[0]
            variable_style = keyword_opts[1]
            variable_values = keyword_opts[2:]

            tmp["style"] = variable_style
            tmp["values"] = variable_values

            variable_list[variable_name] = tmp

        elif keyword in ["bond_style", "angle_style", "dihedral_style"]:
            extra_simulation_data[keyword] = keyword_opts[0]

        elif keyword == "units":
            extra_simulation_data["units"] = keyword_opts[0]

        elif keyword == "special_bonds":
            pass


def get_special_bonds():
    exclusions = {}
    special_bonds_keywords = re.findall("[a-zA-Za-z\/]+", " ".join(keyword_opts)) 
    for idx, exclusions_keyword in enumerate(special_bonds_keywords):
        if exclusions_keyword in lmd.exclusions.keys():
            if (exclusions_keyword == 'lj'):
                exclusions["lj"] = {}
                exclusions["lj"]["scale12"] = float(keyword_opts[idx+1])
                exclusions["lj"]["scale13"] = float(keyword_opts[idx+2])
                exclusions["lj"]["scale14"] = float(keyword_opts[idx+3])
            elif (exclusions_keyword == 'coul'):
                exclusions["coul"] = {}
                exclusions["coul"]["scale12"] = float(keyword_opts[idx+1])
                exclusions["coul"]["scale13"] = float(keyword_opts[idx+2])
                exclusions["coul"]["scale14"] = float(keyword_opts[idx+3])
            elif (exclusions_keyword == 'lj/coul'):
                exclusions["lj"] = {}
                exclusions["lj"]["scale12"] = float(keyword_opts[idx+1])
                exclusions["lj"]["scale13"] = float(keyword_opts[idx+2])
                exclusions["lj"]["scale14"] = float(keyword_opts[idx+3])
                exclusions["coul"] = {}
                exclusions["coul"]["scale12"] = exclusions["lj"]["scale12"]
                exclusions["coul"]["scale13"] = exclusions["lj"]["scale13"]
                exclusions["coul"]["scale14"] = exclusions["lj"]["scale14"]
            else: 
                exclusions["coul"] = lmd.exclusions[exclusion_keyword]["coul"]
                exclusions["lj"] = lmd.exclusions[exclusion_keyword]["lj"]
        elif exclusions_keyword in lmd.exclusions["additional_keywords"]:
            raise Exception("Keyword %s is not currently supported", exclusion_keyword)
    
    if "coul" not in exclusions:
        exclusions["coul"] = lmd.exclusions["default"]["coul"]
    if "lj" not in exclusions:
        exclusions["lj"] = lmd.exclusions["default"]["lj"]
       
    return exclusions
    # dl.set_exclusions(exclusions)
