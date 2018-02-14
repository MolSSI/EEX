"""
Writer for amber

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


def _write_1d(file_handle, data, ncols, fmt):

    remainder_size = data.size % ncols
    if data.size == 0:
        file_handle.write("\n".encode())
    elif remainder_size == 0:
        np.savetxt(file_handle, data.reshape(-1, ncols), fmt=fmt, delimiter="")
    else:
        rem_data = data[-remainder_size:].reshape(1, -1)
        data = data[:-remainder_size].reshape(-1, ncols)
        np.savetxt(file_handle, data, fmt=fmt, delimiter="")
        np.savetxt(file_handle, rem_data, fmt=fmt, delimiter="")
    # print(data.shape, rem_data.shape)

    # Write data to file
    file_handle.flush()


def _write_amber_data(file_handle, data, category):
    fmt_string = amd.data_labels[category][1]
    fmt_data = amd.parse_format(fmt_string)

    file_handle.write(("%%FLAG %s\n" % category).encode())
    file_handle.write((fmt_string + "\n").encode())

    ncols = fmt_data[0]
    fmt = amd.build_format(fmt_data)

    _write_1d(file_handle, np.array(data), ncols, fmt)

def _check_dl_compatibility(dl):
    """
    This function examines a datalayer to determine if it is compatible with Amber.

    Conversions between functional forms and pairwise interaction mixing are performed (if possible).
    """
    print("Checking dl compatibility\n")

    # Loop over force field information. 
    for k, v in amd.forcefield_parameters.items():
        if k != "nonbond":
            terms = dl.list_term_parameters(v["order"])
            for j in terms.values():
                if j[0] != v["form"]:
                    # Will need to insert check to see if these can be easily converted (ex OPLS dihedral <-> charmmfsw)
                    raise Exception("Functional form stored in datalayer is not compatible with Amber")
        else:
            pass




def write_amber_file(dl, filename, inpcrd=None):
    """
    Parameters
    ------------
    dl : eex.DataLayer
        The datalayer containing information about the system to write
    filename : str
        The name of the file to write
    inpcrd : str, optional
        If None, attempts to read the file filename.replace("prmtop", "inpcrd") otherwise passes. #
    """

    ### First get information into Amber pointers. All keys are initially filled with zero.
    # Ones that are currently 0, but should be implemented eventually are marked with

    _check_dl_compatibility(dl)

    output_sizes = {k: 0 for k in amd.size_keys}

    output_sizes['NATOM'] = dl.get_atom_count()  # Number of atoms
    output_sizes["MBONA"] = dl.get_term_count(2, "total")  #  Number of bonds not containing hydrogen
    output_sizes['NBONA'] = output_sizes["MBONA"]  # MBONA + number of constraint bonds (MBONA = NBONA always)
    output_sizes["MTHETA"] = dl.get_term_count(3, "total")  #  Number of angles not containing hydrogen
    output_sizes['NTHETA'] = output_sizes["MTHETA"]  # MTHETA + number of constraint angles (NTHETA = MTHETA always)
    output_sizes["MPHIA"] = dl.get_term_count(4, "total")  #  Number of torsions not containing hydrogen
    output_sizes["NPHIA"] = output_sizes["MPHIA"]
    output_sizes["NUMBND"] = len(dl.list_term_uids(2))  # Number of unique bond types
    output_sizes["NUMANG"] = len(dl.list_term_uids(3))  # Number of unique angle types
    output_sizes["NPTRA"] = len(dl.list_term_uids(4))  # Number of unique torsion types
    output_sizes["NRES"] = len(dl.list_atom_uids("residue_name"))  # Number of residues (not stable)
    output_sizes["NTYPES"] = len(np.unique(dl.get_atoms("atom_type")))  # Number of distinct LJ atom types
    output_sizes["NBONH"] = 0  #  Number of bonds containing hydrogen
    output_sizes["NTHETH"] = 0  #  Number of angles containing hydrogen
    output_sizes["NPHIH"] = 0  #  Number of torsions containing hydrogen
    output_sizes["NPARM"] = 0  #  Used to determine if this is a LES-compatible prmtop (??)
    output_sizes["NNB"] = dl.get_atom_count(
    )  #  Number of excluded atoms - Set to num atoms for our test cases. Amber will not run with 0
    # 0 - no box, 1 - orthorhombic box, 2 - truncated octahedron
    output_sizes["NMXRS"] = 0  #  Number of atoms in the largest residue
    output_sizes["IFCAP"] = 0  # Set to 1 if a solvent CAP is being used
    output_sizes["NUMEXTRA"] = 0  # Number of extra points in the topology file

    ## Needs check for orthorhomibic box (1) or truncated octahedron (2). Currently just 0 or 1
    output_sizes["IFBOX"] = [0 if dl.get_box_size == {} else 1][0]  # Flag indicating whether a periodic box is present

    written_categories = []

    # Figure out size each section should be based on metadata
    label_sizes = {}
    for k, v in amd.data_labels.items():
        if isinstance(v[0], int):
            label_sizes[k] = v[0]
        elif v[0] in list(output_sizes):
            label_sizes[k] = output_sizes[v[0]]
        else:
            # print("%30s %40s %d" % (k, v[0], int(eval(v[0], sizes_dict))))
            label_sizes[k] = int(eval(v[0], output_sizes))

    ### Write title and version information
    f = open(filename, "w")
    f.write('%%VERSION  VERSION_STAMP = V0001.000  DATE = %s  %s\n' % (time.strftime("%x"), time.strftime("%H:%M:%S")))
    f.write("%FLAG TITLE\n%FORMAT(20a4)\n")
    f.write("prmtop generated by MolSSI EEX\n")

    ## Write pointers section
    f.write("%%FLAG POINTERS\n%s\n" % (amd.data_labels["POINTERS"][1]))
    ncols, dtype, width = amd.parse_format(amd.data_labels["POINTERS"][1])
    format_string = "%%%sd" % width

    count = 0
    for k in amd.size_keys:
        f.write(format_string % output_sizes[k])
        count += 1
        if count % ncols == 0:
            f.write("\n")

    f.write("\n")
    f.close()
    written_categories.append("POINTERS")

    ### Write atom properties sections
    file_handle = open(filename, "ab")

    for k in amd.atom_property_names:

        # Get unit type
        utype = None
        if k in amd.atom_data_units:
            utype = amd.atom_data_units[k]

        # Get data
        data = dl.get_atoms(amd.atom_property_names[k], by_value=True, utype=utype).values.ravel()
        _write_amber_data(file_handle, data, k)

        written_categories.append(k)

    ### Handle residues

    # We assume these are sorted WRT to atom and itself at the moment... not great
    res_data = dl.get_atoms(["residue_index", "residue_name"], by_value=True)
    uvals, uidx, ucnts = np.unique(res_data["residue_index"], return_index=True, return_counts=True)

    labels = res_data["residue_name"].iloc[uidx].values
    _write_amber_data(file_handle, labels, "RESIDUE_LABEL")
    written_categories.append("RESIDUE_LABEL")

    starts = np.concatenate(([1], np.cumsum(ucnts) + 1))[:-1]
    _write_amber_data(file_handle, starts, "RESIDUE_POINTER")
    written_categories.append("RESIDUE_POINTER")

    ### Write out term parameters
    for term_type in ["bond", "angle", "dihedral"]:
        uids = sorted(dl.list_term_uids(term_type))

        if len(uids) == 0: continue
        term_md = amd.forcefield_parameters[term_type]

        tmps = {k: [] for k in term_md["column_names"].keys()}
        utype = term_md["units"]
        order = term_md["order"]
        inv_lookup = {v: k for k, v in term_md["column_names"].items()}

        # Build lists of data since AMBER holds this as 1D
        for uid in uids:
            params = dl.get_term_parameter(order, uid, utype=utype)
            for k, v in params[1].items():
                tmps[inv_lookup[k]].append(v)

        # Write out FLAGS
        for k, v in tmps.items():

            _write_amber_data(file_handle, v, k)
            written_categories.append(k)

    ### Handle term data
    hidx = (dl.get_atoms("atomic_number") == 1).values.ravel()

    for term_type, term_name in zip([2, 3, 4], ["bonds", "angles", "dihedrals"]):
        term = dl.get_terms(term_type)

        if term.shape[0] == 0: continue

        # Build up an index of what is in hydrogen or not
        inc_hydrogen_mask = term["atom1"].isin(hidx)
        for n in range(term_type - 1):
            name = "atom" + str(n + 2)
            inc_hydrogen_mask |= term[name].isin(hidx)

        inc_hydrogen = term.loc[inc_hydrogen_mask].values
        without_hydrogen = term.loc[~inc_hydrogen_mask].values

        # Scale by weird AMBER factors
        inc_hydrogen[:, :-1] = (inc_hydrogen[:, :-1] - 1) * 3
        without_hydrogen[:, :-1] = (without_hydrogen[:, :-1] - 1) * 3

        inc_h_name = term_name.upper() + "_INC_HYDROGEN"
        without_h_name = term_name.upper() + "_WITHOUT_HYDROGEN"

        _write_amber_data(file_handle, inc_hydrogen, inc_h_name)
        written_categories.append(inc_h_name)

        _write_amber_data(file_handle, without_hydrogen, without_h_name)
        written_categories.append(without_h_name)

    # Append & forget about SOLVENT_POINTERS section for now
    written_categories.append("SOLVENT_POINTERS")

    # Write box dimensions section (if applicable)
    if output_sizes["IFBOX"] > 0:
        box_dimensions = dl.get_box_size(utype={"a": amd.box_units["length"], "b": amd.box_units["length"],
                                                     "c" : amd.box_units["length"], "alpha": amd.box_units["angle"],
                                                     "beta": amd.box_units["angle"], "gamma": amd.box_units["angle"]})

        write_box = [box_dimensions["beta"], box_dimensions["a"], box_dimensions["b"], box_dimensions["c"]]

        _write_amber_data(file_handle, write_box, "BOX_DIMENSIONS")

    written_categories.append("BOX_DIMENSIONS")

    # Quick fix for radius set  will be one line string description in files prepared by xleap
    _write_amber_data(file_handle, ["Place holder - EEX"], "RADIUS_SET")
    written_categories.append("RADIUS_SET")

    # Handle NB data
    # Relevant headers = NONBOND_PARM_INDEX, LENNARD_JONES_ACOEF, LENNARD_JONES_BCOEF
    stored_atom_types = dl.get_unique_atom_types()
    ntypes = len(stored_atom_types)

    nb_forms = dl.list_stored_nb_types()

    if set(nb_forms) != set(["LJ"]):
        # Write better message here
        raise KeyError("Nonbond forms stored in datalayer are not compatible with Amber - %s" % nb_forms)

    # Get parameters from datalayer using correct amber units
    stored_nb_parameters = dl.list_nb_parameters(
        nb_name="LJ", nb_model="AB", utype=amd.forcefield_parameters["nonbond"]["units"])
    nonbonded_parm_index = np.zeros(ntypes * ntypes)
    lj_a_coeff = []
    lj_b_coeff = []

    # Build a_coeff, b_coeff, and nb_parm lists
    for key, value in stored_nb_parameters.items():
        lj_a_coeff.append(value['A'])
        lj_b_coeff.append(value['B'])
        index_to_nb = ntypes * (key[0] - 1) + key[1]
        index_to_nb2 = ntypes * (key[1] - 1) + key[0]
        nonbonded_parm_index[index_to_nb - 1] = len(lj_a_coeff)
        nonbonded_parm_index[index_to_nb2 - 1] = len(lj_a_coeff)

    _write_amber_data(file_handle, nonbonded_parm_index, "NONBONDED_PARM_INDEX")
    _write_amber_data(file_handle, lj_a_coeff, "LENNARD_JONES_ACOEF")
    _write_amber_data(file_handle, lj_b_coeff, "LENNARD_JONES_BCOEF")

    for category in amd.forcefield_parameters["nonbond"]["column_names"]:
        written_categories.append(category)

    ### Write headers for other sections (file will not work in AMBER without these)
    for k in amd.data_labels:
        if k not in written_categories:
            if label_sizes[k] > 0:
                data = np.zeros(label_sizes[k])
                _write_amber_data(file_handle, data, k)
            else:
                file_handle.write(("%%FLAG %s\n%s\n\n" % (k, amd.data_labels[k][1])).encode())
            written_categories.append(k)

    file_handle.close()

    # Now we need to write out the INPCRD
    if '.prmtop' in filename:
        inpcrd_file = filename.replace('.prmtop', '.inpcrd')
    else:
        inpcrd_file = filename + '.inpcrd'

    file_handle = open(inpcrd_file, "wb")

    xyz = dl.get_atoms("XYZ", utype={"XYZ": "angstrom"})

    file_handle.write("default_name\n".encode())
    file_handle.write(("%6d\n" % xyz.shape[0]).encode())

    _write_1d(file_handle, xyz.values.ravel(), 6, "%12.6f")

    if output_sizes["IFBOX"] > 0:
        box = pd.DataFrame(box_dimensions, index=[0])
        box = box[['a', 'b', 'c', 'alpha', 'beta', 'gamma']]
        _write_1d(file_handle, box.values.ravel(), 6, "%12.6f")



    file_handle.close()

    return 0
