"""
Writer for amber

"""

import time
import pandas as pd
import math
import re
import numpy as np
from collections import Counter

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

    data = data.ravel()

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


def _get_charmm_dihedral_count(dl):
    # TODO this a temporary solution to get the number of dihedrals stored in the datalayer.
    order = 4
    ret = 0
    charmm = amd.forcefield_parameters['dihedral']['form']
    terms = dl.list_term_parameters(order)

    for j in terms.values():
        term_md = eex.metadata.get_term_metadata(order, "forms", j[0])
        # Zip up the parameters
        parameters = {k: v for k, v in zip(term_md["parameters"], j[1:])}
        parameters = eex.form_converters.convert_form(order, parameters, j[0], charmm)

        if isinstance(parameters['K'], float):
            ret += 1
        elif isinstance(parameters['K'], np.ndarray):
            ret += len(parameters['K'])

    return ret


def _check_dl_compatibility(dl):
    """
    This function examines a datalayer to determine if it is compatible with Amber.

    Conversions between functional forms and pairwise interaction mixing are performed (if possible).
    """

    # Loop over force field information - check functional form compatibility
    for k, v in amd.forcefield_parameters.items():
        if k is not "nonbond":
            terms = dl.list_term_parameters(v["order"])

            for j in terms.values():
                term_md = eex.metadata.get_term_metadata(v["order"], "forms", j[0])
                canonical_form = term_md['canonical_form']
                compatible_forms = eex.metadata.get_term_metadata(v["order"], "group")[canonical_form]

                if v['form'] not in compatible_forms:
                    # Will need to insert check to see if these can be easily converted (ex OPLS dihedral <-> charmmfsw)
                    raise TypeError("Functional form %s stored in datalayer is not compatible with Amber.\n" % (j[0]))
        else:
            # Handle nonbonds. Amber must have pair interactions.

            # Grab all pair interactions stored in datalayer
            nb = dl.list_nb_parameters(nb_name="LJ", nb_model="AB", itype="pair")

            # Get the number of atom types. The number of pair interactions should be equal to num_atom_types * (num_atom_types + 1)) / 2
            num_atom_types = len(dl.get_unique_atom_types())

            # This will occur if there are no pair interactions stored in the dl.
            if len(nb.keys()) == 0:

                # Check that the stored mixing rule is compatible with amber. Amber should be able to handle any set of parameters,
                # but if we allow this, it should be something the user has to override somehow (or, they could apply the mixing
                # rule before calling the amber writer.
                if dl.get_mixing_rule() not in amd.mixing_rule:
                    raise TypeError(
                        "Mixing rule %s not compatible with amber (lorentz-berthelot mixing rule)" % dl.get_mixing_rule())

                # Calculate pair interactions according to amber.
                dl.build_LJ_mixing_table()

            # This condition will be met if some pair interactions are stored, but not the correct number.
            elif len(nb.keys()) != (num_atom_types * (num_atom_types + 1)) / 2:
                raise ValueError("Amber compatibility check : Incorrect number of pair interactions\n")

            # Check NB scaling factors are compatible with amber

            scaling_types = eex.metadata.additional_metadata.nb_scaling["scaling_type"]

            for scale_type in scaling_types:

                if scale_type not in dl.store.list_tables():

                    if not dl.get_nb_scaling_factors():
                        raise ValueError("No nonbond scaling (%s) information is set in datalayer" %(scale_type))
                    else:
                        # Build atom-wise scaling list
                        dl.build_scaling_list()

                # Check that NB scaling factors are compatible with amber (ie 1,2 and 1,3 must be 0 (excluded))
                pair_scalings = dl.get_pair_scalings(nb_labels=[scale_type], order=True)

                p12 = pair_scalings[pair_scalings["order"] == 2][scale_type]

                p13 = pair_scalings[pair_scalings["order"] == 3][scale_type]

                if p12.nonzero()[0]:
                    raise ValueError("Nonbond scaling (order=2, %s) is not consistent with Amber. In Amber, 1-2 nonbond "
                                     "interactions are excluded" %(scale_type))

                if p13.nonzero()[0]:
                    raise ValueError("Nonbond scaling (order=3, %s) is not consistent with Amber. In Amber, 1-3 nonbond "
                                     "interactions are excluded" %(scale_type))


    stored_properties = dl.list_atom_properties()
    required_properties = list(amd.atom_property_names.values())

    diff = np.setdiff1d(required_properties, stored_properties)
    natoms = dl.get_atom_count()

    index = np.arange(1, natoms + 1)

    # Build and curate the data
    df = pd.DataFrame({'atom_index': index})
    df.dropna(axis=0, how="any", inplace=True)
    df.set_index('atom_index', inplace=True)

    add_properties = []

    # Fill in default or raise error
    for req in diff:
        if req == 'atom_name':
            atom_names = ['A'] * natoms
            df[req] = atom_names
            add_properties.append(req)

        elif req == 'atomic_number':
            # Just say it's carbon...doesn't seem like this matters too much for amber
            atomic_numbers = [6] * natoms
            df[req] = atomic_numbers
            add_properties.append(req)
        elif req == "mass":
            try:
                dl.get_atoms(properties=["mass"])
            except:
                raise KeyError("No masses stored in datalayer")
        else:
            raise KeyError("Atom property %s is missing from datalayer" % (req))

    # Check for residue_index

    if "residue_index" not in stored_properties:
        # If molecule_index is set, set residue index to this.
        # Otherwise, set all to 1.0
        if "molecule_index" in stored_properties:
            df["residue_index"] = dl.get_atoms(properties=["molecule_index"])
            add_properties.append("residue_index")
        else:
            df["residue_index"] = 1
            add_properties.append("residue_index")

        if "residue_name" not in stored_properties:
            df["residue_name"] = "BLA"

    elif "residue_name" not in stored_properties:
        df["residue_name"] = "BLA"
        add_properties.append("residue_name")

    if len(add_properties) > 0:
        dl.add_atoms(df, by_value=True)


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

    # First get information into Amber pointers. All keys are initially filled with zero.
    # Ones that are currently 0, but should be implemented eventually are marked with

    _check_dl_compatibility(dl)

    dihedral_count = _get_charmm_dihedral_count(dl)

    # Figure out what is hydrogen for the header
    num_H_list = []
    inc_hydrogen = {}
    without_hydrogen = {}
    hidx = (dl.get_atoms("atomic_number") == 1)['atomic_number']
    hidx = hidx[hidx].index
    for term_type, term_name in zip([2, 3, 4], ["bonds", "angles", "dihedrals"]):
        term = dl.get_terms(term_type)
        if term.shape[0] == 0:
            num_H_list.append(0)
            continue

        # Build up an index of what is in hydrogen or not
        inc_hydrogen_mask = term["atom1"].isin(hidx)

        for n in range(term_type - 1):
            name = "atom" + str(n + 2)
            inc_hydrogen_mask |= term[name].isin(hidx)

        num_H_list.append(len(term.loc[inc_hydrogen_mask].values))
        inc_hydrogen[term_name] = term.loc[inc_hydrogen_mask].values
        without_hydrogen[term_name] = term.loc[~inc_hydrogen_mask].values

    output_sizes = {k: 0 for k in amd.size_keys}

    output_sizes['NATOM'] = dl.get_atom_count()  # Number of atoms
    output_sizes["NBONH"] = num_H_list[0]  # Number of bonds containing hydrogen
    output_sizes["MBONA"] = dl.get_term_count(2, "total") - output_sizes["NBONH"]  # Number of bonds not containing hydrogen
    output_sizes['NBONA'] = output_sizes["MBONA"]  # MBONA + number of constraint bonds (MBONA = NBONA always)
    output_sizes["NTHETH"] = num_H_list[1]  # Number of angles containing hydrogen
    output_sizes["MTHETA"] = dl.get_term_count(3, "total") - output_sizes["NTHETH"]  # Number of angles not containing hydrogen
    output_sizes['NTHETA'] = output_sizes["MTHETA"]  # MTHETA + number of constraint angles (NTHETA = MTHETA always)
    output_sizes["NPHIH"] = num_H_list[2]  # Number of torsions containing hydrogen
    output_sizes["MPHIA"] = dl.get_term_count(4, "total") - output_sizes["NPHIH"]  # Number of torsions not containing hydrogen
    output_sizes["NPHIA"] = output_sizes["MPHIA"]
    output_sizes["NUMBND"] = len(dl.list_term_uids(2))  # Number of unique bond types
    output_sizes["NUMANG"] = len(dl.list_term_uids(3))  # Number of unique angle types
    output_sizes["NPTRA"] = dihedral_count  # Number of unique torsion types

    output_sizes["NRES"] = len(dl.list_atom_uids("residue_name"))  # Number of residues (not stable)
    output_sizes["NTYPES"] = len(np.unique(dl.get_atoms("atom_type")))  # Number of distinct LJ atom types

    output_sizes["NPARM"] = 0  # Used to determine if this is a LES-compatible prmtop (??)
    output_sizes["NNB"] = dl.get_atom_count(
    )  # Number of excluded atoms - Set to num atoms for our test cases. Amber will not run with 0
    # 0 - no box, 1 - orthorhombic box, 2 - truncated octahedron
    output_sizes["NMXRS"] = 0  # Number of atoms in the largest residue
    output_sizes["IFCAP"] = 0  # Set to 1 if a solvent CAP is being used
    output_sizes["NUMEXTRA"] = 0  # Number of extra points in the topology file

    # Needs check for orthorhomibic box (1) or truncated octahedron (2). Currently just 0 or 1
    output_sizes["IFBOX"] = [0 if dl.get_box_size() == {} else 1][0]  # Flag indicating whether a periodic box is present

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

    # Write title and version information
    f = open(filename, "w")
    f.write('%%VERSION  VERSION_STAMP = V0001.000  DATE = %s  %s\n' % (time.strftime("%x"), time.strftime("%H:%M:%S")))
    f.write("%FLAG TITLE\n%FORMAT(20a4)\n")
    f.write("prmtop generated by MolSSI EEX\n")

    # Write pointers section
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

    # Write atom properties sections
    file_handle = open(filename, "ab")

    for k in amd.atom_property_names:
        # Get data
        data = dl.get_atoms(amd.atom_property_names[k], by_value=True, utype=amd.atom_data_units).values.ravel()
        _write_amber_data(file_handle, data, k)

        written_categories.append(k)

    # Handle residues

    # We assume these are sorted WRT to atom and itself at the moment... not great
    res_data = dl.get_atoms(["residue_index", "residue_name"], by_value=True)
    uvals, uidx, ucnts = np.unique(res_data["residue_index"], return_index=True, return_counts=True)

    labels = res_data["residue_name"].iloc[uidx].values
    _write_amber_data(file_handle, labels, "RESIDUE_LABEL")
    written_categories.append("RESIDUE_LABEL")

    starts = np.concatenate(([1], np.cumsum(ucnts) + 1))[:-1]
    _write_amber_data(file_handle, starts, "RESIDUE_POINTER")
    written_categories.append("RESIDUE_POINTER")

    # Write out term parameters
    for term_type in ["bond", "angle", "dihedral"]:
        uids = sorted(dl.list_term_uids(term_type))

        if len(uids) == 0:
            continue
        term_md = amd.forcefield_parameters[term_type]

        tmps = {k: [] for k in term_md["column_names"].keys()}
        utype = term_md["units"]
        order = term_md["order"]

        inv_lookup = {v: k for k, v in term_md["column_names"].items()}

        # Build lists of data since AMBER holds this as 1D
        for uid in uids:
            params = dl.get_term_parameter(order, uid, utype=utype, ftype=term_md['form'])

            for k, v in params[1].items():
                tmps[inv_lookup[k]].append(v)
        # Write out FLAGS

        for k, v in tmps.items():

            _write_amber_data(file_handle, v, k)
            written_categories.append(k)

    for term_type, term_name in zip([2, 3, 4], ["bonds", "angles", "dihedrals"]):
        term = dl.get_terms(term_type)

        if term.shape[0] == 0:
            continue

        # Build up an index of what is in hydrogen or not
        inc_hydrogen_mask = term["atom1"].isin(hidx)

        # Scale by weird AMBER factors
        inc_hydrogen[term_name][:, :-1] = (inc_hydrogen[term_name][:, :-1] - 1) * 3
        without_hydrogen[term_name][:, :-1] = (without_hydrogen[term_name][:, :-1] - 1) * 3

        inc_h_name = term_name.upper() + "_INC_HYDROGEN"
        without_h_name = term_name.upper() + "_WITHOUT_HYDROGEN"

        _write_amber_data(file_handle, inc_hydrogen[term_name], inc_h_name)
        written_categories.append(inc_h_name)

        _write_amber_data(file_handle, without_hydrogen[term_name], without_h_name)
        written_categories.append(without_h_name)

    # Handle SOLVENT_POINTERS, ATOMS_PER_MOLECULE and BOX_DIMENSIONS. Only present if IFBOX>0.
    if output_sizes["IFBOX"] > 0:
        # Solvent pointers section
        # There are three numbers here - IPTRES, NSPM, NSPSOL
        # where
        # IPTRES = final residue part of solute, NSPM = total number of molecules, NSPSOL = first solvent molecule
        # Just say everything is solute for now.

        iptres = dl.get_atoms(["residue_index"]).values[-1]
        nspm = len(np.unique(dl.get_atoms(["molecule_index"]).values))

        solvent_pointers = [iptres, nspm, nspm]

        _write_amber_data(file_handle, solvent_pointers, "SOLVENT_POINTERS")

        # Handle atoms per molecule
        molecule_list = dl.get_atoms(["molecule_index"]).values.ravel()
        count_atoms_per_molecule = Counter(molecule_list)
        atoms_per_molecule = []

        for x in range(1, nspm + 1):
            atoms_per_molecule.append(count_atoms_per_molecule[x])

        _write_amber_data(file_handle, atoms_per_molecule, "ATOMS_PER_MOLECULE")

        # Write box dimensions section
        box_dimensions = dl.get_box_size(utype={"a": amd.box_units["length"], "b": amd.box_units["length"],
                                                "c": amd.box_units["length"], "alpha": amd.box_units["angle"],
                                                "beta": amd.box_units["angle"], "gamma": amd.box_units["angle"]})

        write_box = [box_dimensions["beta"], box_dimensions["a"], box_dimensions["b"], box_dimensions["c"]]

        _write_amber_data(file_handle, write_box, "BOX_DIMENSIONS")

    written_categories.append("BOX_DIMENSIONS")
    written_categories.append("SOLVENT_POINTERS")
    written_categories.append("ATOMS_PER_MOLECULE")

    # Quick fix for radius set  will be one line string description in files prepared by xleap
    _write_amber_data(file_handle, ["Place holder - EEX"], "RADIUS_SET")
    written_categories.append("RADIUS_SET")

    # Handle NB data
    # Relevant headers = NONBOND_PARM_INDEX, LENNARD_JONES_ACOEF, LENNARD_JONES_BCOEF
    stored_atom_types = dl.get_unique_atom_types()
    ntypes = len(stored_atom_types)

    nb_forms = dl.list_stored_nb_types()

    # This can be removed if compatibility check inserted at beginning
    if set(nb_forms) != set(["LJ"]):
        # Write better message here
        raise KeyError("Nonbond forms stored in datalayer are not compatible with Amber - %s" % nb_forms)

    # Get parameters from datalayer using correct amber units
    stored_nb_parameters = dl.list_nb_parameters(
        nb_name="LJ", nb_model="AB", utype=amd.forcefield_parameters["nonbond"]["units"], itype="pair")
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

    # Write headers for other sections (file will not work in AMBER without these)
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
