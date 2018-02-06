"""
LAMMPS EEX I/O
"""

import pandas as pd
import math
import numpy as np
import eex

from . import lammps_metadata as lmd

import logging
logger = logging.getLogger(__name__)


def write_lammps_file(dl, filename, utype="real", blocksize=110):

    # handle units
    unit_set = lmd.units_style[utype]


    term_table = lmd.build_term_table(utype)

    nb_term_table = lmd.build_nb_table("real")

    data_file = open(filename, 'w')
    input_file = open(filename+".in", 'w')

    data_file.write("LAMMPS data file generated by MolSSI EEX\n\n")
    input_file.write("# LAMMPS input file generated by MolSSI EEX\n")

    input_file.write("units\t%s\n" %(utype))
    input_file.write("atom_style\tfull\n")
    input_file.write("read_data\t%s\n" %(filename))

    sizes = {}
    sizes["atoms"] = dl.get_atom_count()
    sizes["bonds"] = dl.get_term_count(2, "total")
    sizes["angles"] = dl.get_term_count(3, "total")
    sizes["dihedrals"] = dl.get_term_count(4, "total") # Not qutie right once we do impropers
    sizes["impropers"] = 0
    sizes["atom types"] = len(dl.list_atom_uids("mass"))

    # All the UID's minus the "total" columns
    sizes["bond types"] = len(dl.get_term_count(2)) - 1
    sizes["angle types"] = len(dl.get_term_count(3)) - 1
    sizes["dihedral types"] = len(dl.get_term_count(4)) - 1
    sizes["improper types"] = 0

    # Write header information
    for k in lmd.size_keys:
        data_file.write(" %d %s\n" % (sizes[k], k))
        # data_file.write(' '.join([str(data["sizes"][k]), k, '\n']))

    # Write box information
    box_size = dl.get_box_size(utype={"a": unit_set["[length]"], "b": unit_set["[length]"], "c": unit_set["[length]"],
                                      "alpha": "degree", "beta": "degree", "gamma": "degree"})

    box_center = dl.get_box_center(utype={"x": "angstrom", "y": "angstrom", "z": "angstrom"})

    lo_hi = {}

    lo_hi["x"] = [box_center["x"] - box_size["a"]/2., box_center["x"] + box_size["a"]/2.]
    lo_hi["y"] = [box_center["y"] - box_size["b"]/2., box_center["y"] + box_size["b"] / 2.]
    lo_hi["z"] = [box_center["z"] - box_size["c"] / 2., box_center["z"] + box_size["c"] / 2.]

    for k, v in lo_hi.items():
        data_file.write("% 8.6f% 8.6f %slo %shi\n" % (v[0],v[1], k, k))
    data_file.write('\n')


    data_file.write(("Pair Coeffs\n\n").title())

    param_fmt = "%10.8f"
    # Loop over Pair Coeffs
    nb_forms = dl.list_stored_nb_types()
    for form in nb_forms:
        if not form in nb_term_table:
            raise KeyError("Nonbond forms stored in datalayer are not compatible with Lammps forms - %s" % nb_forms)
        if form == "LJ": 
            stored_nb_parameters = dl.list_nb_parameters(
                    nb_name=form, nb_model="epsilon/sigma", utype={"epsilon": "kcal/mol", "sigma": "angstrom"})

            for key, value in stored_nb_parameters.items():
                data_file.write(("%2d %10.8f %10.8f\n" % (key[0], value['epsilon'], value['sigma'])))
    data_file.write("\n")
        
    # Loop over all of the parameter data
    for param_order, param_type in zip([2, 3, 4], ["bond", "angle", "dihedral"]):
        param_uids = dl.list_term_uids(param_order)

        if len(param_uids) == 0: continue

        data_file.write(("%s Coeffs\n\n" % param_type).title())
        for uid in param_uids:
            param_coeffs = dl.get_term_parameter(param_order, uid)


            term_data = term_table[param_order][param_coeffs[0]]
            param_coeffs = dl.get_term_parameter(param_order, uid, utype=term_data["utype"])

            # Order the data like lammps wants it
            parameters = [param_coeffs[1][k] for k in term_data["parameters"]]

            data_file.write("%2d " % uid)
            data_file.write(" ".join(param_fmt % f for f in parameters))
            data_file.write("\n")

            input_file.write("%s_style\t%s\n" %(param_type, param_coeffs[0]))

        data_file.write("\n")

    # Write out mass data
    data_file.write(" Masses\n\n")
    for idx in dl.list_atom_uids("mass"):
        mass = dl.get_atom_parameter("mass", idx, utype=lmd.get_context(utype, "[mass]"))
        data_file.write("%2d %10.8f\n" % (idx, mass))
    data_file.write('\n')

    # Write out atom data
    data_file.write(" Atoms\n\n")

    atoms = dl.get_atoms(["molecule_index", "atom_type", "charge", "xyz"], by_value=True)
    atoms.index = pd.RangeIndex(start=1, stop=atoms.shape[0] + 1)

    # Build a simple formatter
    def float_fmt(n):
        return "%10.8f" % n
    atoms.to_string(data_file, header=None, float_format=float_fmt)
    data_file.write('\n\n')

    # Write out all of the term data
    for param_order, param_type in zip([2, 3, 4], ["bonds", "angles", "dihedrals"]):
        if sizes[param_type] == 0: continue

        data_file.write((" %s\n\n" % param_type).title())

        # Grab term and reorder
        cols = ["term_index"] + ["atom%s" % d for d in range(1, param_order + 1)]
        term = dl.get_terms(param_type)[cols]
        term.index = pd.RangeIndex(start=1, stop=term.shape[0] + 1)
        # print(term)
        term.to_csv(data_file, header=None, sep=" ")
        data_file.write('\n')

    data_file.close()
    input_file.close()
    return True

