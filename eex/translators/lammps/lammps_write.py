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


def write_lammps_file(dl, data_filename, input_filename, unit_style="real", blocksize=110):

    # handle units
    unit_set = lmd.units_style[unit_style]


    term_table = lmd.build_term_table(unit_style)

    nb_term_table = lmd.build_nb_table("real")

    data_file = open(data_filename, 'w')
    input_file = open(input_filename, 'w')

    data_file.write("LAMMPS data file generated by MolSSI EEX\n\n")
    input_file.write("# LAMMPS input file generated by MolSSI EEX\n")

    input_file.write("units\t%s\n" %(unit_style))
    input_file.write("atom_style\tfull\n")

    sizes = {}
    sizes["atoms"] = dl.get_atom_count()
    sizes["bonds"] = dl.get_term_count(2, "total")
    sizes["angles"] = dl.get_term_count(3, "total")
    sizes["dihedrals"] = dl.get_term_count(4, "total") # Not qutie right once we do impropers
    sizes["impropers"] = 0
    sizes["atom types"] = len(dl.get_unique_atom_types())

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
    if box_center and box_size:
        lo_hi["x"] = [box_center["x"] - box_size["a"]/2., box_center["x"] + box_size["a"]/2.]
        lo_hi["y"] = [box_center["y"] - box_size["b"]/2., box_center["y"] + box_size["b"] / 2.]
        lo_hi["z"] = [box_center["z"] - box_size["c"] / 2., box_center["z"] + box_size["c"] / 2.]

        for k, v in lo_hi.items():
            data_file.write("% 8.6f% 8.6f %slo %shi\n" % (v[0],v[1], k, k))
        data_file.write('\n')

    param_fmt = "%10.8f"
    # Loop over Pair Coeffs
    nb_forms = dl.list_stored_nb_types()

    # Handle nonbonds - This needs to be generalized badly - make it work for now.
    if dl.get_mixing_rule() != '':
        data_file.write(("Pair Coeffs\n\n").title())

        stored_nb_parameters = dl.list_nb_parameters(
            nb_name="LJ", nb_model="epsilon/sigma", utype={"epsilon": unit_set["[energy]"], "sigma": unit_set["[length]"]}, itype="single")
    else:
        data_file.write("PairIJ Coeffs\n\n")
        stored_nb_parameters = dl.list_nb_parameters(
            nb_name="LJ", nb_model="epsilon/sigma", utype={"epsilon": unit_set["[energy]"], "sigma": unit_set["[length]"]}, itype="pair")

    for key, value in stored_nb_parameters.items():
        if key[1] == None:
            data_file.write(("%2d %10.8f %10.8f\n" % (key[0], value['epsilon'], value['sigma'])))
        else:
            data_file.write(("%2d %2d %10.8f %10.8f\n" % (key[0], key[1], value['epsilon'], value['sigma'])))

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

    # Write out mass data - don't use get_atom_parameter since we cannot assume that uid = atom_type
    data_file.write(" Masses\n\n")
    data = dl.get_atoms(["atom_type", "mass"],
                        utype={'atom_type': None, 'mass':lmd.get_context(unit_style, "[mass]")},
                        by_value=True).drop_duplicates()

    data.to_csv(data_file, sep=' ', index=False, header=False)

    # This is not working in Python 3?
    #np.savetxt(data_file, np.array(data), fmt='%2d %10.8f')

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

    input_file.write("read_data\t%s\n" %(data_filename))

    data_file.close()
    input_file.close()
    return True

