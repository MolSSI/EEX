"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import eex_find_files


def test_amber_read():

    fname = eex_find_files.get_example_filename("amber", "water/spce.prmtop")

    dl = eex.datalayer.DataLayer("test_amber_read")

    data = eex.translators.amber.read_amber_file(dl, fname)

    assert data["VERSION"] == "V0001.000"
    #    # Check on the data dictionary
    #    assert data["sizes"]["atoms"] == 600
    #    assert data["sizes"]["bonds"] == 400
    #    assert data["sizes"]["angles"] == 200
    #    assert data["sizes"]["angle types"] == 1
    #
    #    assert data["dimensions"]["xlo"] == -12.362
    #    assert data["dimensions"]["xhi"] == 12.362
    #
    #    # Check Atoms
    atoms = dl.get_atoms(["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index"])
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == set(["H1", "H2", "O"])
    assert np.allclose(np.unique(atoms["charge"]), [-1.54452215E+01, 7.72261074E+00])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [1.008, 16.0])
    assert set(np.unique(atoms["residue_name"])) == set(["WAT"])
    assert np.allclose(np.min(atoms["residue_index"]), 0)
    assert np.allclose(np.max(atoms["residue_index"]), 215)
    # print(atoms["residue_index"])

    # assert False
#
#    # Check Bonds
#    bonds = dl.get_bonds()
#    assert bonds.shape[0] == 400
#    assert np.allclose(np.unique(bonds["bond_type"]), [1])
#
#    # Check Angles
#    angles = dl.get_angles()
#    assert angles.shape[0] == 200
#    assert np.allclose(np.unique(angles["angle_type"]), [1])
