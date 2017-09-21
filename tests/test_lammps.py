"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import eex_find_files


def test_lammps_read():

    fname = eex_find_files.get_example_filename("lammps", "SPCE/data.spce")

    dl = eex.datalayer.DataLayer("test_lammps_read")

    data = eex.translators.lammps.read_lammps_file(dl, fname)

    # Check on the data dictionary
    assert data["sizes"]["atoms"] == 600
    assert data["sizes"]["bonds"] == 400
    assert data["sizes"]["angles"] == 200
    assert data["sizes"]["angle types"] == 1

    assert data["dimensions"]["xlo"] == -12.362
    assert data["dimensions"]["xhi"] == 12.362

    # Check Atoms
    atoms = dl.get_atoms(["atom_type", "charge"])
    assert atoms.shape[0] == 600
    assert np.allclose(np.unique(atoms["atom_type"]), [1, 2])
    assert np.allclose(np.unique(atoms["charge"]), [-0.8476,  0.4238])

    # Check Bonds
    bonds = dl.get_bonds()
    assert bonds.shape[0] == 400
    assert np.allclose(np.unique(bonds["bond_type"]), [1])

    # Check Angles
    angles = dl.get_angles()
    assert angles.shape[0] == 200
    assert np.allclose(np.unique(angles["angle_type"]), [1])

