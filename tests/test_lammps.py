"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import os
import pytest
import pandas as pd
import eex_find_files


@pytest.fixture(scope="module", params=["HDF5", "Memory"])
def spce_dl(request):
    fname = eex_find_files.get_example_filename("lammps", "SPCE", "data.spce")
    dl = eex.datalayer.DataLayer(
        "test_lammps_read", )
    sim_data = {'units': 'real', 'bond_style': 'harmonic', 'angle_style': 'harmonic', 'dihedral_style': 'opls',
                'atom_style': 'full'}
    eex.translators.lammps.read_lammps_data_file(dl, fname, sim_data, blocksize=55)
    yield dl
    dl.close()


def test_lammps_read_data(spce_dl):
    dl = spce_dl

    # Check on the data dictionary
    utype = {"epsilon": "kcal / mol", "sigma": "angstrom"}
    nb_param_atom1 = dl.get_nb_parameter(atom_type=1, nb_model='epsilon/sigma', utype=utype)
    nb_param_atom2 = dl.get_nb_parameter(atom_type=2, nb_model='epsilon/sigma', utype=utype)
    assert np.allclose([nb_param_atom1['epsilon'], nb_param_atom1['sigma']], [0.15524976551, 3.166])
    assert np.allclose([nb_param_atom2['epsilon'], nb_param_atom2['sigma']], [0.0, 0.0])
    assert len(dl.get_unique_atom_types()) == 2
    assert len(dl.list_term_uids()[2]) == 1
    assert len(dl.list_term_uids()[3]) == 1
    assert len(dl.list_term_uids()[4]) == 0
    assert dl.get_atom_count() == 600
    assert dl.get_bond_count() == 400
    assert dl.get_angle_count() == 200
    assert dl.get_dihedral_count() == 0

    #box_size = dl.get_box_size()
    #assert box_size["x"] == pytest.approx(-12.362, 1.e-6)


def test_lammps_read_atoms(spce_dl):
    dl = spce_dl

    # Check Atoms
    atoms = dl.get_atoms(["atom_type", "charge", "mass"])
    assert atoms.shape[0] == 600
    assert np.allclose(np.unique(atoms["atom_type"]), [1, 2])
    assert np.allclose(np.unique(atoms["charge"]), [0, 1])
    assert np.allclose(np.unique(atoms["mass"]), [1, 2])


def test_lammps_read_atoms_value(spce_dl):
    dl = spce_dl

    # Check Atoms
    atoms = dl.get_atoms(["atom_type", "charge", "mass"], by_value=True)
    assert atoms.shape[0] == 600
    assert np.allclose(np.unique(atoms["atom_type"]), [1, 2])
    assert np.allclose(np.unique(atoms["charge"]), [-0.8476, 0.4238])
    assert np.allclose(np.unique(atoms["mass"]), [1.008, 16.000])


def test_lammps_read_bonds(spce_dl):
    dl = spce_dl

    # Check Bonds
    bonds = dl.get_bonds()
    assert bonds.shape[0] == 400
    assert np.allclose(np.unique(bonds["term_index"]), [1])


def test_lammps_read_angles(spce_dl):
    dl = spce_dl

    # Check Angles
    angles = dl.get_angles()
    assert angles.shape[0] == 200
    assert np.allclose(np.unique(angles["term_index"]), [1])


@pytest.mark.parametrize("molecule", [
    "data.trappe_butane_single_molecule",
    "data.trappe_propane_single_molecule",
    "data.trappe_ethane_single_molecule",
])
def test_lammps_writer(molecule):
    sim_data = {'units': 'real', 'bond_style': 'harmonic', 'angle_style': 'harmonic', 'dihedral_style': 'opls',
                'atom_style': 'full'}
    fname = eex_find_files.get_example_filename("lammps", "alkanes", molecule)

    # Read in the data
    dl = eex.datalayer.DataLayer(molecule)
    eex.translators.lammps.read_lammps_data_file(dl, fname, sim_data)

    # Write out the data
    oname = eex_find_files.get_scratch_directory(molecule)
    input_filename = oname + '.in'
    eex.translators.lammps.write_lammps_file(dl, oname, input_filename, unit_style="real")

    # Read in output data
    dl_new = eex.datalayer.DataLayer(molecule)
    eex.translators.lammps.read_lammps_data_file(dl_new, oname, sim_data)
    #assert eex.testing.dl_compare(dl, dl_new)

    ## write tests for unit conversions - requires lammps writer to write lammps in file
