"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import eex_find_files


@pytest.fixture(scope="module", params=["HDF5", "Memory"])
def spce_dl(request):
    fname = eex_find_files.get_example_filename("amber", "water", "spce.prmtop")
    dl = eex.datalayer.DataLayer("test_amber_read", backend=request.param)
    data = eex.translators.amber.read_amber_file(dl, fname, blocksize=55)
    # print("-- Built AMBER SPCE DL --")
    yield (data, dl)
    dl.close()


def test_amber_spce_read_data(spce_dl):
    data, dl = spce_dl
    assert data["VERSION"] == "V0001.000"


# Test AMBER read by_value
def test_amber_spce_read_atoms_value(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms(
        ["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index", "xyz"], by_value=True)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == set(["H1", "H2", "O"])
    assert np.allclose(np.unique(atoms["charge"]), [-0.8476, 0.4238])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [1.008, 16.0])
    assert set(np.unique(atoms["residue_name"])) == set(["WAT"])
    assert np.allclose(np.min(atoms["residue_index"]), 0)
    assert np.allclose(np.max(atoms["residue_index"]), 215)

    assert np.allclose(atoms[["X", "Y", "Z"]].min(), [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(), [19.388333, 19.213448, 19.423807])


# Test AMBER read by_index
def test_amber_spce_read_atoms_index(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms(
        ["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index", "xyz"], by_value=False)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == set(["H1", "H2", "O"])
    assert np.allclose(np.unique(atoms["charge"]), [0, 1])
    assert np.allclose(np.unique(atoms["atomic_number"]), [0, 1])
    assert np.allclose(np.unique(atoms["mass"]), [0, 1])
    assert set(np.unique(atoms["residue_name"])) == set([0])
    assert np.allclose(np.min(atoms["residue_index"]), 0)
    assert np.allclose(np.max(atoms["residue_index"]), 215)
    assert np.allclose(atoms[["X", "Y", "Z"]].min(), [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(), [19.388333, 19.213448, 19.423807])


def test_amber_spce_read_bonds(spce_dl):
    data, dl = spce_dl
    # Test bond df
    bonds = dl.get_bonds()
    assert bonds.shape[0] == 648
    assert set(np.unique(bonds["term_index"])) == set([1, 2])


# def test_amber_spce_parameters(spce_dl):
#     data, dl = spce_dl

#     # Test bond df
#     bonds = dl.get_bonds()
#     assert bonds.shape[0] == 648
#     assert set(np.unique(bonds["term_index"])) == set([1, 2])
