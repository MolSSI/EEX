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
    fname = eex_find_files.get_example_filename("amber", "water/spce.prmtop")
    dl = eex.datalayer.DataLayer("test_amber_read", backend=request.param)
    data = eex.translators.amber.read_amber_file(dl, fname)
    yield (data, dl)
    del dl


def test_amber_spce_read_data(spce_dl):
    data, dl = spce_dl
    assert data["VERSION"] == "V0001.000"


# Test AMBER read by_value
def test_amber_spce_read_atoms_value(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms(
        ["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index"], by_value=True)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == set(["H1", "H2", "O"])
    assert np.allclose(np.unique(atoms["charge"]), [-1.54452215E+01, 7.72261074E+00])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [1.008, 16.0])
    assert set(np.unique(atoms["residue_name"])) == set(["WAT"])
    assert np.allclose(np.min(atoms["residue_index"]), 0)
    assert np.allclose(np.max(atoms["residue_index"]), 215)


# Test AMBER read by_index
def test_amber_spce_read_atoms_value(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms(
        ["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index"], by_value=False)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == set(["H1", "H2", "O"])
    assert np.allclose(np.unique(atoms["charge"]), [0, 1])
    assert np.allclose(np.unique(atoms["atomic_number"]), [0, 1])
    assert np.allclose(np.unique(atoms["mass"]), [0, 1])
    assert set(np.unique(atoms["residue_name"])) == set([0])
    assert np.allclose(np.min(atoms["residue_index"]), 0)
    assert np.allclose(np.max(atoms["residue_index"]), 215)