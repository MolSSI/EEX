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
    assert dl.get_atom_count() == 648
    assert dl.get_atom_count("charge") == 648
    assert dl.get_atom_count("mass") == 648
    assert set(np.unique(atoms["atom_name"])) == {"H1", "H2", "O"}
    assert np.allclose(np.unique(atoms["charge"]), [-0.8476, 0.4238])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [1.008, 16.0])
    assert set(np.unique(atoms["residue_name"])) == {"WAT"}
    assert np.allclose(np.min(atoms["residue_index"]), 1)
    assert np.allclose(np.max(atoms["residue_index"]), 216)

    assert np.allclose(atoms[["X", "Y", "Z"]].min(), [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(), [19.388333, 19.213448, 19.423807])


# Test AMBER read by_index
def test_amber_spce_read_atoms_index(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms(
        ["atom_name", "charge", "atomic_number", "mass", "residue_name", "residue_index", "xyz"], by_value=False)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == {"H1", "H2", "O"}
    assert np.allclose(np.unique(atoms["charge"]), [0, 1])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [0, 1])
    assert set(np.unique(atoms["residue_name"])) == {0}
    assert np.allclose(np.min(atoms["residue_index"]), 1)
    assert np.allclose(np.max(atoms["residue_index"]), 216)
    assert np.allclose(atoms[["X", "Y", "Z"]].min(), [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(), [19.388333, 19.213448, 19.423807])


def test_amber_spce_read_bonds(spce_dl):
    data, dl = spce_dl
    # Test bond df
    bonds = dl.get_bonds()
    assert bonds.shape[0] == 648
    assert dl.get_term_count(2, "total") == 648

    assert set(np.unique(bonds["term_index"])) == {1, 2}
    assert set(dl.get_term_count(2)) == {1, 2, "total"}


# def test_amber_spce_parameters(spce_dl):
#     data, dl = spce_dl

#     # Test bond df
#     bonds = dl.get_bonds()
#     assert bonds.shape[0] == 648
#     assert set(np.unique(bonds["term_index"])) == set([1, 2])


@pytest.mark.parametrize("molecule", [
    "trappe_butane_single_molecule.prmtop",
    "trappe_propane_single_molecule.prmtop",
    "trappe_ethane_single_molecule.prmtop",
])
def test_amber_writer(molecule):
    fname = eex_find_files.get_example_filename("amber", "alkanes", molecule)

    # Read in the data
    dl = eex.datalayer.DataLayer(molecule)
    data = eex.translators.amber.read_amber_file(dl, fname)

    # Write out the data
    oname = eex_find_files.get_scratch_directory(molecule)
    eex.translators.amber.write_amber_file(dl, oname)

    # Read in output data
    dl_new = eex.datalayer.DataLayer(molecule)
    eex.translators.amber.read_amber_file(dl_new, oname)
    assert eex.testing.dl_compare(dl, dl_new)