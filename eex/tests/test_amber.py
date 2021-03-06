"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import copy
from . import eex_find_files


@pytest.fixture(scope="module", params=["HDF5", "Memory"])
def spce_dl(request):
    fname = eex_find_files.get_example_filename("amber", "water",
                                                "spce.prmtop")
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

    atoms = dl.get_atoms([
        "atom_name", "charge", "atomic_number", "mass", "residue_name",
        "residue_index", "xyz"
    ],
                         by_value=True)
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

    assert np.allclose(atoms[["X", "Y", "Z"]].min(),
                       [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(),
                       [19.388333, 19.213448, 19.423807])


# Test AMBER read by_index
def test_amber_spce_read_atoms_index(spce_dl):
    data, dl = spce_dl

    atoms = dl.get_atoms([
        "atom_name", "charge", "atomic_number", "mass", "residue_name",
        "residue_index", "xyz"
    ],
                         by_value=False)
    assert atoms.shape[0] == 648
    assert set(np.unique(atoms["atom_name"])) == {"H1", "H2", "O"}
    assert np.allclose(np.unique(atoms["charge"]), [0, 1])
    assert np.allclose(np.unique(atoms["atomic_number"]), [1.0, 8.0])
    assert np.allclose(np.unique(atoms["mass"]), [0, 1])
    assert set(np.unique(atoms["residue_name"])) == {0}
    assert np.allclose(np.min(atoms["residue_index"]), 1)
    assert np.allclose(np.max(atoms["residue_index"]), 216)
    assert np.allclose(atoms[["X", "Y", "Z"]].min(),
                       [-0.757235, -0.519927, -0.872856])
    assert np.allclose(atoms[["X", "Y", "Z"]].max(),
                       [19.388333, 19.213448, 19.423807])


def test_amber_spce_read_bonds(spce_dl):
    data, dl = spce_dl
    # Test bond df
    bonds = dl.get_bonds()
    assert bonds.shape[0] == 648
    assert dl.get_term_count(2, "total") == 648

    assert set(np.unique(bonds["term_index"])) == {1, 2}
    assert set(dl.get_term_count(2)) == {1, 2, "total"}


def test_amber_read_butane():
    fname = eex_find_files.get_example_filename(
        "amber", "alkanes", "trappe_butane_single_molecule.prmtop")
    dl = eex.datalayer.DataLayer("test_amber_read", backend="memory")
    data = eex.translators.amber.read_amber_file(dl, fname)

    assert dl.get_atom_count() == 4
    assert dl.get_bond_count() == 3

    box_info = dl.get_box_size(
        utype={
            'a': 'angstrom',
            'b': 'angstrom',
            'c': 'angstrom',
            'alpha': 'degree',
            'beta': 'degree',
            'gamma': 'degree'
        })

    ref_box = {
        "a": 100,
        "b": 100,
        "c": 100,
        "alpha": 90,
        "beta": 90,
        "gamma": 90
    }

    assert eex.testing.dict_compare(box_info, ref_box)


@pytest.mark.parametrize("backend", ["HDF5", "Memory"])
@pytest.mark.parametrize("molecule", [
    "trappe_butane_single_molecule.prmtop",
    "trappe_propane_single_molecule.prmtop",
    "trappe_ethane_single_molecule.prmtop",
])
def test_amber_writer(molecule, backend):
    fname = eex_find_files.get_example_filename("amber", "alkanes", molecule)

    # Read in the data
    dl = eex.datalayer.DataLayer(molecule, backend=backend)
    data = eex.translators.amber.read_amber_file(dl, fname)

    # Write out the data
    oname = eex_find_files.get_scratch_directory(molecule)
    eex.translators.amber.write_amber_file(dl, oname)

    # Read in output data
    dl_new = eex.datalayer.DataLayer(molecule)
    eex.translators.amber.read_amber_file(dl_new, oname)

    # Compare the two datalayers
    assert eex.testing.dl_compare(dl, dl_new)

    # Since dl_compare does not yet do nonbonds, get NB from datalayer and compare
    assert (dl.list_stored_nb_types() == ["LJ"])
    assert (dl.list_nb_parameters(nb_name="LJ") == dl_new.list_nb_parameters(
        nb_name="LJ"))


def test_amber_compatibility_NB_number(butane_dl):

    dl = butane_dl(nb=False)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    dl.add_nb_parameter(
        atom_type=1,
        atom_type2=2,
        nb_name="LJ",
        nb_model="epsilon/sigma",
        nb_parameters={
            'sigma': 3.75,
            'epsilon': 0.1947460018
        },
        utype={
            'sigma': 'angstrom',
            'epsilon': 'kcal * mol ** -1'
        })

    with pytest.raises(ValueError):
        eex.translators.amber.write_amber_file(dl, oname)


def test_amber_compatibility_NB_type(butane_dl):

    dl = butane_dl(nb=False)

    dl.add_nb_parameter(
        atom_type=1,
        nb_name="Buckingham",
        nb_model=None,
        nb_parameters=[1.0, 1.0, 1.0])

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    with pytest.raises(KeyError):
        eex.translators.amber.write_amber_file(dl, oname)


def test_amber_compatibility_check_mixing_rule(butane_dl):

    # Get butane topology
    dl = butane_dl(nb=True)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    # Check to make sure that mixing rule will be applied by compatibility check
    eex.translators.amber.write_amber_file(dl, oname)

    # Make sure we're getting pair parameters from DL
    nb_pairs = dl.list_nb_parameters(nb_name="LJ", nb_model="AB", itype="pair")

    assert (len(nb_pairs.keys()) == 3)


def test_amber_compatibility_functional_form(butane_dl):

    dl = butane_dl(ff=False)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    # Add incompatible functional form
    with pytest.raises(TypeError):
        dl.add_term_parameter(2, "fene", {
            "K": 1,
            "R0": 1,
            "epsilon": 1,
            "sigma": 1,
        })
        eex.translators.amber.write_amber_file(dl, oname)


def test_amber_compatibility_scaling(butane_dl):

    dl = butane_dl(scale=False)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    # Set with noncompatible scale13
    scaling_factors = {
        "coul": {
            "scale12": 0.0,
            "scale13": 0.50,
            "scale14": 0.75,
        },
        "vdw": {
            "scale12": 0.0,
            "scale13": 0.0,
            "scale14": 0.75,
        }
    }

    dl.set_nb_scaling_factors(scaling_factors)

    with pytest.raises(ValueError):
        eex.translators.amber.write_amber_file(dl, oname)


def test_amber_compatibility_no_scaling(butane_dl):

    dl = butane_dl(scale=False)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    with pytest.raises(ValueError):
        eex.translators.amber.write_amber_file(dl, oname)
