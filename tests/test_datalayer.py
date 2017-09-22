"""
Tests the datalayer object for EEX
"""

import eex
import pytest
import pandas as pd
import numpy as np

# Set the Seed
np.random.seed(0)

# Any parameters to loop over
_backend_list = ["HDF5", "Memory"]


def _build_atom_df(nmols):

    ncols = nmols * 3
    bond_df = pd.DataFrame(columns=["atom_index", "molecule_index", "atom_type", "charge", "X", "Y", "Z"])
    bond_df["atom_index"] = np.arange(ncols)
    bond_df["molecule_index"] = np.repeat(np.arange(nmols), 3)
    bond_df["atom_type"] = np.tile(["O", "H", "H"], nmols)
    bond_df["charge"] = np.tile([-0.8, 0.4, 0.4], nmols)
    bond_df["X"] = np.random.rand(ncols)
    bond_df["Y"] = np.random.rand(ncols)
    bond_df["Z"] = np.random.rand(ncols)
    return bond_df


@pytest.mark.parametrize("backend", _backend_list)
def test_df_atoms(backend):
    """
    Tests adding bonds as a DataFrame
    """

    dl = eex.datalayer.DataLayer("test_df_bonds", backend=backend)

    tmp_df = _build_atom_df(10)

    # Add a random DataFrame with name "atoms" to check interference with "add_atoms"
    rand_df = pd.DataFrame(np.random.rand(4, 4))
    dl.add_other("atoms", rand_df)

    dl.add_atoms(tmp_df.loc[:5], by_value=True)
    dl.add_atoms(tmp_df.loc[5:], by_value=True)

    dl_df = dl.get_atoms(["molecule_index", "atom_type", "charge", "XYZ"], by_value=True)
    dl_df = dl_df.reset_index()

    # Compare DL df
    tmp_df.equals(dl_df)

    # Compare rand DF
    dl_rand_df = dl.get_other("atoms")
    rand_df.equals(dl_rand_df)


@pytest.mark.parametrize("backend", _backend_list)
def test_list_atoms(backend):
    """
    Tests adding bonds as a list
    """

    dl = eex.datalayer.DataLayer("test_list_bonds", backend=backend)

    tmp_df = _build_atom_df(10)
    for idx, row in tmp_df.iterrows():
        data = list(row)
        dl.add_atoms([data[0], data[1]], property_name="molecule_index", by_value=True)
        dl.add_atoms([data[0], data[2]], property_name="atom_type", by_value=True)
        dl.add_atoms([data[0], data[3]], property_name="charge", by_value=True)
        dl.add_atoms([data[0], data[4], data[5], data[6]], property_name="XYZ", by_value=True)

    dl_df = dl.get_atoms(["molecule_index", "atom_type", "charge", "XYZ"], by_value=True)
    dl_df = dl_df.reset_index()
    tmp_df.equals(dl_df)


def test_register_functional_forms():

    dl = eex.datalayer.DataLayer("test_functional_forms")

    # Add a few functional forms
    dl.register_functional_forms(2, "harmonic", eex.metadata.get_term_metadata(2, "forms", "harmonic"))
    dl.register_functional_forms(2, "fene", eex.metadata.get_term_metadata(2, "forms", "fene"))
    dl.register_functional_forms(3, "harmonic", eex.metadata.get_term_metadata(3, "forms", "harmonic"))

    # Bound ineligable
    with pytest.raises(KeyError):
        dl.register_functional_forms(2, "harmonic", eex.metadata.get_term_metadata(2, "forms", "harmonic"))

    with pytest.raises(KeyError):
        dl.register_functional_forms("turle", "turle", {})

    with pytest.raises(KeyError):
        dl.register_functional_forms(2, "harmonic2", {"form": "(x-x0)**2"})


def test_add_parameters():

    dl = eex.datalayer.DataLayer("test_functional_forms")

    # Add a few functional forms
    two_body_md = eex.metadata.get_term_metadata(2, "forms", "harmonic")
    three_body_md = eex.metadata.get_term_metadata(3, "forms", "harmonic")
    dl.register_functional_forms(2, "harmonic", two_body_md)
    dl.register_functional_forms(3, "harmonic", three_body_md)

    # Check duplicates
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0 + 1.e-10])

    # New parameters
    assert 1 == dl.add_parameters(2, "harmonic", [4.0, 6.0])
    assert 0 == dl.add_parameters(3, "harmonic", [4.0, 6.0])

    # Check uid types
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0], uid=0)
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0 + 1.e-10], uid=0)

    # Add out of order uid
    assert 10 == dl.add_parameters(2, "harmonic", [9.0, 9], uid=10)
    assert 2 == dl.add_parameters(2, "harmonic", [9.0, 10.0])

    # Do we want to allow forced dups?
    assert 15 == dl.add_parameters(2, "harmonic", [22.0, 22.0], uid=15)
    assert 11 == dl.add_parameters(2, "harmonic", [22.0, 22.0], uid=11)

    # Check add by dict
    mdp = two_body_md["parameters"]
    assert 0 == dl.add_parameters(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 5.0})
    assert 1 == dl.add_parameters(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 6.0})
    assert 3 == dl.add_parameters(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 7.0})

    with pytest.raises(KeyError):
        dl.add_parameters(2, "harmonic", {mdp[0]: 4.0, "turtle": 5.0})

    mdp = three_body_md["parameters"]
    assert 0 == dl.add_parameters(3, "harmonic", {mdp[0]: 4.0, mdp[1]: 6.0})

    # Check uid type
    with pytest.raises(TypeError):
        dl.add_parameters(2, "harmonic", [11.0, 9.0], uid="turtle")

    # Validate parameter data
    with pytest.raises(ValueError):
        dl.add_parameters(2, "harmonic", [4.0, 5.0, 3.0])

    with pytest.raises(TypeError):
        dl.add_parameters(2, "harmonic", [11.0, "duck"])

    # Check name errors
    with pytest.raises(KeyError):
        dl.add_parameters("turtle", "harmonic", [4.0, 5.0, 3.0])

    with pytest.raises(KeyError):
        dl.add_parameters(2, "harmonic_abc", [4.0, 5.0])
