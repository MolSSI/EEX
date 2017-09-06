"""
Tests the datalayer object for EEX
"""

import eex
import pytest
import pandas as pd
import numpy as np

_backend_list = ["HDF5", "Memory"]


def _build_atom_df(nmols):

    ncols = nmols * 3
    bond_df = pd.DataFrame(columns=["atom_index", "molecule_index", "atom_type", "charge", "X", "Y", "Z"])
    bond_df["atom_index"] = np.arange(ncols)
    bond_df["molecule_index"] = np.repeat(np.arange(nmols), 3)
    bond_df["atom_type"] = np.tile([1, 2, 2], nmols)
    bond_df["charge"] = np.tile([-0.8, 0.4, 0.4], nmols)
    bond_df["X"] = np.random.rand(ncols)
    bond_df["Y"] = np.random.rand(ncols)
    bond_df["Z"] = np.random.rand(ncols)
    return bond_df


@pytest.mark.parametrize("backend", _backend_list)
def test_df_bonds(backend):
    """
    Tests adding bonds as a DataFrame
    """

    dl = eex.datalayer.DataLayer("test_df_bonds", backend=backend)

    tmp_df = _build_atom_df(10)

    dl.add_atoms(tmp_df.loc[:5])
    dl.add_atoms(tmp_df.loc[5:])

    dl_df = dl.get_atoms(["molecule_index", "atom_type", "charge", "XYZ"])
    dl_df = dl_df.reset_index()

    tmp_df.equals(dl_df)


@pytest.mark.parametrize("backend", _backend_list)
def test_list_bonds(backend):
    """
    Tests adding bonds as a list
    """

    dl = eex.datalayer.DataLayer("test_list_bonds", backend=backend)

    tmp_df = _build_atom_df(10)
    for idx, row in tmp_df.iterrows():
        data = list(row)
        dl.add_atoms([data[0], data[1]], property_name="molecule_index")
        dl.add_atoms([data[0], data[2]], property_name="atom_type")
        dl.add_atoms([data[0], data[3]], property_name="charge")
        dl.add_atoms([data[0], data[4], data[5], data[6]], property_name="XYZ")

    dl_df = dl.get_atoms(["molecule_index", "atom_type", "charge", "XYZ"])
    dl_df = dl_df.reset_index()
    tmp_df.equals(dl_df)
