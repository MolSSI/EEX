"""
Tests the datalayer object for EEX
"""

import eex
import pytest
import pandas as pd
import numpy as np

def build_atom_df(nmols):

    ncols = nmols * 3
    bond_df = pd.DataFrame(columns=["atom_index", "molecule_index", "atom_type", "X", "Y", "Z"])
    bond_df["atom_index"] = np.arange(ncols)
    bond_df["molecule_index"] = np.repeat(np.arange(nmols), 3)
    bond_df["atom_type"] = np.tile([1, 2, 2], nmols)
    bond_df["X"] = np.random.rand(ncols)
    bond_df["Y"] = np.random.rand(ncols)
    bond_df["Z"] = np.random.rand(ncols)
    return bond_df

def test_df_bonds():
    dl = eex.datalayer.DataLayer("test_df_bonds")

    tmp_df = build_atom_df(10)

    dl.add_atoms(tmp_df.loc[:5])
    dl.add_atoms(tmp_df.loc[5:])

    dl_df = dl.read_atoms()

    tmp_df.equals(dl_df)

def test_list_bonds():
    dl = eex.datalayer.DataLayer("test_list_bonds")

    tmp_df = build_atom_df(10)
    for idx, row in tmp_df.iterrows():
        dl.add_atoms(list(row))

    dl_df = dl.read_atoms()
    tmp_df.equals(dl_df)
    
