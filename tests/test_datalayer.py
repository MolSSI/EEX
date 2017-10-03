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
    bond_df.set_index("atom_index", inplace=True)
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

    dl.add_atoms(tmp_df.iloc[:5], by_value=True)
    dl.add_atoms(tmp_df.iloc[5:], by_value=True)

    dl_df = dl.get_atoms(["molecule_index", "atom_type", "charge", "XYZ"], by_value=True)

    # Compare DL df
    assert eex.testing.df_compare(tmp_df, dl_df)

    # Compare rand DF
    dl_rand_df = dl.get_other("atoms")
    assert eex.testing.df_compare(rand_df, dl_rand_df)


def test_add_parameters():

    dl = eex.datalayer.DataLayer("test_add_parameters")
    two_body_md = eex.metadata.get_term_metadata(2, "forms", "harmonic")
    three_body_md = eex.metadata.get_term_metadata(3, "forms", "harmonic")

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


def test_add_parameters_units():

    dl = eex.datalayer.DataLayer("test_add_parameters_units")
    two_body_md = eex.metadata.get_term_metadata(2, "forms", "harmonic")
    three_body_md = eex.metadata.get_term_metadata(3, "forms", "harmonic")

    # Test unit types

    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])

    # Same units
    utype_2b = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    utype_2bl = [utype_2b["K"], utype_2b["R0"]]
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2b)
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2bl)
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0], utype=utype_2b)
    assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0], utype=utype_2bl)

    # Scale by 2
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "0.5 * angstrom"}
    utype_2bl = [utype_2b["K"], utype_2b["R0"]]
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 8.0, "R0": 10.0}, utype=utype_2b)
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 8.0, "R0": 10.0}, utype=utype_2bl)
    assert 0 == dl.add_parameters(2, "harmonic", [8.0, 10.0], utype=utype_2b)
    assert 0 == dl.add_parameters(2, "harmonic", [8.0, 10.0], utype=utype_2bl)

    # Different unit type
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "picometers"}
    assert 0 == dl.add_parameters(2, "harmonic", [8.0, 500.0], utype=utype_2b)

    # Finally a different parameter type all together
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "picometers"}
    assert 1 == dl.add_parameters(2, "harmonic", [8.0, 5.0], utype=utype_2b)

def test_get_parameters():

    dl = eex.datalayer.DataLayer("test_get_parameters")

    # Add a few parameters
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 4.0, "R0": 5.0})
    assert 1 == dl.add_parameters(2, "harmonic", {"K": 6.0, "R0": 7.0})

    parm1 = dl.get_parameters(2, 0)
    assert parm1[0] == "harmonic"
    assert eex.testing.dict_compare(parm1[1], {"K": 4.0, "R0": 5.0})
    assert eex.testing.dict_compare(parm1[1], {"K": 4.0, "R0": 5.0 + 1.e-12})

    parm2 = dl.get_parameters(2, 1)
    assert parm2[0] == "harmonic"
    assert eex.testing.dict_compare(parm2[1], {"K": 6.0, "R0": 7.0})

    with pytest.raises(KeyError):
        dl.get_parameters(2, 1231234123)


def test_get_parameters_units():

    dl = eex.datalayer.DataLayer("test_get_parameters_units")

    utype_2b = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    assert 0 == dl.add_parameters(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2b)
    assert 1 == dl.add_parameters(2, "harmonic", {"K": 6.0, "R0": 7.0}, utype=utype_2b)

    utype_2scale = {"K": "2.0 * (kJ / mol) * angstrom ** -2", "R0": "2.0 * angstrom"}
    parm1 = dl.get_parameters(2, 0, utype=utype_2scale)
    assert parm1[0] == "harmonic"
    assert eex.testing.dict_compare(parm1[1], {"K": 2.0, "R0": 2.5})

    parm2 = dl.get_parameters(2, 0, utype=[utype_2scale["K"], utype_2scale["R0"]])
    assert parm2[0] == "harmonic"
    assert eex.testing.dict_compare(parm2[1], {"K": 2.0, "R0": 2.5})

    with pytest.raises(TypeError):
        dl.get_parameters(2, 0, utype=set([5, 6]))

def test_atom_units():
    """
    Tests adding bonds as a DataFrame
    """

    dl = eex.datalayer.DataLayer("test_df_bonds", backend="memory")
    utypes = {"charge": "elementary_charge", "XYZ": "picometers"}

    tmp_df = _build_atom_df(2)
    tmp_df.index.name = "atom_index"
    dl.add_atoms(tmp_df, by_value=True, utype=utypes)

    assert eex.testing.df_compare(tmp_df, dl.get_atoms("charge", by_value=True), columns="charge")

    # Check multiple output types
    xyz_pm = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "picometers"})
    assert eex.testing.df_compare(xyz_pm, tmp_df[["X", "Y", "Z"]])

    # Check twice incase we accidently scale internal data
    xyz_pm = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "picometers"})
    assert eex.testing.df_compare(xyz_pm, tmp_df[["X", "Y", "Z"]])

    # Check default and specified
    xyz_ang1 = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "angstrom"})
    xyz_ang2 = dl.get_atoms("XYZ", by_value=True)

    assert eex.testing.df_compare(xyz_ang1, xyz_ang2)
    assert eex.testing.df_compare(xyz_ang1, xyz_pm * 0.01)
