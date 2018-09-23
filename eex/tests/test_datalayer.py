"""
Tests the datalayer object for EEX
"""

import eex
import pytest
import pandas as pd
import numpy as np

from eex.testing import df_compare, dict_compare

from eex import nb_converter

# Set the Seed
np.random.seed(0)

# Any parameters to loop over
_backend_list = ["HDF5", "Memory"]


def _build_atom_df(nmols):
    """
    Builds a random testing DataFrame
    """

    ncols = nmols * 3
    bond_df = pd.DataFrame()
    bond_df["atom_index"] = np.arange(ncols)
    bond_df["molecule_index"] = np.repeat(np.arange(nmols), 3)
    bond_df["atom_name"] = np.tile(["O", "H", "H"], nmols)
    bond_df["atom_type"] = np.tile([1, 2, 2], nmols)
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

    # Test count functionality works
    assert dl.get_atom_count() == 30
    assert dl.get_atom_count("charge") == 30
    assert dl.get_atom_count("xyz") == 30

    dl_df = dl.get_atoms(None, by_value=True)

    # Compare DL df
    assert df_compare(tmp_df, dl_df)

    # Compare rand DF
    dl_rand_df = dl.get_other("atoms")
    assert df_compare(rand_df, dl_rand_df)


@pytest.mark.parametrize("backend", _backend_list)
def test_misc(backend):
    dl = eex.datalayer.DataLayer("test_misc_functions", backend=backend)

    atoms = _build_atom_df(2)

    dl.add_atoms(atoms)

    with pytest.raises(AttributeError):
        dl.call_by_string("not_a_function")

    assert(set(dl.list_tables()) == {'atom_name', 'xyz', 'molecule_index', 'charge', 'atom_type'})
    assert(dl.list_other_tables() == [])


def test_add_atom_parameter():
    dl = eex.datalayer.DataLayer("test_add_atom_parameters")

    # Test duplicate and uid add of same
    assert 0 == dl.add_atom_parameter("mass", 5.0)
    assert 0 == dl.add_atom_parameter("mass", 5.0)
    assert 0 == dl.add_atom_parameter("mass", 5.0, uid=0)

    # Once more
    assert 1 == dl.add_atom_parameter("mass", 6.0)
    assert 1 == dl.add_atom_parameter("mass", 6.0, uid=1)

    # Test out of order adds
    assert 2 == dl.add_atom_parameter("mass", 7.0, uid=2)
    assert 5 == dl.add_atom_parameter("mass", 8.0, uid=5)
    assert 3 == dl.add_atom_parameter("mass", 9.0)
    assert 3 == dl.add_atom_parameter("mass", {"mass": 9.0})

    # Test possible raises
    with pytest.raises(KeyError):
        dl.add_atom_parameter("mass", 6.0, uid=0)


def test_add_atom_parameter_units():
    dl = eex.datalayer.DataLayer("test_add_atom_parameters")

    # Test duplicate and uid add of same
    assert 0 == dl.add_atom_parameter("mass", 5.0, utype="kilogram / mol")
    assert pytest.approx(5.0) == dl.get_atom_parameter("mass", 0, utype="kilogram / mol")
    assert pytest.approx(5000.0) == dl.get_atom_parameter("mass", 0, utype="gram / mol")
    assert pytest.approx(5000.0) == dl.get_atom_parameter("mass", 0, utype={"mass": "gram / mol"})


def test_add_term_parameter():
    """
    Test adding parameters to the DL object
    """

    dl = eex.datalayer.DataLayer("test_add_term_parameters")
    two_body_md = eex.metadata.get_term_metadata(2, "forms", "harmonic")
    three_body_md = eex.metadata.get_term_metadata(3, "forms", "harmonic")

    # Check duplicates
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0 + 1.e-10])

    # New parameters
    assert 1 == dl.add_term_parameter(2, "harmonic", [4.0, 6.0])
    assert 0 == dl.add_term_parameter(3, "harmonic", [4.0, 6.0])

    # Check uid types
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0], uid=0)
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0 + 1.e-10], uid=0)

    # Add out of order uid
    assert 10 == dl.add_term_parameter(2, "harmonic", [9.0, 9], uid=10)
    assert 2 == dl.add_term_parameter(2, "harmonic", [9.0, 10.0])

    # Do we want to allow forced dups?
    assert 15 == dl.add_term_parameter(2, "harmonic", [22.0, 22.0], uid=15)
    assert 11 == dl.add_term_parameter(2, "harmonic", [22.0, 22.0], uid=11)

    # Check add by dict
    mdp = two_body_md["parameters"]
    assert 0 == dl.add_term_parameter(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 5.0})
    assert 1 == dl.add_term_parameter(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 6.0})
    assert 3 == dl.add_term_parameter(2, "harmonic", {mdp[0]: 4.0, mdp[1]: 7.0})

    with pytest.raises(KeyError):
        dl.add_term_parameter(2, "harmonic", {mdp[0]: 4.0, "turtle": 5.0})

    mdp = three_body_md["parameters"]
    assert 0 == dl.add_term_parameter(3, "harmonic", {mdp[0]: 4.0, mdp[1]: 6.0})

    # Check uid type
    with pytest.raises(TypeError):
        dl.add_term_parameter(2, "harmonic", [11.0, 9.0], uid="turtle")

    # Validate parameter data
    with pytest.raises(ValueError):
        dl.add_term_parameter(2, "harmonic", [4.0, 5.0, 3.0])

    with pytest.raises(TypeError):
        dl.add_term_parameter(2, "harmonic", [11.0, "duck"])

    # Check name errors
    with pytest.raises(KeyError):
        dl.add_term_parameter("turtle", "harmonic", [4.0, 5.0, 3.0])

    with pytest.raises(KeyError):
        dl.add_term_parameter(2, "harmonic_abc", [4.0, 5.0])


def test_add_term_parameters_units():
    """
    Test adding parameters to the DL object with units
    """

    dl = eex.datalayer.DataLayer("test_add_parameters_units")
    two_body_md = eex.metadata.get_term_metadata(2, "forms", "harmonic")
    three_body_md = eex.metadata.get_term_metadata(3, "forms", "harmonic")

    # Test unit types

    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0])
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0])

    # Same units
    utype_2b = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    utype_2bl = [utype_2b["K"], utype_2b["R0"]]
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2b)
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2bl)
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0], utype=utype_2b)
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0], utype=utype_2bl)

    # Scale by 2
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "0.5 * angstrom"}
    utype_2bl = [utype_2b["K"], utype_2b["R0"]]
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 8.0, "R0": 10.0}, utype=utype_2b)
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 8.0, "R0": 10.0}, utype=utype_2bl)
    assert 0 == dl.add_term_parameter(2, "harmonic", [8.0, 10.0], utype=utype_2b)
    assert 0 == dl.add_term_parameter(2, "harmonic", [8.0, 10.0], utype=utype_2bl)

    # Different unit type
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "picometers"}
    assert 0 == dl.add_term_parameter(2, "harmonic", [8.0, 500.0], utype=utype_2b)

    # Finally a different parameter type all together
    utype_2b = {"K": "0.5 * (kJ / mol) * angstrom ** -2", "R0": "picometers"}
    assert 1 == dl.add_term_parameter(2, "harmonic", [8.0, 5.0], utype=utype_2b)


def test_get_term_parameter():
    """
    Test obtaining parameters from the DL
    """

    dl = eex.datalayer.DataLayer("test_get_parameters")

    # Add a few parameters
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 4.0, "R0": 5.0})
    assert 1 == dl.add_term_parameter(2, "harmonic", {"K": 6.0, "R0": 7.0})

    parm1 = dl.get_term_parameter(2, 0)
    assert parm1[0] == "harmonic"
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0})
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0 + 1.e-12})

    parm2 = dl.get_term_parameter(2, 1)
    assert parm2[0] == "harmonic"
    assert dict_compare(parm2[1], {"K": 6.0, "R0": 7.0})

    parm1 = dl.get_term_parameter(2, 0, ftype='harmonic')
    assert parm1[0] == "harmonic"
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0})
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0 + 1.e-12})

    utype = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    parm1 = dl.get_term_parameter(2, 0, utype=utype)
    assert parm1[0] == "harmonic"
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0})
    assert dict_compare(parm1[1], {"K": 4.0, "R0": 5.0 + 1.e-12})

    with pytest.raises(KeyError):
        dl.get_term_parameter(2, 1231234123)

    with pytest.raises(KeyError):
        dl.get_term_parameter(2, 0, utype=utype, ftype="harmic")


def test_get_term_parameters_units():
    """
    Test obtaining parameters from the DL with units
    """

    dl = eex.datalayer.DataLayer("test_get_parameters_units")

    utype_2b = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2b)
    assert 1 == dl.add_term_parameter(2, "harmonic", {"K": 6.0, "R0": 7.0}, utype=utype_2b)

    utype_2scale = {"K": "2.0 * (kJ / mol) * angstrom ** -2", "R0": "2.0 * angstrom"}
    parm1 = dl.get_term_parameter(2, 0, utype=utype_2scale, ftype="harmonic")
    assert parm1[0] == "harmonic"
    assert dict_compare(parm1[1], {"K": 2.0, "R0": 2.5})

    parm2 = dl.get_term_parameter(2, 0, utype=[utype_2scale["K"], utype_2scale["R0"]], ftype="harmonic")
    assert parm2[0] == "harmonic"
    assert dict_compare(parm2[1], {"K": 2.0, "R0": 2.5})

    with pytest.raises(TypeError):
        dl.get_term_parameter(2, 0, utype={5, 6}, ftype="harmonic")

    with pytest.raises(KeyError):
        utype = {"K": "(kJ / mol) * angstrom ** -2"}
        dl.get_term_parameter(2, 0, utype=utype, ftype="harmonic")


def test_list_parameters():
    """
    Test listing parameters uids
    """

    dl = eex.datalayer.DataLayer("test_list_parameters")

    # Add two-body
    assert 0 == dl.add_term_parameter(2, "harmonic", [4.0, 5.0])
    assert 1 == dl.add_term_parameter(2, "harmonic", [4.0, 6.0])
    assert 9 == dl.add_term_parameter(2, "harmonic", [6.0, 7.0], uid=9)

    # Add three-body
    assert 0 == dl.add_term_parameter(3, "harmonic", [4.0, 5.0])
    assert 1 == dl.add_term_parameter(3, "harmonic", [4.0, 6.0])
    assert 9 == dl.add_term_parameter(3, "harmonic", [6.0, 7.0], uid=9)

    full_list = dl.list_term_uids()
    assert {2, 3, 4} == set(full_list)
    assert {0, 1, 9} == set(full_list[2])
    assert {0, 1, 9} == set(full_list[3])
    assert set() == set(full_list[4])

    bond_list = dl.list_term_uids(2)
    assert {0, 1, 9} == set(bond_list)


def test_atom_units():
    """
    Tests adding bonds as a DataFrame
    """

    dl = eex.datalayer.DataLayer("test_df_bonds", backend="memory")
    utypes = {"charge": "elementary_charge", "XYZ": "picometers"}

    tmp_df = _build_atom_df(2)
    tmp_df.index.name = "atom_index"
    dl.add_atoms(tmp_df, by_value=True, utype=utypes)

    assert df_compare(tmp_df, dl.get_atoms("charge", by_value=True), columns="charge")

    # Check multiple output types
    xyz_pm = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "picometers"})
    assert df_compare(xyz_pm, tmp_df[["X", "Y", "Z"]])

    # Check twice incase we accidently scale internal data
    xyz_pm = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "picometers"})
    assert df_compare(xyz_pm, tmp_df[["X", "Y", "Z"]])

    # Check default and specified
    xyz_ang1 = dl.get_atoms("XYZ", by_value=True, utype={"xyz": "angstrom"})
    xyz_ang2 = dl.get_atoms("XYZ", by_value=True)

    assert df_compare(xyz_ang1, xyz_ang2)
    assert df_compare(xyz_ang1, xyz_pm * 0.01)


def test_box_size():
    dl = eex.datalayer.DataLayer("test_box_size", backend="memory")

    tmp = {"a": 5, "b": 6, "c": 7, "alpha": np.pi / 3, "beta": 0.0, "gamma": 0.0}

    # Normal set/get
    dl.set_box_size(tmp)
    comp = dl.get_box_size()
    assert dict_compare(tmp, comp)

    # Set/get with units
    utype = {"a": "nanometers", "b": "nanometers", "c": "nanometers", "alpha": "radian", "beta": "radian", "gamma": "radian"}
    dl.set_box_size(tmp, utype=utype)
    comp = dl.get_box_size(utype=utype)
    assert np.isclose(comp["a"], tmp["a"])
    assert np.isclose(comp["b"], tmp["b"])
    assert np.isclose(comp["c"], tmp["c"])
    assert np.isclose(comp["alpha"], tmp["alpha"])
    assert np.isclose(comp["beta"], tmp["beta"])
    assert np.isclose(comp["gamma"], tmp["gamma"])

    # Set/get with units
    utype_1 = {"a": "angstroms", "b": "angstroms", "c": "angstroms", "alpha": "degree", "beta": "degree", "gamma": "degree"}
    dl.set_box_size(tmp, utype=utype)
    comp = dl.get_box_size(utype=utype_1)
    assert np.isclose(comp["a"], tmp["a"] * 10)
    assert np.isclose(comp["b"], tmp["b"] * 10)
    assert np.isclose(comp["c"], tmp["c"] * 10)
    assert np.isclose(comp["alpha"], tmp["alpha"] * 180.0 / np.pi)
    assert np.isclose(comp["beta"], tmp["beta"] * 180.0 / np.pi)
    assert np.isclose(comp["gamma"], tmp["gamma"] * 180.0 / np.pi)

    utype_2 = {"a": "miles", "b": "miles", "c": "miles", "alpha": "radians", "beta": "radians", "gamma": "radians"}
    with pytest.raises(AssertionError):
        dl.set_box_size(tmp, utype=utype_2)
        comp = dl.get_box_size()
        assert dict_compare(tmp, comp)


def test_box_center():
    dl = eex.datalayer.DataLayer("test_box_center", backend="memory")

    tmp = {"x": 5, "y": 6, "z": 7}

    # Normal set/get
    dl.set_box_center(tmp)
    comp = dl.get_box_center()
    assert dict_compare(tmp, comp)

    # Check failure
    with pytest.raises(KeyError):
        dl.set_box_center({'x': tmp['x']})

    # Set/get with units
    utype = {"x": "nanometers", "y": "nanometers", "z": "nanometers"}
    dl.set_box_center(tmp, utype=utype)
    comp = dl.get_box_center(utype=utype)
    assert np.isclose(comp["x"], tmp["x"])
    assert np.isclose(comp["y"], tmp["y"])
    assert np.isclose(comp["z"], tmp["z"])

    # Set/get with units
    utype_1 = {"x": "angstroms", "y": "angstroms", "z": "angstroms"}
    dl.set_box_center(tmp, utype=utype)
    comp = dl.get_box_center(utype=utype_1)
    assert np.isclose(comp["x"], tmp["x"] * 10)
    assert np.isclose(comp["y"], tmp["y"] * 10)
    assert np.isclose(comp["z"], tmp["z"] * 10)

    utype_2 = {"x": "miles", "y": "miles", "z": "miles"}
    with pytest.raises(AssertionError):
        dl.set_box_center(tmp, utype=utype_2)
        comp = dl.get_box_center()
        assert dict_compare(tmp, comp)


def test_add_nb_parameter():

    # Create empty data layer
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model="AB", nb_parameters=[1.0, 1.0])
    dl.add_nb_parameter(atom_type=2, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters=[1.0, 1.0])

    # Add AB LJ parameters to data layer - add to two atoms
    dl.add_nb_parameter(atom_type=1, atom_type2=2, nb_name="LJ", nb_model="AB", nb_parameters=[2.0, 2.0])

    # Grab stored test parameters - will need to replace dl._nb_parameters with dl.get_nb_parameter when implemented
    test_parameters = dl._nb_parameters

    assert test_parameters[(1, None)]["parameters"] == {'A': 1.0, 'B': 1.0}
    assert test_parameters[(2, None)]["parameters"] == {'A': 4.0, 'B': 4.0}
    assert test_parameters[(1, 2)]["parameters"] == {'A': 2.0, 'B': 2.0}

    dl.add_nb_parameter(atom_type=1, nb_name="Buckingham", nb_model=None, nb_parameters=[1.0, 1.0, 1.0])
    with pytest.raises(KeyError):
        dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model=None, nb_parameters=[1.0, 1.0])


def test_add_nb_parameter_units():
    # Create empty data layer
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(
        atom_type=1,
        nb_name="LJ",
        nb_model="AB",
        nb_parameters=[1.0, 1.0],
        utype=["kJ * mol ** -1 * nanometers ** 12", "kJ * mol ** -1 * nanometers ** 6"])

    # Add AB LJ parameters to data layer - add to two atoms
    dl.add_nb_parameter(
        atom_type=1,
        nb_name="LJ",
        nb_model="AB",
        nb_parameters=[2.0, 2.0],
        atom_type2=2,
        utype=["kJ * mol ** -1 * nanometers ** 12", "kJ * mol ** -1 * nanometers ** 6"])

    # Test atom types out of order
    dl.add_nb_parameter(
        atom_type=3,
        nb_name="LJ",
        nb_model="AB",
        nb_parameters=[2.0, 2.0],
        atom_type2=1,
        utype=["kJ * mol ** -1 * nanometers ** 12", "kJ * mol ** -1 * nanometers ** 6"])

    # Grab stored test parameters
    test_parameters = dl.list_nb_parameters(nb_name="LJ")

    # make sure this grabs all stored values
    assert (set(list(test_parameters)) == set([(1, None), (1, 2), (1, 3)]))

    # Check conversion
    assert dict_compare(test_parameters[(1, None)], {'A': 1.e12, 'B': 1.e6})
    assert dict_compare(test_parameters[(1, 2)], {'A': 2.e12, 'B': 2.e6})
    assert dict_compare(test_parameters[(1, 3)], {'A': 2.e12, 'B': 2.e6})

    # Check keywords
    test_parameters2 = dl.list_nb_parameters(nb_name="LJ", itype="single")
    assert(list(test_parameters2) == [(1, None)])

    test_parameters3 = dl.list_nb_parameters(nb_name="LJ", itype="pair")
    assert (list(test_parameters3) == [(1, 2), (1, 3)])


def test_get_nb_parameter():
    # Create empty data layer
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model="AB", nb_parameters={'A': 1.0, 'B': 2.0})

    # Add Buckingham parameter to datalayer
    dl.add_nb_parameter(atom_type=2, nb_name="Buckingham", nb_parameters={"A": 1.0, "C": 1.0, "rho": 1.0})

    # The following should raise an error because (1,2) interaction is not set
    with pytest.raises(KeyError):
        dl.get_nb_parameter(atom_type=1, atom_type2=2, nb_model="AB")

    # Test that what returned is expected
    assert dict_compare(dl.get_nb_parameter(atom_type=2), {"A": 1.0, "C": 1.0, "rho": 1.0})
    assert dict_compare(dl.get_nb_parameter(atom_type=1, nb_model="AB"), {'A': 1.0, 'B': 2.0})

    # Test conversion of AB to different forms
    assert dict_compare(
        dl.get_nb_parameter(atom_type=1, nb_model="epsilon/sigma"), {'epsilon': 1.0,
                                                                     'sigma': (1. / 2.)**(1. / 6.)})

    assert dict_compare(dl.get_nb_parameter(atom_type=1, nb_model="epsilon/Rmin"), {'epsilon': 1.0, 'Rmin': 1})

    # Test that correct parameters are pulled from data layer based on name
    assert (set(dl.list_stored_nb_types()) == set(["LJ", "Buckingham"]))

    assert dict_compare(dl.list_nb_parameters(nb_name="LJ"), {(1, None): {'A': 1., 'B': 2.}})

    comp = {(2, None): {'A': 1.0, "C": 1.0, "rho": 1.0}}
    assert dict_compare(dl.list_nb_parameters(nb_name="Buckingham"), comp)

    # Test translation of units
    comp = {(1, None): {'A': 1.e-12, 'B': 2.e-6}}
    result = dl.list_nb_parameters(
        nb_name="LJ", utype={'A': "kJ * mol ** -1 * nanometers ** 12",
                             'B': "kJ * mol ** -1 * nanometers ** 6"})
    assert dict_compare(result, comp)

    # Sigma/epsilon test
    comp = {'sigma': (1. / 2.)**(1. / 6.), 'epsilon': eex.units.conversion_factor('kJ', 'kcal')}
    result = dl.get_nb_parameter(
        atom_type=1, nb_model='epsilon/sigma', utype={'epsilon': 'kcal * mol ** -1',
                                                      'sigma': 'angstrom'})
    assert dict_compare(result, comp)


def test_mixing_rule():
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add mixing rule to datalayer
    dl.set_mixing_rule("geometric")

    # Check failure
    with pytest.raises(ValueError):
        dl.set_mixing_rule("test")

    with pytest.raises(TypeError):
        dl.set_mixing_rule(5)

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 1.0, 'sigma': 2.0})

    # Add Buckingham parameter to datalayer
    dl.add_nb_parameter(atom_type=2, nb_name="Buckingham", nb_parameters={"A": 1.0, "C": 1.0, "rho": 1.0})

    # This should fail because we can not combine LJ and Buckingham parameters.
    with pytest.raises(ValueError):
        dl.mix_LJ_parameters(atom_type1=1, atom_type2=2)

    with pytest.raises(KeyError):
        dl.mix_LJ_parameters(atom_type1=1, atom_type2=3)

    # Overwrite Buckingham parameter
    dl.add_nb_parameter(atom_type=2, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 1.0, 'sigma': 1.0})

    # Apply mixing rule
    dl.mix_LJ_parameters(atom_type1=1, atom_type2=2)

    # Get values from datalayer and check
    params = dl.get_nb_parameter(atom_type=1, atom_type2=2, nb_model="epsilon/sigma")

    ans = {'sigma': 2 ** (1. / 2.), 'epsilon': 1.}

    assert dict_compare(params, ans)


def test_mixing_table():
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add mixing rule to datalayer
    dl.set_mixing_rule("arithmetic")

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 1.0, 'sigma': 2.0})

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=2, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 2.0, 'sigma': 1.0})

    # Apply mixing rule
    dl.build_LJ_mixing_table()

    # Get mixed parameters
    pairIJ = dl.list_nb_parameters(nb_name="LJ", itype="pair", nb_model="epsilon/sigma")

    ans = {
        (1, 1): {
            'sigma': 2.,
            'epsilon': 1,
        },

        (1, 2): {
            'sigma': 1.5,
            'epsilon': (2.) ** (1. / 2.)
        },

        (2, 2): {
            'sigma': 1.,
            'epsilon': 2.
        },
    }

    assert(dict_compare(pairIJ, ans))


def test_nb_scaling():
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add mixing rule to datalayer
    dl.set_mixing_rule("arithmetic")

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 1.0, 'sigma': 2.0})

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=2, nb_name="LJ", nb_model="epsilon/sigma", nb_parameters={'epsilon': 2.0, 'sigma': 1.0})

    # Apply mixing rule
    dl.build_LJ_mixing_table()

    # Build scaling dataframe
    scale_df = pd.DataFrame()
    scale_df["coul_scale"] = [0.0, 0.0, 0.0]
    scale_df["atom_index1"] = [1, 1, 2]
    scale_df["atom_index2"] = [2, 3, 3]
    scale_df["vdw_scale"] = [0.5, 0.5, 0.5]

    # Check adding data
    dl.set_pair_scalings(scale_df)

    # Check function failures
    new_df = scale_df.copy()
    new_df["test"] = [0, 0, 0]

    with pytest.raises(KeyError):
        dl.set_pair_scalings(new_df)

    with pytest.raises(ValueError):
        dl.set_pair_scalings(scale_df[["atom_index1", "atom_index2"]])

    scale_df.drop(['atom_index1'], axis=1, inplace=True)
    with pytest.raises(KeyError):
        dl.set_pair_scalings(scale_df)

    # Retrieve information from dl
    with pytest.raises(KeyError):
        dl.get_pair_scalings(nb_labels=["not_a_label"])

    stored_scalings = dl.get_pair_scalings(order=False)

    for col in stored_scalings.columns:
        assert set(scale_df[col].values) == set(stored_scalings[col].values)


def test_set_nb_scaling_factors():
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add bonds to system
    bond_df = pd.DataFrame()
    bond_data = np.array([[0, 1, 0], [1, 2, 0]])
    bond_columns = ["atom1", "atom2", "term_index"]

    for num, name in enumerate(bond_columns):
        bond_df[name] = bond_data[:, num]

    dl.add_bonds(bond_df)

    # Check dl.query_atom_pair
    assert dl.query_atom_pair(1, 2) == 2

    # Add an angle
    angle_df = pd.DataFrame()
    angle_data = np.array([[0, 1, 2, 0]])
    angle_columns = ["atom1", "atom2", "atom3", "term_index"]

    for num, name in enumerate(angle_columns):
        angle_df[name] = angle_data[:, num]

    dl.add_angles(angle_df)

    # Check that build_scaling_list fails now
    with pytest.raises(ValueError):
        dl.build_scaling_list()

    # Check dl.query_atom_pair
    assert dl.query_atom_pair(0, 2) == 3

    scaling_factors = {
        "coul": {
            "scale12": 0.0,
            "scale13": 0.25,
            "scale14": 0.75,
        },

        "vdw": {
            "scale12": 0.0,
            "scale13": 0.5,
            "scale14": 0.75,
        }
    }

    # Test failures
    with pytest.raises(TypeError):
        dl.set_nb_scaling_factors("not a dictionary")

    # Test adding to dl
    dl.set_nb_scaling_factors(scaling_factors)

    # Retrieve from dl
    stored_scaling = dl.get_nb_scaling_factors()

    assert eex.testing.dict_compare(scaling_factors, stored_scaling)

    # Test build_scaling_list
    dl.build_scaling_list()

    # Retrieve from datalayer
    scaling = dl.get_pair_scalings(nb_labels=["vdw_scale", "coul_scale"], order=True)

    assert(set(scaling['vdw_scale'].values) == set([0, 0.5]))

    assert(set(scaling['coul_scale'].values) == set([0, 0.25]))

    assert(set(scaling['order'].values) == set([2, 3]))

def test_calculate_nb_scaling_factors(butane_dl):

    dl = butane_dl(scale=False)

    scaling_df = pd.DataFrame()

    scaling_df['atom_index1'] = [0, 0, 0, 1, 1, 2]
    scaling_df['atom_index2'] = [1, 2, 3, 2, 3, 3]
    scaling_df["vdw_scale"] = [0, 0.5, 0.75, 0, 0.5, 0]
    scaling_df["coul_scale"] = [0.1, 0.25, 0.5, 0.1, 0.25, 0.1]

    dl.set_pair_scalings(scaling_df)

    dl.calculate_nb_scaling_factors()

    scaling_factors_answer = {
        "vdw" :
            {
                "scale12": 0,
                "scale13": 0.5,
                "scale14": 0.75,
            },
        "coul" :
            {
                "scale12": 0.1,
                "scale13": 0.25,
                "scale14": 0.5,
            }
    }

    scaling_factors_calculated = dl.get_nb_scaling_factors()

    eex.testing.dict_compare(scaling_factors_answer, scaling_factors_calculated)

def test_calculate_nb_scaling_factors2(butane_dl):

    dl = butane_dl(scale=False)

    scaling_df = pd.DataFrame()

    # Set some values so they are not compatible (ie all 1-2 interactions aren't the same scaling factor)
    scaling_df['atom_index1'] = [0, 0, 0, 1, 1, 2]
    scaling_df['atom_index2'] = [1, 2, 3, 2, 3, 3]
    scaling_df["vdw_scale"] = [0, 0.5, 0.75, 0, 0.55, 0]
    scaling_df["coul_scale"] = [0.1, 0.25, 0.5, 0.1, 0.25, 0.1]

    dl.set_pair_scalings(scaling_df)

    with pytest.raises(ValueError):
        dl.calculate_nb_scaling_factors()

def test_calculate_nb_scaling_factors3(butane_dl):

    dl = butane_dl(scale=True)

    scaling_df = pd.DataFrame()
    # Add pair scalings that are inconsistent with what is used already
    scaling_df['atom_index1'] = [0, 0, 0, 1, 1, 2]
    scaling_df['atom_index2'] = [1, 2, 3, 2, 3, 3]
    scaling_df["vdw_scale"] = [0, 0.5, 0.75, 0, 0.5, 0]
    scaling_df["coul_scale"] = [0.1, 0.25, 0.5, 0.1, 0.25, 0.1]

    dl.set_pair_scalings(scaling_df)

    with pytest.raises(ValueError):
        dl.calculate_nb_scaling_factors()

def test_remove_terms(butane_dl):

    dl = butane_dl()

    bonds = dl.get_terms(2)

    assert(not bonds.empty)

    dl.remove_terms(2)

    bonds = dl.get_terms(2)

    assert(bonds.empty)
    assert(dl.get_term_count(2)['total'] == 0)

def test_remove_terms_by_index(butane_dl):

    dl = butane_dl()

    bonds1 = dl.get_terms(2)

    assert(not bonds1.empty)

    # Remove first two bonds
    dl.remove_terms(2, index=[0,1])

    bonds2 = dl.get_terms(2)

    assert(dl.get_term_count(2)['total'] == 1)

    # Check that bonds are what we expect
    assert(sorted(bonds1.loc[2].values) == sorted(bonds2.values[0]))

    # Here, since propagate was set to false. Angles and dihedrals remain unchanged.
    assert(dl.get_term_count(3)['total'] == 2)

    assert (dl.get_term_count(4)['total'] == 1)

    return True

def test_remove_terms_by_index_nonconsecutive(butane_dl):

    dl = butane_dl()

    bonds1 = dl.get_terms(2)

    assert(not bonds1.empty)

    # Remove two non-consectutive bonds
    dl.remove_terms(2, index=[0,2])

    bonds2 = dl.get_terms(2)

    assert(dl.get_term_count(2)['total'] == 1)

    # Check that bonds are what we expect
    assert(sorted(bonds1.loc[1].values) == sorted(bonds2.values[0]))

    # Here, since propagate was set to false. Angles and dihedrals remain unchanged.
    assert(dl.get_term_count(3)['total'] == 2)

    assert (dl.get_term_count(4)['total'] == 1)

    return True

def test_remove_terms_propagate(butane_dl):

    dl = butane_dl()

    bonds = dl.get_terms(2)

    # Assert topology is what we expect.
    assert(not bonds.empty)
    assert(dl.get_term_count(2)['total'] == 3)
    assert(dl.get_term_count(3)['total'] == 2)
    assert (dl.get_term_count(4)['total'] == 1)

    dl.remove_terms(2, propagate=True)

    # Assert all have been removed.
    bonds = dl.get_terms(2)

    assert(bonds.empty)
    assert(dl.get_term_count(2)['total'] == 0)
    assert(dl.get_term_count(3)['total'] == 0)
    assert (dl.get_term_count(4)['total'] == 0)

    return True

def test_remove_terms_by_index_propagate(butane_dl):

    dl = butane_dl()

    bonds = dl.get_terms(2)

    # Assert topology is what we expect.
    assert(not bonds.empty)
    assert(dl.get_term_count(2)['total'] == 3)
    assert(dl.get_term_count(3)['total'] == 2)
    assert (dl.get_term_count(4)['total'] == 1)

    assert(not bonds.empty)

    # Remove one bond - choose to propagate this so dihedral and angle should also be removed.
    dl.remove_terms(2, index=[0], propagate=True)

    assert(dl.get_term_count(2)['total'] == 2)
    assert (dl.get_term_count(3)['total'] == 1)
    assert(dl.get_term_count(4)['total'] == 0)

    return True

def test_remove_terms_by_index_nonconsecutive_propagate(butane_dl):
    dl = butane_dl()

    bonds = dl.get_terms(2)

    # Assert topology is what we expect.
    assert (not bonds.empty)

    #print(dl.get_terms(3))

    assert (dl.get_term_count(2)['total'] == 3)
    assert (dl.get_term_count(3)['total'] == 2)
    assert (dl.get_term_count(4)['total'] == 1)

    # Remove two bonds - choose to propagate this so dihedral and both angles should also be removed.
    dl.remove_terms(2, index=[0, 2], propagate=True)

    assert(dl.get_term_count(2)['total'] == 1)
    assert (dl.get_term_count(3)['total'] == 0)
    assert(dl.get_term_count(4)['total'] == 0)

    return True

def test_remove_and_add_terms(butane_dl):
    dl = butane_dl()

    bonds = dl.get_terms(2)

    assert(not bonds.empty)

    dl.remove_terms(2)

    bonds_new = dl.get_terms(2)

    assert(bonds_new.empty)
    assert(dl.get_term_count(2)['total'] == 0)

    dl.add_bonds(bonds)

    added_bonds = dl.get_terms(2)

    assert (not added_bonds.empty)

def test_remove_and_add_terms_index(butane_dl):
    dl = butane_dl()

    bonds = dl.get_terms(2)

    assert(not bonds.empty)

    dl.remove_terms(2, index=[0])

    assert(dl.get_term_count(2)['total'] == 2)

    dl.add_bonds(bonds.iloc[[0]])

    added_bonds = dl.get_terms(2)

    assert (not added_bonds.empty)

    for col in added_bonds.columns:
        assert(set(added_bonds[col].values) == set(bonds[col].values))

def test_remove_term_parameters(butane_dl):
    dl = butane_dl()

    # Check that term parameter cannot be removed if term is not removed first
    with pytest.raises(ValueError):
        dl.remove_term_parameter(order=3, uid=0)

    # First remove angle
    dl.remove_terms(order=3)

    # Test that this produces error (there is only one term parameter (has uid 0))
    with pytest.raises(KeyError):
        dl.remove_term_parameter(order=3, uid=1)

    # Remove angle type from datalayer
    dl.remove_term_parameter(order=3, uid=0)

    # Assert that there are no stored angle parameters
    assert not dl.list_term_parameters(order=3)

    # Test that second removal produces error
    with pytest.raises(KeyError):
        dl.remove_term_parameter(order=3, uid=0)

def test_remove_term_parameters_two(butane_dl):
    dl = butane_dl()

    # Add second angle
    uid = dl.add_term_parameter(3, "harmonic", {'K': 62.100, 'theta0': 116},
                          utype={'K': 'kcal * mol ** -1 * radian ** -2',
                                 'theta0': 'degree'})

    # First remove angle
    dl.remove_terms(order=3)

    # Remove angle type from datalayer
    dl.remove_term_parameter(order=3, uid=0)

    # Add another angle - what will the uid be?
    uid2 = dl.add_term_parameter(3, "harmonic", {'K': 62.100, 'theta0': 117},

    utype = {'K': 'kcal * mol ** -1 * radian ** -2',
             'theta0': 'degree'})

    print("The uids are", uid, uid2)