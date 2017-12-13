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
    assert eex.testing.df_compare(tmp_df, dl_df)

    # Compare rand DF
    dl_rand_df = dl.get_other("atoms")
    assert eex.testing.df_compare(rand_df, dl_rand_df)


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

    # Multi-unit terms
    assert 0 == dl.add_atom_parameter("lj_coeffs", {"LJ A": 5.0, "LJ B": 6.0})
    assert 0 == dl.add_atom_parameter("lj_coeffs", {"LJ B": 6.0, "LJ A": 5.0})
    assert 1 == dl.add_atom_parameter("lj_coeffs", {"LJ B": 7.0, "LJ A": 5.0})
    assert 1 == dl.add_atom_parameter("lj_coeffs", {"LJ B": 7.0, "LJ A": 5.0}, uid=1)

    # Test possible raises
    with pytest.raises(KeyError):
        dl.add_atom_parameter("mass", 6.0, uid=0)

    with pytest.raises(TypeError):
        assert 0 == dl.add_atom_parameter("lj_coeffs", 5.0)

    with pytest.raises(KeyError):
        dl.add_atom_parameter("lj_coeffs", {"LJA ": 5.0, "LJB": 6.0})

    with pytest.raises(KeyError):
        dl.add_atom_parameter("lj_coeffs", {"LJ B": 7.0, "LJ A": 5.0}, uid=2)


def test_add_atom_parameter_units():
    dl = eex.datalayer.DataLayer("test_add_atom_parameters")

    # Test duplicate and uid add of same
    assert 0 == dl.add_atom_parameter("mass", 5.0, utype="kilogram / mol")
    assert pytest.approx(5.0) == dl.get_atom_parameter("mass", 0, utype="kilogram / mol")
    assert pytest.approx(5000.0) == dl.get_atom_parameter("mass", 0, utype="gram / mol")
    assert pytest.approx(5000.0) == dl.get_atom_parameter("mass", 0, utype={"mass": "gram / mol"})

    # Test multi-unit utypes
    assert 0 == dl.add_atom_parameter(
        "lj_coeffs", {"LJ A": 5.0,
                      "LJ B": 6.0},
        utype={"LJ A": "kcal/mol * angstrom ** -12",
               "LJ B": "kcal/mol * angstrom ** -6"})

    assert 0 == dl.add_atom_parameter(
        "lj_coeffs", {"LJ A": 10.0,
                      "LJ B": 12.0},
        utype={"LJ B": "0.5 * kcal/mol * angstrom ** -6",
               "LJ A": "0.5 * kcal/mol * angstrom ** -12"})


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
    assert eex.testing.dict_compare(parm1[1], {"K": 4.0, "R0": 5.0})
    assert eex.testing.dict_compare(parm1[1], {"K": 4.0, "R0": 5.0 + 1.e-12})

    parm2 = dl.get_term_parameter(2, 1)
    assert parm2[0] == "harmonic"
    assert eex.testing.dict_compare(parm2[1], {"K": 6.0, "R0": 7.0})

    with pytest.raises(KeyError):
        dl.get_term_parameter(2, 1231234123)


def test_get_term_parameters_units():
    """
    Test obtaining parameters from the DL with units
    """

    dl = eex.datalayer.DataLayer("test_get_parameters_units")

    utype_2b = {"K": "(kJ / mol) * angstrom ** -2", "R0": "angstrom"}
    assert 0 == dl.add_term_parameter(2, "harmonic", {"K": 4.0, "R0": 5.0}, utype=utype_2b)
    assert 1 == dl.add_term_parameter(2, "harmonic", {"K": 6.0, "R0": 7.0}, utype=utype_2b)

    utype_2scale = {"K": "2.0 * (kJ / mol) * angstrom ** -2", "R0": "2.0 * angstrom"}
    parm1 = dl.get_term_parameter(2, 0, utype=utype_2scale)
    assert parm1[0] == "harmonic"
    assert eex.testing.dict_compare(parm1[1], {"K": 2.0, "R0": 2.5})

    parm2 = dl.get_term_parameter(2, 0, utype=[utype_2scale["K"], utype_2scale["R0"]])
    assert parm2[0] == "harmonic"
    assert eex.testing.dict_compare(parm2[1], {"K": 2.0, "R0": 2.5})

    with pytest.raises(TypeError):
        dl.get_term_parameter(2, 0, utype={5, 6})


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


def test_box_size():
    dl = eex.datalayer.DataLayer("test_box_size", backend="memory")

    tmp = {"x": [-5, 5], "y": [-6, 6], "z": [-7, 7]}

    # Normal set/get
    dl.set_box_size(tmp)
    comp = dl.get_box_size()
    eex.testing.dict_compare(tmp, comp)

    # Set/get with units
    dl.set_box_size(tmp, utype="miles")
    comp = dl.get_box_size(utype="miles")
    eex.testing.dict_compare(tmp, comp)

    with pytest.raises(AssertionError):
        dl.set_box_size(tmp, utype="miles")
        comp = dl.get_box_size()
        eex.testing.dict_compare(tmp, comp)


def test_add_nb_parameter():

    # Create empty data layer
    dl = eex.datalayer.DataLayer("test_add_nb_parameters", backend="memory")

    # Create system with three molecules
    atom_sys = _build_atom_df(3)

    # Add atomic system to datalayer
    dl.add_atoms(atom_sys)

    # Add AB LJ parameters to data layer - add to single atom
    dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_form="AB", nb_parameters=[1.0, 1.0])
    dl.add_nb_parameter(atom_type=2, nb_name="LJ", nb_form="epsilon/sigma", nb_parameters=[1.0, 1.0])

    # Add AB LJ parameters to data layer - add to two atoms
    dl.add_nb_parameter(atom_type=1, atom_type2=2, nb_name="LJ", nb_form="AB", nb_parameters=[2.0, 2.0])

    # Grab stored test parameters - will need to replace dl._nb_parameters with dl.get_nb_parameter when implemented
    test_parameters = dl._nb_parameters

    assert test_parameters[(1, )]["parameters"] == {'A': 1.0, 'B': 1.0}
    assert test_parameters[(2, )]["parameters"] == {'A': 4.0, 'B': 4.0}
    assert test_parameters[(1, 2)]["parameters"] == {'A': 2.0, 'B': 2.0}

    dl.add_nb_parameter(atom_type=1, nb_name="Buckingham", nb_form=None, nb_parameters=[1.0, 1.0, 1.0])
    with pytest.raises(KeyError):
        dl.add_nb_parameter(atom_type=1, nb_name="LJ", nb_form=None, nb_parameters=[1.0, 1.0])



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
        nb_form="AB",
        nb_parameters=[1.0, 1.0],
        utype=["kJ * mol ** -1 * nanometers ** 12", "kJ * mol ** -1 * nanometers ** 6"])

    # Add AB LJ parameters to data layer - add to two atoms
    dl.add_nb_parameter(
        atom_type=1,
        nb_name="LJ",
        nb_form="AB",
        nb_parameters=[2.0, 2.0],
        atom_type2=2,
        utype=["kJ * mol ** -1 * nanometers ** 12", "kJ * mol ** -1 * nanometers ** 6"])

    # Grab stored test parameters - will need to replace dl._nb_parameters with dl.get_nb_parameter when implemented
    test_parameters = dl._nb_parameters

    # Check conversion
    eex.testing.dict_compare(test_parameters[(1, )]['parameters'], {'A': 1.e12, 'B': 1.e6})
    eex.testing.dict_compare(test_parameters[(1, 2)]['parameters'], {'A': 2.e12, 'B': 2.e6})

def test_get_nb_parameter():
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
        nb_form="AB",
        nb_parameters={'A': 1.0, 'B': 2.0})

    # Add Buckingham parameter to datalayer
    dl.add_nb_parameter(atom_type=2, nb_name="Buckingham", nb_parameters={"A": 1.0, "C": 1.0, "rho": 1.0})

    # The following should raise an error because (1,2) interaction is not set
    with pytest.raises(KeyError):
        dl.get_nb_parameter(atom_type=1, atom_type2=2, nb_form="AB")

    # The following should raise an error because nb_form is not set for LJ interaction
    with pytest.raises(KeyError):
        dl.get_nb_parameter(atom_type=1)

    # The following should raise an error because units are set for A, but not for B
    #with pytest.raises(KeyError):
    #    dl.get_nb_parameter(atom_type=1, nb_form="AB", utype={'A': 'kJ * mol ** -1 * nanometers ** 12'})

    # Test that what returned is expected
    assert(dl.get_nb_parameter(atom_type=2) == {"A": 1.0, "C": 1.0, "rho": 1.0} )
    assert(dl.get_nb_parameter(atom_type=1, nb_form="AB") == {'A': 1.0, 'B': 2.0} )

    # Test conversion of AB to different forms
    assert(dl.get_nb_parameter(atom_type=1, nb_form="epsilon/sigma")  == {'epsilon': 1.0,'sigma': (1./2.) ** (1./6.)})

    assert (dl.get_nb_parameter(atom_type=1, nb_form="epsilon/Rmin") == {'epsilon': 1.0, 'Rmin': 1 })


    assert(set(dl.list_stored_nb_types()) == set(["LJ", "Buckingham"]))
    print(dl.list_stored_nb_types())
    dl.list_nb_parameters()
    assert False

