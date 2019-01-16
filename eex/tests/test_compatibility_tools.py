"""

Tests datalayer compatbility functions

"""
import eex
import pytest

from eex.translators.amber.amber_ff import term_data as amber_term_data

""" What tests do we need?

1. Passes when everything is compatible
2. Fails when not compatible
    - Bonds
    - Angles
    - Dihedrals
3. Performs conversion when applicable
    - Dihedrals
    
------------
NB Tests (different format that n body terms)
* Metadata format must be finalized. 

1. Fails with incompatible NB type
2. Returns correct - pair or single 


"""

def test_compatibility_no_changes(butane_dl):

    # Initialize datalayer
    dl = butane_dl()

    # Grab details from datalayer
    term_3_type, term_3 = dl.get_term_parameter(order=3, uid=0)
    term_2_type, term_2 = dl.get_term_parameter(order=2, uid=0)

    # Run compatibility check function
    eex.metadata.check_term_compatibility(dl, amber_term_data)

    # Grab new details from datalayer
    term_3_new_type, term_3_new = dl.get_term_parameter(order=3, uid=0)
    term_2_new_type, term_2_new = dl.get_term_parameter(order=2, uid=0)

    # Assert that keywords are the same.
    assert term_2_type == term_2_new_type
    assert term_3_type == term_3_new_type

    # Assert that parameters are unchanged
    eex.testing.dict_compare(term_3_new, term_3)
    eex.testing.dict_compare(term_2_new, term_2)

def test_compatibility_term_2(butane_dl):

    dl = butane_dl(ff=False)

    # Add bond type which is incompatible with amber
    dl.add_term_parameter(2, "fene", {'K': 300.9, 'R0': 1.540, 'epsilon': 0, 'sigma': 0}, uid=0)

    # Run compatibility check
    with pytest.raises(TypeError):
        eex.metadata.check_term_compatibility(dl, amber_term_data)


def test_compatibility_term_3(butane_dl):

    dl = butane_dl(ff=False)

    # Add angle type which is incompatible with amber
    dl.add_term_parameter(3, "cosine", {'K': 62.100}, uid=0)

    # Run compatibility check
    with pytest.raises(TypeError):
        eex.metadata.check_term_compatibility(dl, amber_term_data)

def test_compatibility_term_4(butane_dl):

    dl = butane_dl(ff=False)

    # Add term parameter which is incompatible with amber
    dl.add_term_parameter(4, "quadratic", {'K': 1, 'phi0': 0}, uid=0)

    # Run compatibility check
    with pytest.raises(TypeError):
        eex.metadata.check_term_compatibility(dl, amber_term_data)


def test_compatibility_conversion(butane_dl):

    dl = butane_dl()

    # Add term parameter which must be converted to be compatible with amber (opls)
    dl.add_term_parameter(4, "opls", {'K_1': 1.41103414, 'K_2': -0.27101489,
                                      'K_3': 3.14502869, 'K_4': 0}, uid=0, utype={'K_1': 'kcal * mol ** -1',
                                                                                  'K_2': 'kcal * mol ** -1',
                                                                                  'K_3': 'kcal * mol ** -1',
                                                                                  'K_4': 'kcal * mol ** -1'})

    # Run compatibility check
    eex.metadata.check_term_compatibility(dl, amber_term_data)

    term_4_new = dl.list_term_parameters(4)

    for parameter in term_4_new.values():
        assert parameter[0] == list(amber_term_data[4].keys())[0]


"""
def test_amber_compatibility_NB_number(butane_dl):

    dl = butane_dl(nb=False)

    oname = eex_find_files.get_scratch_directory("dl_compatibility.prmtop")

    dl.add_nb_parameter(atom_type=1, atom_type2=2, nb_name="LJ",
                        nb_model="epsilon/sigma", nb_parameters={'sigma': 3.75, 'epsilon': 0.1947460018},
                        utype={'sigma': 'angstrom', 'epsilon': 'kcal * mol ** -1'})

    with pytest.raises(ValueError):
        eex.translators.amber.write_amber_file(dl, oname)

def test_amber_compatibility_NB_type(butane_dl):

    dl = butane_dl(nb=False)

    dl.add_nb_parameter(atom_type=1, nb_name="Buckingham", nb_model=None, nb_parameters=[1.0, 1.0, 1.0])

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
        dl.add_term_parameter(2, "fene", {"K": 1, "R0": 1, "epsilon": 1, "sigma": 1,})
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
"""
