"""
Contains tests for the utility functions within EEX
"""

import eex
import pytest
import numpy as np

@pytest.mark.parametrize("form", ["epsilon/sigma", "epsilon/Rmin", "AB"])
def test_convert_same_LJ_coeffs(form):

    var_names = {
    'AB': ['A', 'B'],
    'epsilon/sigma': ['epsilon', 'sigma'],
    'epsilon/Rmin': ['epsilon', 'Rmin'],
    }
    coeffs = dict.fromkeys(var_names[form])
    coeffs[var_names[form][0]] = 2.0
    coeffs[var_names[form][1]] = 4.0
    new_coeffs = eex.nb_converter.convert_LJ_coeffs(coeffs, form, form)
    assert np.allclose(new_coeffs[var_names[form][0]], coeffs[var_names[form][0]])
    assert np.allclose(new_coeffs[var_names[form][1]], coeffs[var_names[form][1]])

@pytest.mark.parametrize("form", ["epsilon/sigma", "epsilon/Rmin"])
@pytest.mark.parametrize("coeffs", [{'A': 3.48, 'B': 0.0}, {'A': 0.0, 'B': 148.0}])
def test_convert_impossible_LJ_coeffs(form, coeffs):
    with pytest.raises(ZeroDivisionError):
        eex.nb_converter.convert_LJ_coeffs(coeffs, "AB", form)

#@pytest.mark.parametrize("end", ["LJ", "Rmin", "AB"])
#@pytest.mark.parametrize("origin", ["LJ", "Rmin", "AB"])
#def test_convert_LJ_coeffs(origin, end):
#    coeffs = [3.48, 148.0]
#    new_coeffs = eex.nb_converter.convert_LJ_coeffs(coeffs, origin, end)



def test_fuzzy_list_match():

    assert (False, None) == eex.utility.fuzzy_list_match("thing", [])
    assert (True, "data") == eex.utility.fuzzy_list_match("data tmp", ["other", "data", "else"])


def test_find_lowest_hole():

    assert 0 == eex.utility.find_lowest_hole([])
    assert 1 == eex.utility.find_lowest_hole([0])
    assert 2 == eex.utility.find_lowest_hole([0, 1, 3, 4])


def test_hash():

    # Quick hash
    assert '3e88c7ece2cbec43538d8898f8a7601e' == eex.utility.hash([5])
    assert 'a1b8364f4dd47beb4cc0dd364d78d05d' == eex.utility.hash(["something", 5])

    # Nested hash
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, [5, 6]])
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, np.array([5, 6])])
    assert eex.utility.hash([[5], [6]]) == eex.utility.hash(np.array([[5], [6]]))

    # Test numerical hashes
    assert eex.utility.hash([5.0]) == eex.utility.hash([5])
    assert eex.utility.hash([np.int(5)]) == eex.utility.hash([np.float64(5)])
    assert eex.utility.hash([5.0, 6]) == eex.utility.hash([5, 6.0])

    # Test dicts
    assert eex.utility.hash({"k": 5, "d": 6}) == eex.utility.hash({"d": 6.0, "k": 5.0 + 1.e-14})
    assert eex.utility.hash({"k": 5, "d": {"a": 6}}) == eex.utility.hash({"d": {"a": 6.0}, "k": 5.0 + 1.e-14})
    assert eex.utility.hash({"k": 5, "d": {"b": 6}}) != eex.utility.hash({"d": {"a": 6.0}, "k": 5.0 + 1.e-14})

    # Test tolerances
    assert eex.utility.hash(5) == eex.utility.hash(5.0)
    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-10])
    assert eex.utility.hash([5.0, 6.0]) == eex.utility.hash([5.0 + 1.e-10, 6.0 + 1.e-10])

    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-9])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-8])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-7])

    # Change tolerances
    assert eex.utility.hash([5.0], rtol=9) == eex.utility.hash([5.0 + 1.e-10], rtol=9)
    assert eex.utility.hash([5.0], rtol=9) != eex.utility.hash([5.0 + 1.e-9], rtol=9)

    # Fun gotchas
    assert eex.utility.hash([5.0], rtol=8) != eex.utility.hash([5.0], rtol=9)

    assert eex.utility.hash(("k", {"a": 5, "b": 6})) == eex.utility.hash(("k", {"b": 6, "a": 5}))
