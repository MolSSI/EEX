"""
Contains tests for the NB functions of EEX
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
    coeffs = {var_names[form][0]: 2.0, var_names[form][1]: 4.0}
    new_coeffs = eex.nb_converter.convert_LJ_coeffs(coeffs, form, form)
    assert np.allclose(new_coeffs[var_names[form][0]], coeffs[var_names[form][0]])
    assert np.allclose(new_coeffs[var_names[form][1]], coeffs[var_names[form][1]])


@pytest.mark.parametrize("form", ["epsilon/sigma", "epsilon/Rmin"])
@pytest.mark.parametrize("coeffs", [{'A': 3.48, 'B': 0.0}, {'A': 0.0, 'B': 148.0}])
def test_convert_impossible_LJ_coeffs(form, coeffs):
    with pytest.raises(ZeroDivisionError):
        eex.nb_converter.convert_LJ_coeffs(coeffs, "AB", form)


## Test LJ mixing

@pytest.mark.parametrize("mixing_rule", ["lorentz-berthelot", "geometric", "sixthpower"])
def test_LJ_mixing(mixing_rule):
    coeff1 = {'epsilon': 1, 'sigma': 1}
    coeff2 = {'epsilon': 2, 'sigma': 2}


    # Expected answers
    mixed_coeffs = {"lorentz-berthelot": {'epsilon': (2) ** (1./2.), 'sigma': 1.5},
                    "geometric": {'epsilon': (2) ** (1./2.), 'sigma': (2) ** (1./2.)},
                    "sixthpower": {'epsilon': (2 ** (1./2.)) / (4), 'sigma': 32. **(1./6.)},
                    }



    new_coeffs = eex.nb_converter.mix_LJ(coeff1, coeff2, origin="epsilon/sigma", mixing_rule=mixing_rule, final='epsilon/sigma')
    print(new_coeffs)

    assert(eex.testing.dict_compare(new_coeffs, mixed_coeffs[mixing_rule]))

