"""
Contains tests for the dihedral functional form conversions
"""

import eex
import pytest
import numpy as np
import pytest


def _rb_eval(rb_coeffs):
    phi = np.linspace(-10, 10, 100)
    e_charmm = np.zeros(len(phi))

    for k, v in rb_coeffs.items():
        A_0 = rb_coeffs['A_0']
        A_1 = rb_coeffs['A_1']
        A_2 = rb_coeffs['A_2']
        A_3 = rb_coeffs['A_3']
        A_4 = rb_coeffs['A_4']
        A_5 = rb_coeffs['A_5']
        e_rb = A_0 + A_1 * (np.cos(phi)) + A_2 * (np.cos(phi)) ** 2 + A_3 * (np.cos(phi)) ** 3 + A_4 * (np.cos(phi)) ** (4) + A_5 * (np.cos(phi)) ** (5)
    return e_rb


def _opls_eval(opls_coeffs):
    x = np.linspace(-10, 10, 100)
    e_opls = 0.5 * (opls_coeffs['K_1'] * (1.0 + np.cos(x))
                    + opls_coeffs['K_2'] * (1.0 - np.cos(2 * x))
                    + opls_coeffs['K_3'] * (1.0 + np.cos(3 * x))
                    + opls_coeffs['K_4'] * (1.0 - np.cos(4 * x)))
    return e_opls


def _charmm_eval(charmm_coeffs):
    x = np.linspace(-10, 10, 100)
    e_charmm = np.zeros(len(x))

    K = charmm_coeffs['K']
    n = charmm_coeffs['n']
    d = charmm_coeffs['d']

    for j in range(0, len(n)):
        e_charmm += K[j] * (1 + np.cos(x * n[j] - d[j]))
    return e_charmm


# @pytest.mark.parametrize("coeffs", [(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#                                     ])
# def test_rb_to_charmm_exceptions(coeffs):

#     rb_coeffs = dict()

#     for idx, c in enumerate(coeffs):
#         rb_coeffs['A_' + str(idx)] = c

#     with pytest.raises(ValueError):
#         charmm_coeffs = eex.form_converters.convert_form(4, rb_coeffs, 'RB', 'charmmfsw')


# @pytest.mark.parametrize("knd", [(0.0, 1, np.pi),
#                                  (4.0, 0, np.pi),
#                                  (4.0, 0, -np.pi),
#                                  (4.0, 6, 0.0)])
# def test_charmm_to_rb_exceptions(knd):

#     K, n, d = knd
#     charmm_coeffs = dict()
#     charmm_coeffs[n] = {'K': K, 'n': n, 'd': d}
#     with pytest.raises(ValueError):
#         rb_coeffs = eex.form_converters.convert_form(4, charmm_coeffs, 'charmmfsw', 'RB')


# @pytest.mark.parametrize("K", [-4.0])
# @pytest.mark.parametrize("n", [1])
# @pytest.mark.parametrize("d", [0.0])
# def test_charmm_to_rb_loop(K, n, d):
#     charmm_coeffs = dict()
#     charmm_coeffs[1] = {'K': K, 'n': n, 'd': d}
#     # print(charmm_coeffs)
#     rb_coeffs = eex.form_converters.convert_form(4, charmm_coeffs, 'charmmfsw', 'RB')
#     e_charmm1 = _charmm_eval(charmm_coeffs)
#     e_rb = _rb_eval(rb_coeffs)
#     # assert np.allclose(e_charmm1, e_rb)
#     charmm_coeffs = eex.form_converters.convert_form(4, rb_coeffs, 'RB', 'charmmfsw')
#     # print(charmm_coeffs)
#     # e_charmm2 = _charmm_eval(charmm_coeffs)
#     # assert np.allclose(e_charmm2, e_rb)
#     # assert np.allclose(e_charmm1, e_charmm2)


# def test_opls_to_rb_loop():

#     const = 10.0 * (2 * np.random.random(4) - 1.0)
#     opls_coeffs = {'K_1': const[0], 'K_2': const[1], 'K_3': const[2], 'K_4': const[3]}
#     rb_coeffs = eex.form_converters.convert_form(4, opls_coeffs, 'opls', 'RB')
#     e_opls = _opls_eval(opls_coeffs)
#     e_rb = _rb_eval(rb_coeffs)
#     assert np.allclose(e_opls, e_rb)
#     opls_coeffs = eex.form_converters.convert_form(4, rb_coeffs, 'RB', 'opls')
#     e_opls = _opls_eval(opls_coeffs)
#     assert np.allclose(e_opls, e_rb)


# def _generate_params(form_type):
#     if form_type == 'opls':
#         const = 10.0 * (2 * np.random.random(4) - 1.0)
#         coeffs = {'K_1': const[0], 'K_2': const[1], 'K_3': const[2], 'K_4': const[3]}
#     if form_type == 'charmmfsw':
#         n_orders = np.random.randint(5)
#         coeffs = dict()
#         for k in range(1, n_orders):
#             coeffs[k] = dict()
#             coeffs[k]['n'] = k
#             coeffs[k]['K'] = 10 * (2 * np.random.random() - 1.0)
#             coeffs[k]['d'] = np.round(np.random.random()) * np.pi
#     if form_type == 'RB':
#         coeffs = dict()
#         for idx in range(0, 6):
#             coeffs['A_' + str(idx)] = 10.0 * (2 * np.random.random() - 1.0)
#         coeffs['A_5'] = 0.0
#     return coeffs


# @pytest.mark.parametrize("type_1", ['opls', 'charmmfsw', 'RB'])
# @pytest.mark.parametrize("type_2", ['opls', 'charmmfsw', 'RB'])
def test_quick():

    # type_1_coeffs = _generate_params(type_1)
    # type_2_coeffs = eex.form_converters.convert_form(4, type_1_coeffs, type_1, type_2)

    opls_coeffs = {'K_1': 5.9037668417599996, 'K_2': -1.1339262997600001, 'K_3': 13.15880003896, 'K_4': 0.0}
    # {'K': array([-0.56696315,  2.95188342,  0.56696315,  6.57940002]), 'n': array([ 0.,  1.,  2.,  3.]), 'd': array([ 0.,  0.,  0.,  0.])}

    charmm_coeffs = eex.form_converters.convert_form(4, opls_coeffs, 'opls', 'charmmfsw')

    # print(rb_coeffs)
    e_opls = _opls_eval(opls_coeffs)
    # e_rb = _rb_eval(rb_coeffs)
    # print(e_rb[0])
    # print(e_opls[0])

    print(charmm_coeffs)
    e_charmm = _charmm_eval(charmm_coeffs)
    print(e_opls[0] - e_charmm[0])

# assert np.allclose(e_opls, e_rb)
# charmm_coeffs = eex.form_converters.convert_form(4, rb_coeffs, 'RB', 'charmmfsw')
# e_charmm = _charmm_eval(charmm_coeffs)
# assert np.allclose(e_rb, e_charmm)


# @pytest.mark.parametrize("coeffs", [(0.0, 0.0, 0.0, 0.0),
#                                     ])
# def test_opls_to_rb_exceptions(coeffs):

#     opls_coeffs = dict()

#     for idx, c in enumerate(coeffs):
#         opls_coeffs['K_' + str(idx + 1)] = c

#     with pytest.raises(ValueError):
#         rb_coeffs = eex.form_converters.convert_form(4, opls_coeffs, 'opls', 'RB')


# @pytest.mark.parametrize("coeffs", [(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#                                     ])
# def test_rb_to_opls_exceptions(coeffs):

#     rb_coeffs = dict()

#     for idx, c in enumerate(coeffs):
#         rb_coeffs['A_' + str(idx)] = c

#     with pytest.raises(ValueError):
#         opls_coeffs = eex.form_converters.convert_form(4, rb_coeffs, 'RB', 'opls')
