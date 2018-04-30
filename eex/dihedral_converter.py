"""
Converts various dihedral forms to other equivalents.
"""
import numpy as np


def _alternating_signs(n):
    """
        Return 1 if n is zero
        Returns -1 if n is odd
        Returns 1 if n is even
    """

    if (n < 0.0):
        raise ValueError("Argument has to be zero or a positive integer")

    return (-1.0) ** n


def _cosnx(n):
    """
        Converts Cos(n*x) to a polynomial of cosines whose coefficients are given by a
        Chebyshev polynomial
    """
    p = np.polynomial.chebyshev.Chebyshev.basis(n)
    return p.convert(kind=np.polynomial.Polynomial)


def _OPLS_to_RB(coeffs):
    """
        RB form: "A_0 + A_1 * (cos(phi)) + A_2 * (cos(phi)) ** 2 + A_3 * (cos(phi)) ** 3 + A_4 * (cos(phi)) ** (4) + A_5 * (cos(phi)) ** (5)"

        opls form:"0.5*K_1*(1+cos(phi)) + 0.5 * K_2 * (1-cos(2*phi)) + 0.5 * K_3 * (1+cos(3*phi)) + 0.5 * K_4 * (1-cos(4*phi))"
    """

    k_1 = coeffs['K_1']
    k_2 = coeffs['K_2']
    k_3 = coeffs['K_3']
    k_4 = coeffs['K_4']

    if(k_1 == 0.0 and k_2 == 0.0 and k_3 == 0.0 and k_4 == 0.0):
        raise ValueError('OPLS to RB dihedral conversion not possible. All the coefficients of this dihedral are zero.')

    ret = dict()
    # Note: c1 and c3 are the negative of what is defined on equation 4.64 of Gromacs Manual 4.6.1
    ret['A_0'] = k_2 + 0.5 * (k_1 + k_3)
    ret['A_1'] = 0.5 * (k_1 - 3.0 * k_3)
    ret['A_2'] = -k_2 + 4.0 * k_4
    ret['A_3'] = 2.0 * k_3
    ret['A_4'] = -4.0 * k_4
    ret['A_5'] = 0.0

    return ret


def _RB_to_OPLS(coeffs):

    a_0 = coeffs['A_0']
    a_1 = coeffs['A_1']
    a_2 = coeffs['A_2']
    a_3 = coeffs['A_3']
    a_4 = coeffs['A_4']
    a_5 = coeffs['A_5']

    ret = dict()

    if(a_0 == 0.0 and a_1 == 0.0 and a_2 == 0.0 and a_3 == 0.0 and a_4 == 0.0 and a_5 == 0.0):
        raise ValueError('OPLS to RB dihedral conversion not possible. All the coefficients of this dihedral are zero.')

    if (a_5 != 0.0 and a_1 + a_2 + a_3 + a_4 != 0.0):
        raise ValueError("RB to OPLS dihedral conversion not possible. This RB dihedral is inconsistent with OPLS style")

    # note - f1 and f3 are opposite sign as expected in GROMACS, probably because of angle conventions.
    ret['K_1'] = 2.0 * a_1 + 3.0 * a_3 / 2.0
    ret['K_2'] = -a_2 - a_4
    ret['K_3'] = a_3 / 2.0
    ret['K_4'] = -a_4 / 4.0
    return ret


def _CHARMM_to_RB(coeffs):

    ret = dict()

    n = coeffs['n']
    d = coeffs['d']
    k = coeffs['K']

    if n > 5:
        raise ValueError("CHARMM to RB dihedral conversion not possible. The multiplicity value n in the CHARMM dihedral must be less than 5")

    if n < 0:
        raise ValueError("CHARMM to RB dihedral conversion not possible. The multiplicity must be a positive integer.")
    # Make sure d is either 0.0, pi, 2pi, 3pi, 4pi.
    # TODO: Need to generalize for other values of d.

    if n == 0 and d % np.pi == 0:
        raise ValueError("CHARMM to RB dihedral conversion not possible. The CHARMM dihedral energy is zero.")

    if k == 0:
        raise ValueError("CHARMM to RB dihedral conversion not possible. The CHARMM dihedral energy is zero.")

    if not np.isclose(d % np.pi, 0.0):
        raise ValueError("CHARMM to RB dihedral conversion not possible. The CHARMM phase shift d is not a multiple of PI.")

    div = np.abs(d / np.pi)
    p = list(_alternating_signs(div) * _cosnx(n))
    ret = dict()
    ret["A_0"] = 0.0
    ret["A_1"] = 0.0
    ret["A_2"] = 0.0
    ret["A_3"] = 0.0
    ret["A_4"] = 0.0
    ret["A_5"] = 0.0

    for i in range(0, len(p)):
        ret['A_' + str(i)] = p[i] * k
    ret["A_0"] = ret["A_0"] + k
    return ret


def _RB_to_CHARMM(coeffs):

    a_0 = coeffs['A_0']
    a_1 = coeffs['A_1']
    a_2 = coeffs['A_2']
    a_3 = coeffs['A_3']
    a_4 = coeffs['A_4']
    a_5 = coeffs['A_5']

    # TODO probably we need a more robust way of doing this
    rb_coeffs = np.array([a_0, a_1, a_2, a_3, a_4, a_5])

    if(np.all(rb_coeffs == 0)):
        raise ValueError('All the coefficients of this dihedral are zero.')

    nz = np.nonzero(rb_coeffs)[0]
    n = np.max(nz)
    p = list(_cosnx(n))
    idx = np.nonzero(rb_coeffs)[0][-1]
    if idx == 0:
        k = rb_coeffs[idx] / (p[idx] + 1)
    else:
        k = rb_coeffs[idx] / p[idx]
    ret = dict()
    ret['n'] = n
    ret['K'] = np.abs(k)
    if k < 0.0:
        ret['d'] = np.pi
    else:
        ret['d'] = 0.0

    return ret
