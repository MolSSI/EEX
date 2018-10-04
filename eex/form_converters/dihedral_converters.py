"""
Converts various dihedral forms to other equivalents.
"""
import numpy as np
from .converter_registry import register_converter
from scipy.special import comb


def _alternating_signs(n):
    """
        Return 1 if n is zero
        Returns -1 if n is odd
        Returns 1 if n is even
    """

    if (n < 0.0):
        raise ValueError("Argument has to be zero or a positive integer")

    return (-1.0)**n


def _cosnx(n):
    """
        Converts Cos(n*x) to a polynomial of cosines whose coefficients are given by a
        Chebyshev polynomial
    """

    p = np.polynomial.chebyshev.Chebyshev.basis(n)
    return p.convert(kind=np.polynomial.Polynomial)


def _cosx_n(n):
    ret = np.zeros(n + 1)
    if (n % 2 == 0):
        for k in range(0, n // 2):
            idx = n - 2 * k
            ret[idx] = comb(n, k)
        ret = ret * 2 / (2**n)
        ret[0] += 1 / (2**n) * comb(n, n // 2)
    else:

        for k in range(0, (n - 1) // 2 + 1):
            idx = n - 2 * k
            ret[idx] = comb(n, k)
        ret = ret * 2 / (2**n)

    return ret


@register_converter(order=4)
def _RB_to_RB(coeffs):
    if not isinstance(coeffs, dict):
        raise TypeError(
            "RB to RB dihedral conversion requires an input dictionary")

    return coeffs


@register_converter(order=4)
def _opls_to_RB(coeffs):
    """
        RB form: "A_0 + A_1 * (cos(phi)) + A_2 * (cos(phi)) ** 2 + A_3 * (cos(phi)) ** 3 + A_4 * (cos(phi)) ** (4) + A_5 * (cos(phi)) ** (5)"

        opls form:"0.5*K_1*(1+cos(phi)) + 0.5 * K_2 * (1-cos(2*phi)) + 0.5 * K_3 * (1+cos(3*phi)) + 0.5 * K_4 * (1-cos(4*phi))"
    """

    if not isinstance(coeffs, dict):
        raise TypeError(
            "OPLS to RB dihedral conversion requires an input dictionary")

    k_1 = coeffs['K_1']
    k_2 = coeffs['K_2']
    k_3 = coeffs['K_3']
    k_4 = coeffs['K_4']

    if (k_1 == 0.0 and k_2 == 0.0 and k_3 == 0.0 and k_4 == 0.0):
        raise ValueError(
            'OPLS to RB dihedral conversion not possible. All the coefficients of this dihedral are zero.'
        )

    ret = dict()
    # Note: c1 and c3 are the negative of what is defined on equation 4.64 of Gromacs Manual 4.6.1
    ret['A_0'] = k_2 + 0.5 * (k_1 + k_3)
    ret['A_1'] = 0.5 * (k_1 - 3.0 * k_3)
    ret['A_2'] = -k_2 + 4.0 * k_4
    ret['A_3'] = 2.0 * k_3
    ret['A_4'] = -4.0 * k_4
    ret['A_5'] = 0.0

    return ret


@register_converter(order=4)
def _RB_to_opls(coeffs):

    if not isinstance(coeffs, dict):
        raise TypeError(
            "RB to OPLS dihedral conversion requires an input dictionary")

    a_0 = coeffs['A_0']
    a_1 = coeffs['A_1']
    a_2 = coeffs['A_2']
    a_3 = coeffs['A_3']
    a_4 = coeffs['A_4']
    a_5 = coeffs['A_5']

    ret = dict()

    if (a_0 == 0.0 and a_1 == 0.0 and a_2 == 0.0 and a_3 == 0.0 and a_4 == 0.0
            and a_5 == 0.0):
        raise ValueError(
            'OPLS to RB dihedral conversion not possible. All the coefficients of this dihedral are zero.'
        )

    if (a_5 != 0.0 and a_1 + a_2 + a_3 + a_4 != 0.0):
        raise ValueError(
            "RB to OPLS dihedral conversion not possible. This RB dihedral is inconsistent with OPLS style"
        )

    # note - f1 and f3 are opposite sign as expected in GROMACS, probably because of angle conventions.
    ret['K_1'] = 2.0 * a_1 + 3.0 * a_3 / 2.0
    ret['K_2'] = -a_2 - a_4
    ret['K_3'] = a_3 / 2.0
    ret['K_4'] = -a_4 / 4.0
    return ret


@register_converter(order=4)
def _charmmfsw_to_RB(coeffs):

    if not isinstance(coeffs, dict):
        raise TypeError(
            "CHARMM to RB dihedral conversion requires an input dictionary")

    tmp = np.zeros(6)
    n = coeffs['n']
    d = coeffs['d']
    k = coeffs['K']

    if n > 5:
        raise ValueError(
            "CHARMM to RB dihedral conversion not possible. The multiplicity value n in the CHARMM dihedral must be less than 5"
        )

    if n < 0:
        raise ValueError(
            "CHARMM to RB dihedral conversion not possible. The multiplicity must be a positive integer."
        )
    # Make sure d is either 0.0, pi, 2pi, 3pi, 4pi.
    # TODO: Need to generalize for other values of d.

    if n == 0 and d % np.pi == 0:
        raise ValueError(
            "CHARMM to RB dihedral conversion not possible. The CHARMM dihedral energy is zero."
        )

    if k == 0:
        raise ValueError(
            "CHARMM to RB dihedral conversion not possible. The CHARMM dihedral energy is zero."
        )

    if not np.isclose(d % np.pi, 0.0):
        raise ValueError(
            "CHARMM to RB dihedral conversion not possible. The CHARMM phase shift d is not a multiple of PI."
        )

    div = np.abs(d / np.pi)
    p = list(_alternating_signs(div) * _cosnx(n))

    for i in range(0, len(p)):
        tmp[i] = p[i] * k
    tmp[0] += k

    ret = dict()

    for idx in range(0, len(tmp)):
        ret['A_' + str(idx)] = tmp[idx]
    return ret


@register_converter(order=4)
def _RB_to_charmmfsw(coeffs):
    # TODO Need more efficient way to do this conversion

    if not isinstance(coeffs, dict):
        raise TypeError(
            "RB to CHARMM dihedral conversion requires an input dictionary")

    a_0 = coeffs['A_0']
    a_1 = coeffs['A_1']
    a_2 = coeffs['A_2']
    a_3 = coeffs['A_3']
    a_4 = coeffs['A_4']
    a_5 = coeffs['A_5']

    # TODO probably we need a more efficient way of doing this
    rb_coeffs = np.array([a_0, a_1, a_2, a_3, a_4, a_5])

    if (np.all(rb_coeffs == 0)):
        raise ValueError('All the coefficients of this dihedral are zero.')

    nz = np.nonzero(rb_coeffs)
    charmm_coeffs = np.zeros(np.size(rb_coeffs))
    for order in nz[0]:
        tmp = rb_coeffs[order] * _cosx_n(order)
        tmp.resize(charmm_coeffs.shape)
        charmm_coeffs += tmp

    nz_max = np.max(nz)
    ret = dict()
    ret['K'] = np.zeros(nz_max + 1)
    ret['n'] = np.zeros(nz_max + 1)
    ret['d'] = np.zeros(nz_max + 1)

    flag = False
    if charmm_coeffs[0] < 0:
        flag = True
        charmm_coeffs = -charmm_coeffs
        ret['K'][0] = -charmm_coeffs[0] + np.sum(np.abs(charmm_coeffs[1:]))

    else:
        ret['K'][0] = charmm_coeffs[0] - np.sum(np.abs(charmm_coeffs[1:]))

    ret['K'][0] = 0.5 * ret['K'][0]

    for order, coeff in enumerate(charmm_coeffs):
        if order is 0:
            continue

        if np.isclose(coeff, 0.0):
            continue
        ret['n'][order] = order

        if flag is True:
            ret['K'][order] = -charmm_coeffs[order]
        else:
            ret['K'][order] = charmm_coeffs[order]

        if coeff < 0:
            if flag is True:
                ret['K'][order] = charmm_coeffs[order]
            else:
                ret['K'][order] = np.abs(charmm_coeffs[order])
            ret['d'][order] = np.pi
        else:
            ret['d'][order] = 0.0
    return ret
