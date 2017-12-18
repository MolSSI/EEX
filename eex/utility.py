"""
A EEX utility folder
"""

import os
import hashlib

import numpy as np

## LJ Conversions
def _ab_to_ab(coeffs):
    """
    Convert AB representation to AB representation of the LJ potential
    """
    return {'A': coeffs['A'], 'B': coeffs['B']}

def _epsilonsigma_to_ab(coeffs):
    """
    Convert epsilon/sigma representation to AB representation of the LJ
    potential
    """
    A = 4.0 * coeffs['epsilon'] * coeffs['sigma'] ** 12.0
    B = 4.0 * coeffs['epsilon'] * coeffs['sigma'] ** 6.0
    return {"A": A, "B": B}

def _ab_to_epsilonsigma(coeffs):
    """
    Convert AB representation to epsilon/sigma representation of the LJ
    potential
    """
    if (coeffs['A'] == 0.0 and coeffs['B'] == 0.0):
        return {"sigma": 0.0, "epsilon": 0.0}

    try:
        sigma = (coeffs['A'] / coeffs['B']) ** (1.0 / 6.0)
        epsilon = coeffs['B'] ** 2.0 / (4.0 * coeffs['A'])

    except ZeroDivisionError:
        raise ZeroDivisionError("Lennard Jones functional form conversion not possible, division by zero found.")

    return {"sigma": sigma, "epsilon": epsilon}

def _rminepsilon_to_ab(coeffs):
    """
    Convert rmin/epsilon representation to AB representation of the LJ
    potential
    """
    A = coeffs['epsilon'] * coeffs['Rmin'] ** 12.0
    B = 2 * coeffs['epsilon'] * coeffs['Rmin'] ** 6.0
    return {"A": A, "B": B}

def _ab_to_rminepsilon(coeffs):
    """
    Convert AB representation to Rmin/epsilon representation of the LJ potential
    """

    if (coeffs['A'] == 0.0 and coeffs['B'] == 0.0):
        return {"sigma": 0.0, "epsilon": 0.0}

    try:
        Rmin = (2.0 * coeffs['A'] / coeffs['B'])**(1.0 / 6.0)
        Eps = coeffs['B']**2.0 / (4.0 * coeffs['A'])

    except ZeroDivisionError:
        raise ZeroDivisionError("Lennard Jones functional form conversion not possible, division by zero found.")

    return {"Rmin": Rmin, "epsilon": Eps}


_conversion_matrix = {
    'AB': (['A', 'B'], _ab_to_ab, _ab_to_ab),
    'epsilon/sigma': (['epsilon', 'sigma'], _epsilonsigma_to_ab, _ab_to_epsilonsigma),
    'epsilon/Rmin': (['epsilon', 'Rmin'], _rminepsilon_to_ab, _ab_to_rminepsilon),
}


def convert_LJ_coeffs(coeffs, origin, final):

    difference = set([origin, final]) - set(_conversion_matrix.keys())
    if (difference):
        raise KeyError("Conversion cannot be made since %s is not in conversion matrix %s" %
                       (difference, _conversion_matrix.keys()))

    difference = set(coeffs.keys()) - set(_conversion_matrix[origin][0])
    if (difference):
        raise KeyError("The key %s in the coefficient dictionary is not in the list of allowed keys %s" %
                       (difference, _conversion_matrix[origin][0]))

    internal = _conversion_matrix[origin][1](coeffs)
    external = _conversion_matrix[final][2](internal)
    return external

## LJ Combining Rules
def _mix_LJ(coeff_i, coeff_j, origin, final):
    # Calculate interactions between two atom types based on Lorentz-Berthelot mixing rules

    # First convert from input form to internal AB representation
    internal_coeff_i = convert_LJ_coeffs(coeff_i, origin=origin, final="AB" )
    internal_coeff_j = convert_LJ_coeffs(coeff_j, origin=origin, final="AB" )

    # Convert from internal AB representation to epsilon/sigma
    sigma_epsilon_i = convert_LJ_coeffs(internal_coeff_i, origin="AB", final="epsilon/sigma")
    sigma_epsilon_j = convert_LJ_coeffs(internal_coeff_j, origin="AB", final="epsilon/sigma")

    # Calculate new parameters based on mixing rules
    ## Call functions here

    # Convert from epsilon-sigma to AB, then to final specified form. Double conversion is necessary because of
    # form of conversion matrix.
    convert_params_temp = convert_LJ_coeffs(new_params, origin="epsilon/sigma", final=final)
    convert_params = _conversion_matrix[final][2](convert_params_temp)

    return convert_params

def _lorentz_berthelot(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = (sigma_epsilon_i['sigma'] + sigma_epsilon_j['sigma']) / 2
    new_params['epsilon'] = (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['sigma']) ** (1. / 2.)

    return new_params

def _geometric(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = ( sigma_epsilon_i['sigma'] * sigma_epsilon_j['sigma'] ) ** (1./2.)
    new_params['epsilon'] = (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['epsilon']) ** (1. / 2.)

    return new_params

def _sixthpower(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = ((sigma_epsilon_i['sigma'] ** 6 + sigma_epsilon_j['sigma'] ** 6) / 2 ) ** (1./6.)
    new_params['epsilon'] = ( 2 * (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['epsilon']) ** (1./2.) *
                              sigma_epsilon_i['sigma'] ** 3 * sigma_epsilon_j['sigma'] ** 3 )
    return new_params


def fuzzy_list_match(line, ldata):
    """
    Searches for a line in a list of lines and returns the match if found.

    Examples
    --------

    >>> tmp = fuzzy_list_match("data tmp", ["other", "data", "else"])
    >>> print(tmp)
    (True, "data")

    >>> tmp = fuzzy_list_match("thing", ["other", "else"])
    >>> print(tmp)
    (False, None)

    """

    for match in ldata:
        if match in line:
            return True, match

    return False, None


def read_lines(filename, nlines=-1, start=0):
    """
    Reads the first nlines of a file with a `start` offset. Care is taken
    """

    if not os.path.isfile(filename):
        raise OSError("Could not find file '%s'" % filename)

    ret_data = []
    with open(filename, "r") as infile:

        # Advance to start
        for num in range(start):
            next(infile)

        # Read in the data
        if nlines == -1:
            for line in infile:
                ret_data.append(line.strip())
        else:
            for num in range(nlines):
                try:
                    ret_data.append(next(infile).strip())
                except StopIteration:
                    break

    return ret_data


def find_lowest_hole(data):
    """
    Finds the next lowest value in a list


    >>> find_lowest([0, 1, 3, 4])
    2
    """

    possible_values = set(range(len(data) + 1))
    new_key = min(possible_values - set(data))
    return new_key


def _build_hash_string(data, float_fmt):

    ret = []
    if isinstance(data, (str)):
        ret.append(data)
        ret.append(", ")
    elif isinstance(data, (int, float, np.int, np.int32, np.int64, np.float, np.float32, np.float64)):
        ret.append(float_fmt % data)
        ret.append(", ")
    elif isinstance(data, (tuple, list, np.ndarray)):
        ret.append(" (")
        for item in data:
            ret.append(_build_hash_string(item, float_fmt))
        ret.append("), ")
    elif isinstance(data, dict):
        ret.append("{")
        for k in sorted(data):
            ret.append(k)
            ret.append(": ")
            ret.append(_build_hash_string(data[k], float_fmt))
        ret.append("}, ")
    else:
        raise TypeError("hash: data type not understood: '%s'" % type(data))

    return "".join(ret)


def hash(data, rtol=8):
    """
    A special hashing function to deal with floating point numbers.
    """

    # Build up formatters
    float_fmt = "%." + str(rtol) + "f"

    serialized_data = _build_hash_string(data, float_fmt).encode()
    return hashlib.md5(serialized_data).hexdigest()
