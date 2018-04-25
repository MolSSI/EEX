"""
Converts various NB forms to other equivalents. In addtion, programs combining rules
"""

from .metadata import mixing_rules

## LJ Conversions
def _LJ_ab_to_ab(coeffs):
    """
    Convert AB representation to AB representation of the LJ potential
    """
    return {'A': coeffs['A'], 'B': coeffs['B']}


def _LJ_epsilonsigma_to_ab(coeffs):
    """
    Convert epsilon/sigma representation to AB representation of the LJ
    potential
    """
    A = 4.0 * coeffs['epsilon'] * coeffs['sigma']**12.0
    B = 4.0 * coeffs['epsilon'] * coeffs['sigma']**6.0
    return {"A": A, "B": B}


def _LJ_ab_to_epsilonsigma(coeffs):
    """
    Convert AB representation to epsilon/sigma representation of the LJ
    potential
    """
    if (coeffs['A'] == 0.0 and coeffs['B'] == 0.0):
        return {"sigma": 0.0, "epsilon": 0.0}

    try:
        sigma = (coeffs['A'] / coeffs['B'])**(1.0 / 6.0)
        epsilon = coeffs['B']**2.0 / (4.0 * coeffs['A'])

    except ZeroDivisionError:
        raise ZeroDivisionError("Lennard Jones functional form conversion not possible, division by zero found.")

    return {"sigma": sigma, "epsilon": epsilon}


def _LJ_rminepsilon_to_ab(coeffs):
    """
    Convert rmin/epsilon representation to AB representation of the LJ
    potential
    """
    A = coeffs['epsilon'] * coeffs['Rmin']**12.0
    B = 2 * coeffs['epsilon'] * coeffs['Rmin']**6.0
    return {"A": A, "B": B}


def _LJ_ab_to_rminepsilon(coeffs):
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


_LJ_conversion_matrix = {
    'AB': (['A', 'B'], _LJ_ab_to_ab, _LJ_ab_to_ab),
    'epsilon/sigma': (['epsilon', 'sigma'], _LJ_epsilonsigma_to_ab, _LJ_ab_to_epsilonsigma),
    'epsilon/Rmin': (['epsilon', 'Rmin'], _LJ_rminepsilon_to_ab, _LJ_ab_to_rminepsilon),
}


def convert_LJ_coeffs(coeffs, origin, final):

    difference = set([origin, final]) - set(_LJ_conversion_matrix.keys())
    if (difference):
        raise KeyError("Conversion cannot be made since %s is not in conversion matrix %s" %
                       (difference, _LJ_conversion_matrix.keys()))

    difference = set(coeffs.keys()) - set(_LJ_conversion_matrix[origin][0])
    if (difference):
        raise KeyError("The key %s in the coefficient dictionary is not in the list of allowed keys %s" %
                       (difference, _LJ_conversion_matrix[origin][0]))

    internal = _LJ_conversion_matrix[origin][1](coeffs)
    external = _LJ_conversion_matrix[final][2](internal)
    return external


## LJ combining rules

def _lorentz_berthelot(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = (sigma_epsilon_i['sigma'] + sigma_epsilon_j['sigma']) * 0.5
    new_params['epsilon'] = (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['epsilon'])**(1. / 2.)

    return new_params

def _arithmetic(sigma_epsilon_i, sigma_epsilon_j):
    return _lorentz_berthelot(sigma_epsilon_i, sigma_epsilon_j)


def _geometric(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = (sigma_epsilon_i['sigma'] * sigma_epsilon_j['sigma'])**(1. / 2.)
    new_params['epsilon'] = (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['epsilon'])**(1. / 2.)

    return new_params


def _sixthpower(sigma_epsilon_i, sigma_epsilon_j):
    new_params = {}

    new_params['sigma'] = ((sigma_epsilon_i['sigma']**6. * sigma_epsilon_j['sigma'] ** 6.) / 2.) ** (1./6.)

    new_params['epsilon'] = (2 * (sigma_epsilon_i['epsilon'] * sigma_epsilon_j['epsilon']) ** (1./2.) * sigma_epsilon_i['sigma'] ** 3.
                            * sigma_epsilon_j['sigma'] ** 3.) /(sigma_epsilon_i['sigma'] ** 6. * sigma_epsilon_j['sigma'] ** 6.)
    return new_params

def mix_LJ(coeff_i, coeff_j, mixing_rule, origin="AB", final="AB"):
    # Calculate interactions between two atom types based on specified mixing rules

    # First convert from input form to internal AB representation
    internal_coeff_i = convert_LJ_coeffs(coeff_i, origin=origin, final="AB")
    internal_coeff_j = convert_LJ_coeffs(coeff_j, origin=origin, final="AB")

    # Convert from internal AB representation to epsilon/sigma
    sigma_epsilon_i = convert_LJ_coeffs(internal_coeff_i, origin="AB", final="epsilon/sigma")
    sigma_epsilon_j = convert_LJ_coeffs(internal_coeff_j, origin="AB", final="epsilon/sigma")

    # Calculate new parameters based on mixing rules
    mixing_rule = mixing_rule.lower()
    new_params = LJ_mixing_functions[mixing_rule](sigma_epsilon_i, sigma_epsilon_j)

    # Convert from epsilon-sigma to AB, then to final specified form. Double conversion is necessary because of
    # form of conversion matrix.
    convert_params_temp = convert_LJ_coeffs(new_params, origin="epsilon/sigma", final="AB")
    convert_params = convert_LJ_coeffs(convert_params_temp, origin="AB", final=final)

    return convert_params

# Build mixing rules conversion

LJ_mixing_functions = {}

for mix in mixing_rules:
    internal_function = "_" + mix

    LJ_mixing_functions[mix] = eval(internal_function)
