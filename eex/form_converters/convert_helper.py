
"""
A helper function to access EEX converteres easier
"""

from .bond_converters import convert_bond_coeffs
from .angle_converters import convert_angle_coeffs
from .dihedral_converters import convert_dihedral_coeffs


def sanitize_term_order_name(order):
    if isinstance(order, str):
        order = order.lower()

    if order in [2, "two", "bond", "bonds"]:
        return 2
    elif order in [3, "three", "angle", "angles"]:
        return 3
    elif order in [4, "four", "dihedral", "dihedrals"]:
        return 4
    else:
        raise KeyError("EEX: Term order name '%s' not recognized." % str(order))


def convert_form(order, coeffs, origin, final):
    order = sanitize_term_order_name(order)

    if order == 2:
        ret = convert_bond_coeffs(coeffs, origin, final)
    elif order == 3:
        ret = convert_angle_coeffs(coeffs, origin, final)
    elif order == 4:
        ret = convert_dihedral_coeffs(coeffs, origin, final)
    else:
        raise KeyError("EEX: Term order %s not recognized." % str(order))

    return ret
    # # Do we have a name?
    # if name:
    #     tmpdata = tmpdata[name]
    # else:
