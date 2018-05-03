
"""
A helper function to access EEX converteres easier
"""
from .converter_registry import REGISTERED_CONVERTERS
from .. import metadata
from . import dihedral_converters
from . import angle_converters
from . import bond_converters
from .convert_helper import sanitize_term_order_name


def convert_form(order, coeffs, origin, final):

    order = sanitize_term_order_name(order)

    term_md = metadata.get_term_metadata(order, "forms")

    # TODO temporary solution
    to_canonical = '_' + origin + '_to_' + term_md[origin]["canonical_form"]
    from_canonical = '_' + term_md[final]["canonical_form"] + '_to_' + final
    if order == 2:
        ret = convert_bond_coeffs(coeffs, origin, final)
    elif order == 3:
        ret = convert_angle_coeffs(coeffs, origin, final)
    elif order == 4:
        canonical_coeffs = REGISTERED_CONVERTERS[order][to_canonical](coeffs)
        ret = REGISTERED_CONVERTERS[order][from_canonical](canonical_coeffs)
    else:
        raise KeyError("EEX: Term order %s not recognized." % str(order))

    return ret
