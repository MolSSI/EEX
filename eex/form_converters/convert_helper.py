
"""
A helper function to access EEX converteres easier
"""
from .converter_registry import REGISTERED_CONVERTERS
from .. import metadata
from . import dihedral_converters
from . import angle_converters
from . import bond_converters


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
    # # Do we have a name?
    # if name:
    #     tmpdata = tmpdata[name]
    # else:


# TODO temporary solution, need something better than this
# _dihedral_conversion_matrix = dict()

# for k, v in four_body_metadata['forms'].items():
#     to_canonical = '_' + k + '_to_' + v['canonical_form']
#     from_canonical = '_' + v['canonical_form'] + '_to_' + k
#     try:
#         _dihedral_conversion_matrix[k] = (v['parameters'], locals()[to_canonical], locals()[from_canonical])
#     except:
#         _dihedral_conversion_matrix[k] = (v['parameters'], _convert_nothing, _convert_nothing)


# def convert_dihedral_coeffs(coeffs, origin, final):

#     difference = set([origin, final]) - set(_dihedral_conversion_matrix.keys())
#     if (difference):
#         raise KeyError("Conversion cannot be made since %s is not in conversion matrix %s" %
#                        (difference, _dihedral_conversion_matrix.keys()))

#     difference = set(coeffs.keys()) - set(_dihedral_conversion_matrix[origin][0])
#     if (difference):
#         raise KeyError("The key %s in the coefficient dictionary is not in the list of allowed keys %s" %
#                        (difference, _dihedral_conversion_matrix[origin][0]))

#     internal = _dihedral_conversion_matrix[origin][1](coeffs)
#     external = _dihedral_conversion_matrix[final][2](internal)
#     return external
