"""
Converts various angle forms to other equivalents.
"""
import numpy as np
from ..metadata import three_body_metadata


def _convert_nothing(coeffs):
    return coeffs


# TODO temporary solution, need something better than this
_angle_conversion_matrix = dict()

for k, v in three_body_metadata['forms'].items():
    to_canonical = '_' + k + '_to_' + v['canonical_form']
    from_canonical = '_' + v['canonical_form'] + '_to_' + k
    try:
        _angle_conversion_matrix[k] = (v['parameters'], locals()[to_canonical], locals()[from_canonical])
    except:
        _angle_conversion_matrix[k] = (v['parameters'], _convert_nothing, _convert_nothing)


def convert_angle_coeffs(coeffs, origin, final):

    difference = set([origin, final]) - set(_angle_conversion_matrix.keys())
    if (difference):
        raise KeyError("Conversion cannot be made since %s is not in conversion matrix %s" %
                       (difference, _angle_conversion_matrix.keys()))

    difference = set(coeffs.keys()) - set(_angle_conversion_matrix[origin][0])
    if (difference):
        raise KeyError("The key %s in the coefficient dictionary is not in the list of allowed keys %s" %
                       (difference, _angle_conversion_matrix[origin][0]))

    internal = _angle_conversion_matrix[origin][1](coeffs)
    external = _angle_conversion_matrix[final][2](internal)
    return external
