"""
Converts various angle forms to other equivalents.
"""
import numpy as np
from .converter_registry import register_converter


@register_converter(order=3)
def _convert_nothing(coeffs):
    return coeffs
