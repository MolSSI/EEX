"""
A file for unit helpers
"""

import re
import pint

from .ureg import ureg
from ..metadata import default_contexts

__all__ = ["conversion_dict", "conversion_factor", "convert_contexts"]

def conversion_dict(from_units, to_units):
    """
    Builds a standard conversion dictionary between two dictionaries of units.
    """

    diff = set(from_units) ^ set(to_units)
    if len(diff):
        raise KeyError("conversion_dict: Not all units found, missing '%s'" % str(diff))

    ret = {}
    for key in from_units:
        ret[key] = conversion_factor(from_units[key], to_units[key])

    return ret


def conversion_factor(base_unit, conv_unit):
    """
    Provides the conversion factor from one unit to another

    Examples
    --------

    >>> conversion_factor("meter", "picometer")
    1e-12

    >>> conversion_factor(1, "feet", "meter")
    0.30479999999999996

    """

    # Add a little magic incase the incoming values have scalars

    factor = 1.0

    # First make sure we either have Units or Quantities
    if isinstance(base_unit, str):
        base_unit = ureg.parse_expression(base_unit)

    if isinstance(conv_unit, str):
        conv_unit = ureg.parse_expression(conv_unit)

    # Need to play with prevector if we have Quantities
    if isinstance(base_unit, pint.quantity._Quantity):
        factor *= base_unit.magnitude
        base_unit = base_unit.units

    if isinstance(conv_unit, pint.quantity._Quantity):
        factor /= conv_unit.magnitude
        conv_unit = conv_unit.units

    return ureg.convert(factor, base_unit, conv_unit)


def convert_contexts(context, converter=default_contexts):
    """
    Converts a context string to the

    Parameters
    ----------
    context : str
        The required context as a str
    converter : {dict}, optional
        A string to string dict that to convert a context to a specific value. Defaults to the
        EEX internal parameters.

    Returns
    -------
    q : pint.Quantity
        The converted unit as a Pint base class

    Examples
    --------
    >>> convert_context("[energy]")
    <Quantity(1.0, 'kilocalorie / mole')>

    >>> convert_context("[energy]", converter={"[energy]": "calorie / mol"})
    <Quantity(1.0, 'calorie / mole')>

    """

    if not isinstance(context, str):
        raise TypeError("Can only convert strings to units.")

    ucontexts = set(re.findall(r"\[([A-Za-z0-9_]+)\]", context))

    for uc in ucontexts:

        # Add brackets back after re stripped them
        uc = "[" + uc + "]"

        if uc not in list(converter):
            raise KeyError("Input context '%s' not understood." % uc)

        # Adds parenthesis just in case
        tmp_conv = "(" + converter[uc] + ")"
        context = context.replace(uc, tmp_conv)

    return ureg.Quantity(context)
