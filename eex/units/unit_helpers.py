"""
A file for unit helpers
"""

import re

from .ureg import ureg
from ..metadata import default_contexts


__all__ = ["conversion_factor", "convert_contexts"]

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

    return ureg.convert(1, base_unit, conv_unit)


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

        if uc not in list(default_contexts):
            raise KeyError("Input context '%s' not understood.")

        # Adds parenthesis just in case
        tmp_conv = "(" + default_contexts[uc] + ")"
        context = context.replace(uc, tmp_conv)

    return ureg.Quantity(context)

