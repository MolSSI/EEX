"""
A helper function to access EEX metadata easier
"""

from .atom_fields import atom_metadata
from .two_body_terms import two_body_metadata
from .three_body_terms import three_body_metadata


def get_term_data(order, name, field=None):
    if order in [2, "two"]:
        tmpdata = two_body_metadata[name]
    elif order in [3, "three"]:
        tmpdata = three_body_metadata[name]
    else:
        raise KeyError("get_term_data: Order %s not recognized." % str(order))

    if field is None:
        return tmpdata
    else:
        return tmpdata[field]
