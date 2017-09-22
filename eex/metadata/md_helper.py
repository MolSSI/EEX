"""
A helper function to access EEX metadata easier
"""

from .atom_fields import atom_metadata
from .two_body_terms import two_body_metadata
from .three_body_terms import three_body_metadata

def sanitize_term_order_name(order):
    if isinstance(order, str):
        order = order.lower()

    if order in [2, "two", "bond", "bonds"]:
        return 2
    elif order in [3, "two", "angle", "angles"]:
        return 3
    elif order in [4, "four", "dihedral", "dihedrals"]:
        return 4
    else:
        raise KeyError("EEX: Term order name '%s' not recognized." % str(order))


def get_term_metadata(order, name=None, field=None):
    order = sanitize_term_order_name(order)

    if order == 2:
        tmpdata = two_body_metadata
    elif order == 3:
        tmpdata = three_body_metadata
    else:
        raise KeyError("EEX: Term order %s not recognized." % str(order))

    # Do we have a name?
    if name:
        tmpdata = tmpdata[name]
    else:
        return tmpdata

    # Do we have a field
    if field:
        return tmpdata[field]
    else:
        return tmpdata