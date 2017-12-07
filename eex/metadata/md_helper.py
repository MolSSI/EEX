"""
A helper function to access EEX metadata easier
"""

from .atom_metadata import atom_metadata
from .two_body_terms import two_body_metadata
from .three_body_terms import three_body_metadata
from .four_body_terms import four_body_metadata
from .nb_terms import nb_metadata


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


def get_atom_metadata(propery_name, field_name=None):

    tmp = atom_metadata[property_name]
    if field_name:
        return tmp[field_name]
    else:
        return tmp


def get_term_metadata(order, name=None, field=None):
    order = sanitize_term_order_name(order)

    if order == 2:
        tmpdata = two_body_metadata
    elif order == 3:
        tmpdata = three_body_metadata
    elif order == 4:
        tmpdata = four_body_metadata
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


def get_nb_metadata(name=None, form=None):

    # Check for name
    if name:
        # Not sure of the reason for the inconsistency with above functions - here must use "forms" key
        tmp_data = nb_metadata["forms"][name]
    else:
        return nb_metadata

    if form:
        return tmp_data[form]
    else:
        return tmp_data
