"""
Provides helper functions that validate the functional data.
"""

term_requied_fields = ['variables', 'store_name', 'store_indices', 'styles']
term_style_required_fields = ['form', 'parameters', 'units', 'description']

from eex import units


def validate_term_style_dict(name, term_style):
    """
    Checks an individual
    """

    # First check if the required fields are present
    missing_fields = set(term_style_required_fields) - set(list(term_style))
    if len(missing_fields):
        raise KeyError("Validate term style: term '%s' is missing the %s fields." % (name, list(missing_fields)))

    # Make sure each parameter has an associated unit
    missing_param_units = set(term_style["parameters"]) - set(list(term_style["units"]))
    if len(missing_param_units):
        raise KeyError("Validate term style: term '%s' is missing the %s parameter units." %
                       (name, list(missing_param_units)))

    # Validate the unit contexts
    for k, v in term_style["units"].items():
        try:
            units.ureg.check(v)
        except:
            raise KeyError("Validate term style: term '%s' has unknown dimension '%s'" % (name, v))

    return True
