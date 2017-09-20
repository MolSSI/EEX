"""
Provides helper functions that validate the functional data.
"""

term_requied_fields = ['variables', 'store_name', 'store_indices', 'styles']
functional_form_required_fields = ['form', 'parameters', 'units', 'description']

from eex import units


def validate_functional_form_dict(name, functional_form):
    """
    Checks an individual
    """

    # First check if the required fields are present
    missing_fields = set(functional_form_required_fields) - set(list(functional_form))
    if len(missing_fields):
        raise KeyError("Validate term style: term '%s' is missing the %s fields." % (name, list(missing_fields)))

    # Make sure each parameter has an associated unit
    missing_param_units = set(functional_form["parameters"]) - set(list(functional_form["units"]))
    if len(missing_param_units):
        raise KeyError("Validate term style: term '%s' is missing the %s parameter units." %
                       (name, list(missing_param_units)))

    # Validate the unit contexts
    for k, v in functional_form["units"].items():
        try:
            units.ureg.check(v)
        except:
            raise KeyError("Validate term style: term '%s' has unknown dimension '%s'" % (name, v))

    return True
