"""
Provides helper functions that validate the functional data.
"""

from eex import units

term_requied_fields = ['variables', 'store_name', 'store_indices', 'forms']
functional_form_required_fields = ['form', 'parameters', 'units', 'description']


def validate_term_dict(name, functional_form, parameters):
    """
    Validates a dictionary term
    """

    if isinstance(parameters, (list, tuple)):
        if len(parameters) != len(functional_form["parameters"]):
            raise ValueError("Validate term dict: Nunber of parameters passed is %d, expected %d for terms %s" %
                             (len(parameters), len(functional_form["parameters"]), name))
        params = list(parameters)
    elif isinstance(parameters, dict):
        params = []
        for key in functional_form["parameters"]:
            try:
                params.append(parameters[key])
            except KeyError:
                raise KeyError("Validate term dict: Did not find expected key '%s' from term '%s'." % (key, name))

    for value in params:
        if not isinstance(value, (int, float)):
            raise TypeError("Validate term dict: Parameters must be floats, found type %s." % type(value))

    # Cast any ints to floats
    return list(map(float, params))


def validate_functional_form_dict(name, functional_form):
    """
    Checks an individual functional form for corretness
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
