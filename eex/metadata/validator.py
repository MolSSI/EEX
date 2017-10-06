"""
Provides helper functions that validate the functional data.
"""

from eex import units

term_requied_fields = ['variables', 'store_name', 'store_indices', 'forms']
functional_form_required_fields = ['form', 'parameters', 'units', 'description']


def validate_term_dict(name, functional_form, parameters, utype=None):
    """
    Validates a dictionary term
    """

    if isinstance(parameters, (list, tuple)):
        if len(parameters) != len(functional_form["parameters"]):
            raise ValueError("Validate term dict: Nuber of parameters passed is %d, expected %d for terms %s" %
                             (len(parameters), len(functional_form["parameters"]), name))
        params = list(parameters)
    elif isinstance(parameters, dict):
        params = []
        for key in functional_form["parameters"]:
            try:
                params.append(parameters[key])
            except KeyError:
                raise KeyError("Validate term dict: Did not find expected key '%s' from term '%s'." % (key, name))
    else:
        raise TypeError("Validate term dict: Parameter type '%s' not understood" % str(type(parameters)))

    for value in params:
        if not isinstance(value, (int, float)):
            raise TypeError("Validate term dict: Parameters must be floats, found type %s." % type(value))

    # Deal with units
    if utype is not None:
        if isinstance(utype, (list, tuple)):
            if len(utype) != len(functional_form["utype"]):
                raise ValueError("Validate term dict: Number of units passed is %d, expected %d for terms %s" %
                                 (len(utype), len(functional_form["utype"]), name))
            form_units = list(utype)
        elif isinstance(utype, dict):
            form_units = []
            for key in functional_form["parameters"]:
                try:
                    form_units.append(utype[key])
                except KeyError:
                    raise KeyError("Validate term dict: Did not find expected key '%s' from term '%s'." % (key, name))
        else:
            raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

        # Convert units to internal
        for x, key in enumerate(functional_form["parameters"]):
            cf = units.conversion_factor(form_units[x], functional_form["utype"][key])
            params[x] *= cf


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
            units.convert_contexts(v)
        except:
            raise KeyError("Validate term style: term '%s' has unknown dimension '%s'" % (name, v))

    return True


def validate_units(unit_data, context=None):
    """
    Attempts to parse a Pint unit expression
    """

    unit_data = units.ureg.parse_expression(unit_data)

    # Make sure the units are valid
    if context is not None:
        units.ureg.check(context)(lambda x: x)(unit_data)

    return unit_data
