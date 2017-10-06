"""
Contains functions to help translate functional forms
"""

def register_functional_forms(self, order, name, form=None, utype=None):
    """
    Registers a single functional forms with the DL object.

    Parameters
    ----------
    order : int
        The order of the functional form (2, 3, 4, ...)
    name : str
        The name of the functional form you are adding
    form : dict
        The metadata for the incoming form. Follows the term order descriptions.
    utype : {dict, pint.Unit}
        The unit type of the functional form.

    Returns
    -------
    pass : bool
        If the operation was sucessfull or not.


    Examples
    --------

    # Register form by passing in explicit data
    form_metadata = {
        "form": "K*(r-R0) ** 2",
        "parameters": ["K", "R0"],
        "units": {
            "K": "kcal * mol ** 2",
            "R0": "picometers"
        },
        "description": "This is a harmonic bond"
    }

    dl.register_functional_form(2, "custom_harmonic", form=form_metadata)

    # Register form by using EEX build-in forms and setting units
    dl.register_functional_form(2, "harmonic", units={"K": "(kcal / mol) / angstrom ** 2", "R0": "picometers"})
    """

    user_order = order
    order = metadata.sanitize_term_order_name(order)

    if order not in self._functional_forms:
        raise KeyError("DataLayer:register_functional_forms: Did not understand order key '%s'." % str(user_order))

    # We are using an internal form
    if (form is None) and (utype is None):
        raise KeyError(
            "DataLayer:register_functional_forms: Must either pass in form (external form) or units (internal form)."
        )

    elif form is None:
        if utype is None:
            raise KeyError("DataLayer:register_functional_forms: Must pass in units if using a EEX built-in form.")

        try:
            form = metadata.get_term_metadata(order, "forms", name)
        except KeyError:
            raise KeyError(
                "DataLayer:register_functional_forms: Could not find built in form of order `%s` and name `%s" %
                (str(user_order), str(name)))

        form["units"] = utype

    # We are using an external form
    else:
        if name in self._functional_forms[order]:
            raise KeyError(
                "DataLayer:register_functional_forms: Key '%s' has already been registered." % str(name))
        # Pass validator later

    assert metadata.validate_functional_form_dict(name, form)

    # Make sure the data is valid and add
    self._functional_forms[order][name] = form
