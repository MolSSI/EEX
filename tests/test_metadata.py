"""
Validates the EEX metadata
"""

import eex
import pytest
import numpy as np
import random
import numexpr as ne

np.random.seed(0)
random.seed(0)

# Test the two-order forms
term_dict = {}
term_dict["two"] = eex.metadata.two_body_metadata
term_dict["three"] = eex.metadata.three_body_metadata
term_dict["four"] = eex.metadata.four_body_metadata

form_list = [("two", form) for form in list(term_dict["two"]["forms"])]
form_list += [("three", form) for form in list(term_dict["three"]["forms"])]
form_list += [("four", form) for form in list(term_dict["four"]["forms"])]


@pytest.mark.parametrize("order,form", form_list)
def test_style_metadata(order, form):
    md = term_dict[order]["forms"][form]
    assert eex.metadata.validate_functional_form_dict(form, md)


# form_list = [("two", form) for form in list(term_dict["two"]["forms"])]
# form_list += [("three", form) for form in list(term_dict["three"]["forms"])]


@pytest.mark.parametrize("order,form", form_list)
def test_evaluate_metadata(order, form):
    md = term_dict[order]["forms"][form]

    global_dict = {k: (np.random.rand(4) + 1) for k in list(term_dict[order]["variables"])}
    local_dict = {k: random.random() for k in md["parameters"]}
    for ival in ["n", "phase"]:
        if ival in local_dict:
            local_dict[ival] = int(local_dict[ival])

    if form == "oxdna/fene":
        local_dict["R0"] = np.mean(global_dict["r"])

    # First compare that the input is valid
    expr = eex.energy_eval.evaluate_form(md["form"], local_dict, evaluate=False)

    expr_names = set(expr.input_names)
    dict_names = set(list(local_dict) + list(global_dict) + ["PI"])
    missing_names = expr_names - dict_names
    if len(missing_names):
        raise KeyError("For term %s %s missing keys %s in either expression or parameters" % (order, form,
                                                                                              str(missing_names)))

    # Looks like we have a few singularities, wait on that
    assert np.all(np.isfinite(eex.energy_eval.evaluate_form(md["form"], local_dict, global_dict)))