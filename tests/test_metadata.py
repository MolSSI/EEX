"""
Validates the EEX metadata
"""

import eex
import pytest


# Test the two-body styles
two_bm = eex.metadata.two_body_metadata
@pytest.mark.parametrize("functional_form", list(two_bm["styles"]))
def test_two_body_style_metadata(functional_form):
    assert eex.metadata.validate_functional_form_dict(functional_form, two_bm["styles"][functional_form])

# Test the three-body styles
three_bm = eex.metadata.three_body_metadata
@pytest.mark.parametrize("functional_form", list(three_bm["styles"]))
def test_three_body_style_metadata(functional_form):
    assert eex.metadata.validate_functional_form_dict(functional_form, three_bm["styles"][functional_form])
