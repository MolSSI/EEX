"""
Validates the EEX metadata
"""

import eex
import pytest

# Test the two-body forms
two_bm = eex.metadata.two_body_metadata


@pytest.mark.parametrize("functional_form", list(two_bm["forms"]))
def test_two_body_style_metadata(functional_form):
    assert eex.metadata.validate_functional_form_dict(functional_form, two_bm["forms"][functional_form])


# Test the three-body forms
three_bm = eex.metadata.three_body_metadata


@pytest.mark.parametrize("functional_form", list(three_bm["forms"]))
def test_three_body_style_metadata(functional_form):
    assert eex.metadata.validate_functional_form_dict(functional_form, three_bm["forms"][functional_form])


# Test the four-body forms
four_bm = eex.metadata.four_body_metadata


@pytest.mark.parametrize("functional_form", list(four_bm["forms"]))
def test_four_body_style_metadata(functional_form):
    assert eex.metadata.validate_functional_form_dict(functional_form, four_bm["forms"][functional_form])
