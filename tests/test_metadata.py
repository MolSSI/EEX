"""
Validates the EEX metadata
"""

import eex
import pytest


# Test the two-body styles
two_bm = eex.metadata.two_body_metadata
@pytest.mark.parametrize("term_style", list(two_bm["styles"]))
def test_two_body_style_metadata(term_style):
    assert eex.metadata.validate_term_style_dict(term_style, two_bm["styles"][term_style])

# Test the three-body styles
three_bm = eex.metadata.three_body_metadata
@pytest.mark.parametrize("term_style", list(three_bm["styles"]))
def test_three_body_style_metadata(term_style):
    assert eex.metadata.validate_term_style_dict(term_style, three_bm["styles"][term_style])
