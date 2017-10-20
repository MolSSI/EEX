"""
Contains tests for the utility functions within EEX
"""

import eex
import pytest
import numpy as np


def test_fuzzy_list_match():

    assert (False, None) == eex.utility.fuzzy_list_match("thing", [])
    assert (True, "data") == eex.utility.fuzzy_list_match("data tmp", ["other", "data", "else"])


def test_find_lowest_hole():

    assert 0 == eex.utility.find_lowest_hole([])
    assert 1 == eex.utility.find_lowest_hole([0])
    assert 2 == eex.utility.find_lowest_hole([0, 1, 3, 4])

def test_hash():

    # Quick hash
    assert '5523b81a4d1a4b1d3c4f5dc3dda4598d' == eex.utility.hash([5])
    assert '932107fbc6d59956139515c8ac2cdeed' == eex.utility.hash(["something", 5])

    # Nested hash
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, [5, 6]])
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, np.array([5, 6])])
    assert eex.utility.hash([[5], [6]]) == eex.utility.hash(np.array([[5], [6]]))

    # Test numerical hashes
    assert eex.utility.hash([5.0]) == eex.utility.hash([5])
    assert eex.utility.hash([np.int(5)]) == eex.utility.hash([np.float64(5)])
    assert eex.utility.hash([5.0, 6]) == eex.utility.hash([5, 6.0])

    # Test tolerances
    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-10])
    assert eex.utility.hash([5.0, 6.0]) == eex.utility.hash([5.0 + 1.e-10, 6.0 + 1.e-10])

    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-9])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-8])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-7])

    # Change tolerances
    assert eex.utility.hash([5.0], ftol=9) == eex.utility.hash([5.0 + 1.e-10], ftol=9)
    assert eex.utility.hash([5.0], ftol=9) != eex.utility.hash([5.0 + 1.e-9], ftol=9)

    # Fun gotchas
    assert eex.utility.hash([5.0], ftol=8) != eex.utility.hash([5.0], ftol=9)


