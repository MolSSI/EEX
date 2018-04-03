"""
Contains tests for the utility functions within EEX
"""

import eex
import pytest
import numpy as np


def test_compute_lattice_constants():
    length = 10.0
    bsize = {'x': length, 'y': length, 'z': length}
    tilt_factors = {'xy': 0.0, 'xz': 0.0, 'yz': 0.0}
    lattice_const = eex.translators.lammps.lammps_utility.compute_lattice_constants(bsize, tilt_factors)

    assert np.isclose(lattice_const['alpha'], np.pi/2.0)
    assert np.isclose(lattice_const['beta'], np.pi/2.0)
    assert np.isclose(lattice_const['gamma'], np.pi/2.0)

def test_fuzzy_list_match():

    assert (False, None) == eex.utility.fuzzy_list_match("thing", [])
    assert (True, "data") == eex.utility.fuzzy_list_match("data tmp", ["other", "data", "else"])


def test_find_lowest_hole():

    assert 0 == eex.utility.find_lowest_hole([])
    assert 1 == eex.utility.find_lowest_hole([0])
    assert 2 == eex.utility.find_lowest_hole([0, 1, 3, 4])


def test_hash():

    # Quick hash
    assert '3e88c7ece2cbec43538d8898f8a7601e' == eex.utility.hash([5])
    assert 'a1b8364f4dd47beb4cc0dd364d78d05d' == eex.utility.hash(["something", 5])

    # Nested hash
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, [5, 6]])
    assert eex.utility.hash([5, (5, 6)]) == eex.utility.hash([5, np.array([5, 6])])
    assert eex.utility.hash([[5], [6]]) == eex.utility.hash(np.array([[5], [6]]))

    # Test numerical hashes
    assert eex.utility.hash([5.0]) == eex.utility.hash([5])
    assert eex.utility.hash([np.int(5)]) == eex.utility.hash([np.float64(5)])
    assert eex.utility.hash([5.0, 6]) == eex.utility.hash([5, 6.0])

    # Test dicts
    assert eex.utility.hash({"k": 5, "d": 6}) == eex.utility.hash({"d": 6.0, "k": 5.0 + 1.e-14})
    assert eex.utility.hash({"k": 5, "d": {"a": 6}}) == eex.utility.hash({"d": {"a": 6.0}, "k": 5.0 + 1.e-14})
    assert eex.utility.hash({"k": 5, "d": {"b": 6}}) != eex.utility.hash({"d": {"a": 6.0}, "k": 5.0 + 1.e-14})

    # Test tolerances
    assert eex.utility.hash(5) == eex.utility.hash(5.0)
    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-10])
    assert eex.utility.hash([5.0, 6.0]) == eex.utility.hash([5.0 + 1.e-10, 6.0 + 1.e-10])

    assert eex.utility.hash([5.0]) == eex.utility.hash([5.0 + 1.e-9])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-8])
    assert eex.utility.hash([5.0]) != eex.utility.hash([5.0 + 1.e-7])

    # Change tolerances
    assert eex.utility.hash([5.0], rtol=9) == eex.utility.hash([5.0 + 1.e-10], rtol=9)
    assert eex.utility.hash([5.0], rtol=9) != eex.utility.hash([5.0 + 1.e-9], rtol=9)

    # Fun gotchas
    assert eex.utility.hash([5.0], rtol=8) != eex.utility.hash([5.0], rtol=9)

    assert eex.utility.hash(("k", {"a": 5, "b": 6})) == eex.utility.hash(("k", {"b": 6, "a": 5}))
