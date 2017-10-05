"""
Tests the energy expression evaluation
"""

import eex
import pytest
import numpy as np

np.random.seed(0)


def test_distance():

    assert pytest.approx(1.0) == eex.energy_eval.compute_distance([0, 0, 0], [0, 0, 1])[0]
    assert pytest.approx(0.0) == eex.energy_eval.compute_distance([0, 0, 0], [0, 0, 0])[0]

    tmp1 = np.random.rand(20, 3) * 4
    tmp2 = np.random.rand(20, 3) * 4
    np_dist = np.linalg.norm(tmp1 - tmp2, axis=1)
    ee_dist = eex.energy_eval.compute_distance(tmp1, tmp2)
    assert np.allclose(np_dist, ee_dist)


def test_angle():
    def _test_angle(p1, p2, p3, value, degrees=True):
        tmp = eex.energy_eval.compute_angle(p1, p2, p3, degrees=degrees)
        assert pytest.approx(value) == tmp

    # Check all 90 degree domains
    p1 = [5, 0, 0]
    p2 = [0, 0, 0]
    p3 = [0, 2, 0]
    p4 = [0, 0, 4]
    _test_angle(p1, p2, p3, 90)
    _test_angle(p3, p2, p1, 90)

    _test_angle(p1, p2, p4, 90)
    _test_angle(p4, p2, p1, 90)

    _test_angle(p3, p2, p4, 90)
    _test_angle(p4, p2, p3, 90)
    _test_angle(p4, p2, p3, np.pi / 2, degrees=False)

    # Zero angle
    p5 = [6, 0, 0]
    _test_angle(p1, p2, p1, 0)
    _test_angle(p1, p2, p5, 0)
    _test_angle(p1, p2, p5, 0, degrees=False)

    # Linear
    p6 = [-5, 0, 0]
    p7 = [-4, 0, 0]
    _test_angle(p1, p2, p6, 180)
    _test_angle(p1, p2, p6, 180)
    _test_angle(p1, p2, p7, np.pi, degrees=False)


def test_dihedral():
    def _test_dihedral(p1, p2, p3, p4, value, degrees=True):
        tmp = eex.energy_eval.compute_dihedral(p1, p2, p3, p4, degrees=degrees)
        assert pytest.approx(value) == tmp

    p1 = [0, 0, 0]
    p2 = [0, 2, 0]
    p3 = [2, 2, 0]
    p4 = [2, 0, 0]
    p5 = [2, 4, 0]

    # Cis check
    _test_dihedral(p1, p2, p3, p4, 0)

    # Trans check
    _test_dihedral(p1, p2, p3, p5, 180)

    # 90 phase checks
    p6 = [2, 2, -2]
    p7 = [2, 2, 2]
    _test_dihedral(p1, p2, p3, p6, -90)
    _test_dihedral(p1, p2, p3, p7, 90)

    p8 = [0, 4, 0]
    _test_dihedral(p8, p2, p3, p6, 90)
    _test_dihedral(p8, p2, p3, p7, -90)

    # Linear checks
    _test_dihedral(p1, p2, p3, [3, 2, 0], 0)
    _test_dihedral(p1, p2, p3, [3, 2 + 1.e14, 0], 180)
    _test_dihedral(p1, p2, p3, [3, 2 - 1.e14, 0], 0)


def test_evaluate():
    def _test_evaluate(value, expr, local, info=None):
        tmp = eex.energy_eval.evaluate_form(expr, local, global_dict=info)
        assert np.allclose(value, tmp)

    local_dict = {"a": np.arange(5), "b": np.arange(5) * 2 + 1}
    global_dict = {"a": np.arange(4), "b": np.arange(4) * 2 + 1}

    _test_evaluate(local_dict["a"] + 5, "a + 5", local_dict)
    _test_evaluate(local_dict["a"] + 5, "a + 5", local_dict, global_dict)
    _test_evaluate(global_dict["a"] + 5, "a + 5", {}, global_dict)

    # Logarithms and exponentials
    _test_evaluate(local_dict["b"], "log(exp(b))", local_dict)
    _test_evaluate(np.log(local_dict["b"]), "log(b)", local_dict)
    _test_evaluate(local_dict["b"], "log10(10**b)", local_dict)
    _test_evaluate(np.log10(local_dict["b"]), "log10(b)", local_dict)

    # Sums
    _test_evaluate(np.sum(local_dict["a"] ** 2), "sum(a ** 2)", local_dict)


