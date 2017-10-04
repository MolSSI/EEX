"""
Tests the energy expression evaluation
"""

import eex
import pytest
import numpy as np

import eex.energy_eval as ee


def test_distance():

    assert 1.0 == ee.compute_distance([0, 0, 0], [0, 0, 1])[0]
    assert 0.0 == ee.compute_distance([0, 0, 0], [0, 0, 0])[0]


def test_angle():

    # Check all 90 degree domains
    p1 = [5, 0, 0]
    p2 = [0, 0, 0]
    p3 = [0, 2, 0]
    p4 = [0, 0, 4]
    assert 90.0 == ee.compute_angle(p1, p2, p3, degrees=True)[0]
    assert 90.0 == ee.compute_angle(p3, p2, p1, degrees=True)[0]

    assert 90.0 == ee.compute_angle(p1, p2, p4, degrees=True)[0]
    assert 90.0 == ee.compute_angle(p4, p2, p1, degrees=True)[0]

    assert 90.0 == ee.compute_angle(p3, p2, p4, degrees=True)[0]
    assert 90.0 == ee.compute_angle(p4, p2, p3, degrees=True)[0]

    # Zero angle
    p5 = [6, 0, 0]
    assert 0.0 == ee.compute_angle(p1, p2, p1, degrees=True)[0]
    assert 0.0 == ee.compute_angle(p1, p2, p5, degrees=True)[0]
    assert 0.0 == ee.compute_angle(p1, p2, p5)[0]

    # Linear
    p6 = [-5, 0, 0]
    p7 = [-4, 0, 0]
    assert 180.0 == ee.compute_angle(p1, p2, p6, degrees=True)[0]
    assert 180.0 == ee.compute_angle(p1, p2, p7, degrees=True)[0]
    assert np.pi == ee.compute_angle(p1, p2, p7)[0]
