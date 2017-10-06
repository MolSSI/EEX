"""
Unit tests for EEX
"""
import pytest

import eex
from eex.units import convert_contexts, conversion_factor

Q_ = eex.ureg.Quantity


def test_default_contexts():

    assert convert_contexts("[energy] ** 2") == Q_("kilojoules ** 2 / mol ** 2")
    assert convert_contexts("[charge] ** 2") == Q_("elementary_charge ** 2")
    assert convert_contexts("[mass] ** 2") == Q_("g ** 2")
    assert convert_contexts("([mass] / [substance]) ** 2") == Q_("g ** 2 * mol ** -2")


def test_complex_contexts():
    assert convert_contexts("[energy] ** 2 * [length] ** -4") == Q_("kilojoules ** 2  mol ** -2 * angstrom ** -4")


def test_conversion_factor_magnitude():

    assert pytest.approx(10) == conversion_factor("10 * meter", "meter")
    assert pytest.approx(0.1) == conversion_factor("meter", "10 * meter")

    assert pytest.approx(2) == conversion_factor("10 * meter", "5 * meter")
    assert pytest.approx(0.5) == conversion_factor("5 * meter", "10 * meter")

    assert pytest.approx(2) == conversion_factor("10 * kcal/mol", "5 * kcal/mol")
    assert pytest.approx(0.5) == conversion_factor("5 * kcal/mol", "10 * kcal/mol")

    assert pytest.approx(2.8839179174400003e-18, 1.e-20) == conversion_factor("18 * elementary_charge", "coulomb")

    assert pytest.approx(1.e12) == conversion_factor("meter", "picometer")
