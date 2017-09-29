"""
Unit tests for EEX
"""
import pytest

import eex
from eex.units import convert_contexts

Q_ = eex.units.ureg.Quantity


def test_default_contexts():

    assert convert_contexts("[energy] ** 2") == Q_("kilojoules ** 2 / mol ** 2")
    assert convert_contexts("[charge] ** 2") == Q_("elementary_charge ** 2")
    assert convert_contexts("[mass] ** 2") == Q_("g ** 2")
    assert convert_contexts("([mass] / [substance]) ** 2") == Q_("g ** 2 * mol ** -2")


def test_complex_contexts():
    assert eex.units.convert_contexts("[energy] ** 2 * [length] ** -4") == Q_(
        "kilojoules ** 2  mol ** -2 * angstrom ** -4")

    assert eex.units.convert_contexts("[energy] ** 2 * [length] ** -4") == Q_(
        "kilojoules ** 2  mol ** -2 * angstrom ** -4")
