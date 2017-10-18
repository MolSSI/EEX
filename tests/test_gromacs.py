"""
Tests for GROMACS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import eex_find_files

@pytest.fixture(scope="module")
def propane_dl():
    dl = eex.datalayer.DataLayer("test_gromacs_read")
    gro_folder = eex_find_files.get_example_filename("gromacs", "alkanes", "propane")
    eex.translators.gromacs.read_gromacs_gro_file(dl, gro_folder)
    return dl


def test_gromacs_read_data(propane_dl):
    dl = propane_dl

    box_size = dl.get_box_size()
    assert box_size["x"][0] == pytest.approx(-2.60288, 1.e-6)
    assert box_size["x"][1] == pytest.approx(2.60288, 1.e-6)