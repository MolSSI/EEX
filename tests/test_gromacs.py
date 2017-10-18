"""
Tests for GROMACS IO
"""
import eex
import numpy as np
import pytest
import pandas as pd
import eex_find_files


@pytest.fixture(scope="module")
def nbutane_dl():
    dl = eex.datalayer.DataLayer("test_gromacs_read")
    gro_folder = eex_find_files.get_example_filename("gromacs", "alkanes", "nbutane")
    eex.translators.gromacs.read_gromacs_gro_file(dl, gro_folder)
    return dl


def test_gromacs_read_conf(nbutane_dl):
    dl = nbutane_dl

    box_size = dl.get_box_size()
    assert box_size["x"][0] == pytest.approx(-2.5737, 1.e-6)
    assert box_size["x"][1] == pytest.approx(2.5737, 1.e-6)

    data = dl.get_atoms(None)
    assert data.shape[0] == 4
    assert dl.get_atom_count() == 4

    assert np.allclose(data["atomic_number"], [6, 6, 6, 6])
    assert np.allclose(data[["X", "Y", "Z"]].min(axis=0), [-0.147, -0.046, -0.153])
    assert np.allclose(data[["X", "Y", "Z"]].max(axis=0), [0.0, 0.16, 0.0])
