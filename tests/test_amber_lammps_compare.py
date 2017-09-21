"""
Test for LAMMPS/Amber comparison
"""

import eex
import numpy as np
import pytest
import eex_find_files


@pytest.fixture(scope="module", params=["HDF5", "Memory"])
def butane_amber_dl(request):
    fname = eex_find_files.get_example_filename("amber", "alkanes/butane_single_molecule.prmtop")
    dl = eex.datalayer.DataLayer("test_amber_read", backend=request.param)
    data = eex.translators.amber.read_amber_file(dl, fname)
    yield (data, dl)
    del dl

@pytest.fixture(scope="module", params=["HDF5", "Memory"])
def butane_lammps_dl(request):
    fname = eex_find_files.get_example_filename("lammps", "C4/data.butane")
    dl = eex.datalayer.DataLayer("test_amber_read", backend=request.param)
    data = eex.translators.lammps.read_lammps_file(dl, fname)
    yield (data, dl)
    del dl


# Test bond data

def test_lammmps_amber_bonds(butane_lammps_dl, butane_amber_dl):
    amber_data, amber_dl = butane_amber_dl
    lammps_data, lammps_dl = butane_lammps_dl

    amber_bonds = amber_dl.get_bonds()
    lammps_bonds = lammps_dl.get_bonds()

    assert amber_bonds.shape[0] == lammps_bonds.shape[0]
    assert set(np.unique(amber_bonds["bond_type"])) == set(np.unique(lammps_bonds["bond_type"]))

"""
# Test atom data
def test_lammps_amber_atoms(butane_lammps_dl, butane_amber_dl):
    amber_data, amber_dl = butane_amber_dl
    lammps_data, lammps_dl = butane_lammps_dl

    amber_atoms = amber_dl.get_atoms(properties=["mass"])
    lammps_atoms = lammps_dl.get_atoms(properties=["mass"])
"""