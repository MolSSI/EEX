"""
Test for LAMMPS/Amber comparison
"""

import eex
import numpy as np
import pytest
import eex_find_files

_alkane_molecules = ["ethane", "propane", "butane", "propane"]

def build_dl(program, molecule):
    if program.lower() == "amber":
        file_name = "trappe_%s_single_molecule.prmtop" % molecule
        fname = eex_find_files.get_example_filename("amber", "alkanes", file_name)
        dl = eex.datalayer.DataLayer("test_amber")
        data = eex.translators.amber.read_amber_file(dl, fname)
        return (data, dl)
    elif program.lower() == "lammps":
        file_name = "data.trappe_%s_single_molecule" % molecule
        fname = eex_find_files.get_example_filename("lammps", "alkanes", file_name)
        dl = eex.datalayer.DataLayer("test_lammps")
        data = eex.translators.lammps.read_lammps_file(dl, fname)
        return (data, dl)
    else:
        raise KeyError("Program %s not understood" % program)


@pytest.fixture(scope="module")
def ethane_bench(request):
    data, dl = build_dl("lammps", "ethane")
    yield dl
    dl.close()

@pytest.fixture(scope="module")
def propane_bench(request):
    data, dl = build_dl("lammps", "propane")
    yield dl
    dl.close()

@pytest.mark.parametrize("program", ["amber"])
def test_ethane(ethane_bench, program):

    test_dl = build_dl(program, "ethane")[1]
    assert eex.testing.dl_compare(ethane_bench, test_dl)

@pytest.mark.parametrize("program", ["amber"])
def test_propane(propane_bench, program):

    test_dl = build_dl(program, "propane")[1]
    assert eex.testing.dl_compare(propane_bench, test_dl)



