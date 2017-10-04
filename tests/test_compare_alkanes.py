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

# Loop over alkane molecules
@pytest.fixture(scope="module", params=_alkane_molecules[:2])
def lammps_bench(request):
    molecule = request.param
    data, dl = build_dl("lammps", molecule)
    yield (dl, molecule)
    dl.close()

# Loop over programs
@pytest.mark.parametrize("program", ["amber"])
def test_alkane(lammps_bench, program):
    bench_dl, molecule = lammps_bench

    test_dl = build_dl(program, molecule)[1]
    assert eex.testing.dl_compare(bench_dl, test_dl)


