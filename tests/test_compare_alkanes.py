"""
Test for LAMMPS/Amber comparison
"""

import eex
import numpy as np
import pytest
import eex_find_files

_alkane_molecules = ["ethane", "propane", "butane"]


def build_dl(program, molecule):
    if program.lower() == "amber":
        file_name = "trappe_%s_single_molecule.prmtop" % molecule
        fname = eex_find_files.get_example_filename("amber", "alkanes", file_name)
        dl = eex.datalayer.DataLayer("test_amber")
        data = eex.translators.amber.read_amber_file(dl, fname)
        return dl
    elif program.lower() == "lammps":
        file_name = "in.%s" % molecule
        fname = eex_find_files.get_example_filename("lammps", "alkanes", file_name)
        dl = eex.datalayer.DataLayer("test_lammps")
        eex.translators.lammps.read_lammps_input_file(dl, fname)
        return dl
    else:
        raise KeyError("Program %s not understood" % program)

def build_dl2(program, fname):
    if program.lower() == "amber":
        dl = eex.datalayer.DataLayer("test_amber")
        eex.translators.amber.read_amber_file(dl, fname)
        return dl
    elif program.lower() == "lammps":
        dl = eex.datalayer.DataLayer("test_lammps")
        fname = fname + ".in"
        eex.translators.lammps.read_lammps_input_file(dl, fname)
        return dl
    else:
        raise KeyError("Program %s not understood" % program)



def write_dl(program, dl, oname):
    if program.lower() == "amber":
        oname = oname+'.prmtop'
        eex.translators.amber.write_amber_file(dl, oname)
        return oname
    elif program.lower() == "lammps":
        infile = oname + ".in"
        eex.translators.lammps.write_lammps_file(dl, oname, infile)
        return oname
    else:
        raise KeyError("Program %s not understood" % program)

# Loop over alkane molecules
@pytest.fixture(scope="module", params=_alkane_molecules)
def lammps_bench(request):
    molecule = request.param
    dl = build_dl("lammps", molecule)
    energy = dl.evaluate()
    yield (molecule, dl, energy)
    dl.close()


# Loop over programs
@pytest.mark.parametrize("program", ["amber"])
def test_alkane(lammps_bench, program):
    molecule, bench_dl, bench_energy = lammps_bench

    test_dl = build_dl(program, molecule)
    #assert eex.testing.dl_compare(bench_dl, test_dl)
    assert eex.testing.dict_compare(bench_energy, test_dl.evaluate())

# write test that loads in amber --> lammps --> EEX compare datalayers

@pytest.mark.parametrize("program1", ["amber"]) #build dl
@pytest.mark.parametrize("program2", ["amber", "lammps"]) #write dl
@pytest.mark.parametrize("molecule", _alkane_molecules)
def test_translation(program1, program2, molecule):

    oname = eex_find_files.get_scratch_directory("test_output")

    original_dl = build_dl(program1, molecule)
    original_energy = original_dl.evaluate()

    original_box = original_dl.get_box_size()

    oname = write_dl(program2, original_dl, oname)

    new_dl = build_dl2(program2, oname)
    new_energy = new_dl.evaluate()

    new_box = new_dl.get_box_size()

    assert(eex.testing.dict_compare(original_box, new_box))

    assert(eex.testing.dict_compare(original_energy, new_energy))
