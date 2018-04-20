"""
Test for LAMMPS/Amber comparison
"""

import eex
import numpy as np
import pytest
from . import eex_find_files

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
    """
    This test does two translations, and two energy comparisons

        program1 -> program2 -> program1

        At each stage, energies and atom metadata are compared.

    """
    oname = eex_find_files.get_scratch_directory("test_output")

    # Load data into 'program1'
    original_dl = build_dl(program1, molecule)
    original_energy = original_dl.evaluate()

    original_box = original_dl.get_box_size()

    original_atoms = original_dl.get_atoms(properties= ["mass", "atom_type", "charge", "xyz"], by_value=True)

    # Write new dl as 'program2'
    oname = write_dl(program2, original_dl, oname)

    # Read written dl
    new_dl = build_dl2(program2, oname)
    new_energy = new_dl.evaluate()

    new_box = new_dl.get_box_size()

    new_atoms = new_dl.get_atoms(properties=["mass", "atom_type", "charge", "xyz"], by_value=True)

    # Test that box parameters are the same
    assert(eex.testing.dict_compare(original_box, new_box))

    # Test that the system energy is the same
    assert(eex.testing.dict_compare(original_energy, new_energy))

     # Test that atom metadata is the same
    assert (eex.testing.dict_compare(original_atoms.to_dict(), new_atoms.to_dict()))

    # Translate back to original program
    oname = write_dl(program1, new_dl, oname)

    # Read written dl
    new_dl2 = build_dl2(program1, oname)
    new_energy2 = new_dl2.evaluate()
    new_atoms2 = new_dl2.get_atoms(properties=["mass", "atom_type", "charge", "xyz"], by_value=True)

    assert (eex.testing.dict_compare(original_energy, new_energy2))
    assert (eex.testing.dict_compare(original_atoms.to_dict(), new_atoms2.to_dict()))

    # Test NB parameters

    original_dl_nb_parameters = original_dl.list_nb_parameters(
        nb_name="LJ", nb_model="AB", itype="pair")

    new_dl2_nb_parameters = new_dl2.list_nb_parameters(
        nb_name="LJ", nb_model="AB", itype="pair")

    assert (eex.testing.dict_compare(original_dl_nb_parameters, new_dl2_nb_parameters))

