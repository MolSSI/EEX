"""
Tests to compare EEX energies to reference energies computed by amber
"""

import eex
import numpy as np
import pandas as pd
import pytest
import eex_find_files
import glob
import eex_build_dl

_test_directories = ["alkanes", "alcohols", "cyclic"]

_energy_types = {"two-body" : "bond", "three-body": "angle", "four-body": "dihedral"}

# Loop over amber test directories
@pytest.fixture(scope="module", params=_test_directories)
def get_references(request):
    test_directory = request.param
    reference_file = eex_find_files.get_example_filename("amber", test_directory, "energies.csv")
    reference_energies = pd.read_csv(reference_file, header=0)

    reference_directory = eex_find_files.get_example_filename("amber", test_directory)
    files = glob.glob(reference_directory + "/*.prmtop")
    yield (reference_energies, files, test_directory)

def test_references(get_references):
    ref_energy, file_list, test_directory = get_references
    for f in file_list:
        molecule = str(f.split("/")[-1]).split(".")[0]
        data, dl = eex_build_dl.build_dl("amber",test_directory, molecule)
        dl_energies = dl.evaluate(energy_unit = 'kcal * mol ** -1')

        reference = ref_energy.loc[ref_energy['molecule'] == molecule]

        for k in _energy_types:
            k_reference = reference[_energy_types[k]].get_values()[0]

            # Test against reference values in CSV - absolute toloerance of 0.001 is used since amber only
            # outputs energies to third decimal place
            assert dl_energies[k] == pytest.approx(k_reference, abs=0.001)
