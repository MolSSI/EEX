"""
Tests to compare EEX energies to reference energies computed by amber
"""

import os
import eex
import numpy as np
import pandas as pd
import pytest
import eex_find_files
import glob
import eex_build_dl

# Build out file list
_test_directories = ["alkanes", "alcohols", "cyclic"]
_test_systems = []

for test_dir in _test_directories:
    print(test_dir)
    test_dir = eex_find_files.get_example_filename("amber", test_dir)
    systems = glob.glob(test_dir + "/*.prmtop")
    for s in systems:
        path, file_name = os.path.split(s)
        # Grab test dir and system name
        _test_systems.append((path, file_name))

# List current energy tests
_energy_types = {"two-body" : "bond", "three-body": "angle", "four-body": "dihedral"}

def test_references(amber_references):

    test_dir, system_name = amber_references
    molecule = str(system_name.split("/")[-1]).split(".")[0]

    data, dl = eex_build_dl.build_dl("amber", test_dir, molecule)
    dl_energies = dl.evaluate(utype='kcal * mol ** -1')

    reference_file = eex_find_files.get_example_filename("amber", test_dir, "energies.csv")
    reference_energies = pd.read_csv(reference_file, header=0)

    reference = reference_energies.loc[reference_energies['molecule'] == molecule]

    for k in _energy_types:
        k_reference = reference[_energy_types[k]].get_values()[0]

        # Test against reference values in CSV -- absolute toloerance of 0.001 is used since amber only
        # outputs energies to third decimal place
        success = dl_energies[k] == pytest.approx(k_reference, abs=0.001)
        if not success:
            raise ValueError("AMBER Test failed, energy type %s.\n"
                             "   Test description: %s" % (k, molecule))


# Loop over amber test directories
@pytest.fixture(scope="module", params=_test_systems)
def amber_references(request):
    test_dir, test_system = request.param
    return (test_dir, test_system)
