"""
Tests for LAMMPS IO
"""
import eex
import numpy as np
import os
import pytest
import pandas as pd
import eex_find_files

def test_gromacs_read_gro():
    dl = eex.datalayer.DataLayer("test_gromacs_read")
    fname = eex_find_files.get_example_filename("gromacs", "alkanes", "nbutane", "conf.gro")
    eex.translators.gromacs.read_gromacs_gro_file(dl, fname)
