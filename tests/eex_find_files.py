"""
Obtains filepaths from examples for testing.
"""

import os

test_dir = os.path.dirname(os.path.abspath(__file__))
examples_dir = os.path.join(test_dir, "..", "examples")


def get_example_filename(*filename):

    program = filename[0].lower()
    # Make sure we have written these programs
    if program not in ["amber", "lammps"]:
        raise KeyError("Examples for program %s not found!" % program)

    # Make sure file exists
    fname = os.path.join(examples_dir, *filename)
    if not os.path.exists(fname):
        raise OSError("File %s not found!" % fname)

    return fname


def get_scratch_directory(filename):

    scr_dir = os.path.join(test_dir, "scratch")
    if not os.path.exists(scr_dir):
        os.mkdir(scr_dir)

    return os.path.join(scr_dir, filename)
