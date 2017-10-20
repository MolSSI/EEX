"""
A EEX utility folder
"""

import os
import hashlib

import numpy as np


def fuzzy_list_match(line, ldata):
    """
    Searches for a line in a list of lines and returns the match if found.

    Examples
    --------

    >>> tmp = fuzzy_list_match("data tmp", ["other", "data", "else"])
    >>> print(tmp)
    (True, "data")

    >>> tmp = fuzzy_list_match("thing", ["other", "else"])
    >>> print(tmp)
    (False, None)

    """

    for match in ldata:
        if match in line:
            return True, match

    return False, None


def read_lines(filename, nlines, start=0):
    """
    Reads the first nlines of a file with a `start` offset. Care is taken
    """

    if not os.path.isfile(filename):
        raise OSError("Could not find file '%s'" % filename)

    ret_data = []
    with open(filename, "r") as infile:

        # Advance to start
        for num in range(start):
            next(infile)

        # Read int eh data
        for num in range(nlines):
            try:
                ret_data.append(next(infile).strip())
            except StopIteration:
                break

    return ret_data


def find_lowest_hole(data):
    """
    Finds the next lowest value in a list


    >>> find_lowest([0, 1, 3, 4])
    2
    """

    possible_values = set(range(len(data) + 1))
    new_key = min(possible_values - set(data))
    return new_key


def hash(data, rtol=8):
    """
    A special hashing function to deal with floating point numbers.
    """

    # Build up formatters
    float_fmt = "%." + str(rtol) + "f"

    m = hashlib.md5()
    for d in data:
        if isinstance(d, (str)):
            m.update(d.encode())
        elif isinstance(d, (int, float, np.int, np.int32, np.int64, np.float, np.float32, np.float64)):
            m.update((float_fmt % d).encode())
        elif isinstance(d, (tuple, list, np.ndarray)):
            m.update(hash(d).encode())
        else:
            raise TypeError("hash: Data type '%s' not understood." % type(d))

    return m.hexdigest()
