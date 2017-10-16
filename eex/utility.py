"""
A EEX utility folder
"""

import os


def line_fuzzy_list(line, ldata):
    """
    Searches for a line in a list of lines and returns the match if found.

    Examples
    --------

    >>> tmp = line_fuzzy_list("thing", ["other", "something words", "else"])
    >>> print(tmp)
    (True, "something words")

    >>> tmp = line_fuzzy_list("thing", ["other", "else"])
    >>> print(tmp)
    (False None)

    """

    for match in ldata:
        if match in line:
            return (True, match)

    return (False, None)


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
