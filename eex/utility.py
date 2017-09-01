"""
A EEX utility folder
"""


def line_fuzzy_list(line, ldata):

    for match in ldata:
        if match in line:
            return (True, match)

    return (False, None)
