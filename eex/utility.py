"""
A EEX utility folder
"""

import os
import hashlib
from . import units
import numpy as np
from subprocess import PIPE, Popen

def canonicalize_energy_names(energy_dict, canonical_keys):
    """Adjust the keys in energy_dict to the canonical names.

    Parameters
    ----------
    energy_dict : dict
    engine : str

    Returns
    -------
    normalized : dict

    """
    # TODO: Look into creating an `EnergyDict` class.
#    normalized = OrderedDict.fromkeys(canonical_energy_names,
#                                      0 * units.kilojoules_per_mole)
    ret = dict()
    for key, energy in energy_dict.items():
        canonical_key = canonical_keys.get(key)
        if canonical_key is None:
            continue
        elif isinstance(canonical_key, list):
            for k in canonical_key:
                ret[k] = energy
        else:
            ret[canonical_key] = energy

#    if 'Non-bonded' in canonical_keys:
#        normalized['nonbonded'] = energy_dict['Non-bonded']
#    else:
#        normalized['nonbonded'] = normalized['vdw total'] + normalized['coulomb total'] + normalized['h-bond']

#    # could be a problem here since desmond calls UB angle as stretch = bond
#    normalized['bonded'] = (normalized['bond'] + normalized['angle'] + normalized['dihedral'])
#
    return ret


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def run_subprocess(cmd, stdout_path, stderr_path, stdin=None):
    """
        General method to run MM codes. Taken from InterMol.
    """
    proc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=stdin)
   
    with open(stdout_path, 'a') as stdout, open(stderr_path, 'a') as stderr:
        stdout.write(out)
        stderr.write(err)

    if proc.returncode != 0:
        raise OSError("Command %s failed. Exit code %d. Error %s. Check file %s" % (cmd, proc.returncode, str(err), stderr_path))
    return out

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


def read_lines(filename, nlines=-1, start=0):
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

        # Read in the data
        if nlines == -1:
            for line in infile:
                ret_data.append(line.strip())
        else:
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


def _build_hash_string(data, float_fmt):

    ret = []
    if isinstance(data, (str)):
        ret.append(data)
        ret.append(", ")
    elif isinstance(data, (int, float, np.int, np.int32, np.int64, np.float, np.float32, np.float64)):
        ret.append(float_fmt % data)
        ret.append(", ")
    elif isinstance(data, (tuple, list, np.ndarray)):
        ret.append(" (")
        for item in data:
            ret.append(_build_hash_string(item, float_fmt))
        ret.append("), ")
    elif isinstance(data, dict):
        ret.append("{")
        for k in sorted(data):
            ret.append(k)
            ret.append(": ")
            ret.append(_build_hash_string(data[k], float_fmt))
        ret.append("}, ")
    else:
        raise TypeError("hash: data type not understood: '%s'" % type(data))

    return "".join(ret)


def hash(data, rtol=8):
    """
    A special hashing function to deal with floating point numbers.
    """

    # Build up formatters
    float_fmt = "%." + str(rtol) + "f"

    serialized_data = _build_hash_string(data, float_fmt).encode()
    return hashlib.md5(serialized_data).hexdigest()
