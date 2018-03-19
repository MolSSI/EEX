import numpy as np
import eex
import os
from . import lammps_metadata as lmd


def compute_lattice_constants(bsize, tilt_factors):

    for key in ["x", "y", "z"]:
        if key.lower() not in bsize and key.upper() not in bsize:
            raise KeyError("Could not find key '%s'." % key)

    for key in ["xy", "xz", "yz"]:
        if key.lower() not in tilt_factors and key.upper() not in tilt_factors:
            raise KeyError("Could not find key '%s'." % key)

    lx = bsize['x']
    ly = bsize['y']
    lz = bsize['z']

    xy = tilt_factors['xy']
    xz = tilt_factors['xz']
    yz = tilt_factors['yz']

    a = lx
    b = np.sqrt(np.power(ly, 2) + np.power(xy, 2))
    c = np.sqrt(np.power(lz, 2) + np.power(xz, 2) + np.power(yz, 2))

    if np.isclose(b * c, 0.0):
        raise ZeroDivisionError("One of the box sizes is zero")

    cos_alpha = (xy *  xz + ly * yz) / (b * c)
    cos_beta = xz / c
    cos_gamma = xy / b

    alpha = np.arccos(cos_alpha)
    beta = np.arccos(cos_beta)
    gamma = np.arccos(cos_gamma)

    return {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

def get_energies(input_file=None, lmp_path=None, unit_style=None):
    """Evaluate energies of LAMMPS files. Based on InterMol

    Args:
        input_file = path to input file (expects data file in same folder)
        lmp_path = path to LAMMPS binaries
    """

    for exe in ['lammps', 'lmp_mpi', 'lmp_serial', 'lmp_openmpi',
                'lmp_mac_mpi']:
        if eex.utility.which(exe):
            LMP_PATH = exe
            break
    else:
        LMP_PATH = None

    if LMP_PATH is None and lmp_path is None:
        raise OSError('Unable to find LAMMPS executables.')
    else:
        if LMP_PATH is not None and lmp_path is None:
            lmp_path = LMP_PATH
        if LMP_PATH is not None and lmp_path is not None:
            if LMP_PATH not in lmp_path:
                raise OSError('Specified different LAMMPS executables.')

    if not os.path.isfile(input_file):
        raise OSError("Could not find file '%s'" % input_file)


    saved_path = os.getcwd()
    directory, input_file = os.path.split(os.path.abspath(input_file))
    os.chdir(directory)
    stdout_path = os.path.join(directory, 'lammps_stdout.txt')
    stderr_path = os.path.join(directory, 'lammps_stderr.txt')
    cmd = [lmp_path, '-in', input_file]
    out = eex.utility.run_subprocess(cmd, 'lammps', stdout_path, stderr_path)
    os.chdir(saved_path)

    log_unit_style = _extract_unit_style(out)
    log_prop_table = lmd.build_prop_table(log_unit_style)

    if unit_style is not None:
        if unit_style not in lmd.units_style:
            raise KeyError(" '%s' unit style not recognized" % unit_style)
    else:
        unit_style = _extract_unit_style(out)

    prop_table = lmd.build_prop_table(unit_style)

    ret = _group_energy_terms(out)

#    for key in ret:
#        cf = eex.units.conversion_factor(log_prop_table[key]['utype'], prop_table[key]['utype'])
#        ret[key] *= cf

    return eex.utility.canonicalize_energy_names(ret, {k:v['canonical'] for (k, v) in log_prop_table.items()})

def _extract_unit_style(stdout):
    """Parse LAMMPS stdout to extract unit style. """
    lines = stdout.split("\n")
    line_nbr = ['Unit style' in item for item in lines].index(True)
    unit_style = lines[line_nbr].split()[-1]
    return unit_style

def _group_energy_terms(stdout):
    """Parse LAMMPS stdout to extract and group the energy terms in a dict. """

    lines = stdout.split("\n")
    line_nbr = [item.startswith('Step') for item in lines].index(True)
    k = lines[line_nbr].split()[1:]
    v = [float(item) for item in lines[line_nbr + 1].split()[1:]]
    return dict(zip(k, v))
