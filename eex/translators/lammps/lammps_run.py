import eex
import os
import .lammps_metadata

for _exe in ['lammps', 'lmp_mpi', 'lmp_serial', 'lmp_openmpi',
            'lmp_mac_mpi']:
    if eex.utility.which(_exe):
        _LMP_PATH = _exe
        break
else:
    _LMP_PATH = None


def run(input_file, lmp_path=None):
    """Evaluate energies of LAMMPS files. Based on InterMol

    Args:
        input_file = path to input file (expects data file in same folder)
        lmp_path = path to LAMMPS binaries
    """
    if _LMP_PATH is None and lmp_path is None:
        raise OSError('Unable to find LAMMPS executables.')
    else:
        if _LMP_PATH is not None and lmp_path is None:
            lmp_path = _LMP_PATH
        if _LMP_PATH is not None and lmp_path is not None:
            if _LMP_PATH not in lmp_path:
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
    print(_group_energy_terms(out))

def _group_energy_terms(stdout):
    """Parse LAMMPS stdout to extract and group the energy terms in a dict. """

    lines = stdout.split("\n")
    line_nbr = [item.startswith('Step') for item in lines].index(True)
    k = lines[line_nbr].split()[1:]
    v = lines[line_nbr + 1].split()[1:]
    return dict(zip(k, v))
