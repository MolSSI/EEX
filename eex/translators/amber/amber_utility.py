import numpy as np
import eex
import os
from . import amber_metadata as amd

def get_energies(prmtop=None, crd=None, input_file=None, amb_path=None):
    """Evaluate energies of AMBER files. Based on InterMol

    Args:
        input_file = path to input file (expects data file in same folder)
        lmp_path = path to LAMMPS binaries
    """

    if not prmtop:
        #prmtop = os.path.join(directory, 'parm.prmtop')
        raise OSError('Cannot find %s Amber parameter file' % prmtop)
    if not crd:
        #crd = os.path.join(directory, 'ener.edr')
        raise OSError('Cannot find %s Amber parameter file' % crdtop)
    directory, _ = os.path.split(os.path.abspath(prmtop))

    mdout = os.path.join(directory, 'amber.out')
    stdout_path = os.path.join(directory, 'amber_stdout.txt')
    stderr_path = os.path.join(directory, 'amber_stderr.txt')

    # Did they give a path, or the name of the file?
    is_last_bin = os.path.basename(os.path.normpath(amb_path))
    if is_last_bin == 'sander':
        amber_bin = amb_path
    else:
        amber_bin = os.path.join(amb_path, 'sander')

    if not eex.utility.which(amber_bin):
        raise OSError('Unable to find AMBER executable (sander).')

    # Run sander.
    cmd = [amber_bin, '-i', input_file, '-c', crd, '-p', prmtop, '-o', mdout, '-O']
    _ = eex.utility.run_subprocess(cmd, stdout_path, stderr_path)

    ret = _group_energy_terms(mdout)

    # TODO: Unit conversion

    return eex.utility.canonicalize_energy_names(ret, amd.to_canonical)



def _group_energy_terms(mdout):
    """Parse AMBER output file and group the energy terms in a dict. """

    with open(mdout) as f:
        all_lines = f.readlines()

    # Find where the energy information starts.
    for i, line in enumerate(all_lines):
        if line[0:8] == '   NSTEP':
            startline = i
            break
    else:
        raise AmberError('Unable to detect where energy info starts in AMBER '
                         'output file: {}'.format(mdout))

    # Strange ranges for amber file data.
    ranges = [[1, 24], [26, 49], [51, 77]]

    e_out = dict()
    potential = 0
    for line in all_lines[startline+3:]:
        if '=' in line:
            for i in range(3):
                r = ranges[i]
                term = line[r[0]:r[1]]
                if '=' in term:
                    energy_type, energy_value = term.split('=')
                    energy_value = float(energy_value)
                    potential += energy_value
                    energy_type = energy_type.rstrip()
                    e_out[energy_type] = energy_value
        else:
            break
    e_out['ENERGY'] = potential
#    eex.utility.canonicalize_energy_names(e_out)
    return e_out


