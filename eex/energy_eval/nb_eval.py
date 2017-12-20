"""
Computes the electrostatics of a expression
"""

import itertools

import numpy as np


def _coulomb_energy(qij, R):

    return np.sum(qij / R)


def _fast_coulomb_sum(coords, charges, shift, efunc):

    energy = 0.0

    for i in range(coords.shape[0]):
        tmp = coords + (shift - coords[i])
        dR = np.sqrt(np.einsum('ij,ij->i', tmp, tmp))
        energy += efunc(charges[i] * charges, dR)

    return 0.5 * energy


def lattice_sum(coords, charges, boxlength, nboxes, func="coulomb", return_shells=False, spherical_truncation=True):
    """
    Computes the direct latice sum electostatic energy.

    Parameters
    ----------
    coords : array_like
        A (N, 3) dimensional array of cartesian coordinates
    charges : array_like
        A (N,) dimensional array of charges with the order matching the `coords` array.
    boxlength : array_like
        A (3,) dimensional array of box sizes (x, y, z)
    nboxes : int
        The number of boxes in each direction to sum over
    spherical_truncation : bool, optional
        Truncate in roughly a spherical pattern from the home box, otherwise use a 3D box of nboxes in length.
    func : str, optional
        The type of operator to sum over
    return_shells : bool, optional
        Return the energy of each shell, if false just returns the total energy.
    """

    # Setup boxes
    halfbox = nboxes // 2
    boxlist = list(range(-halfbox, halfbox + 1))

    # Setup empty energy dict
    energy = {"home": 0.0}
    energy.update({k: 0.0 for k in range(1, int(halfbox**1.5) + 2)})

    # Handle 'home' box, only use lower triangular of pairs
    for i in range(1, coords.shape[0]):
        tmp = coords[:i] - coords[i]
        dR = np.sqrt(np.einsum('ij,ij->i', tmp, tmp))
        energy["home"] += _coulomb_energy(charges[i] * charges[:i], dR)

    # Image boxes
    for nvecs in itertools.product(boxlist, boxlist, boxlist):
        shell = int(np.linalg.norm(nvecs))
        # Skip the home box
        if shell == 0: continue

        # We only want the sphere, not the box
        if spherical_truncation and (shell > nboxes):
            continue

        shiftvec = np.array(nvecs) * boxlength
        energy[shell] += _fast_coulomb_sum(coords, charges, shiftvec, _coulomb_energy)

    # Sum up the energy
    energy["total"] = sum(v for k, v in energy.items())

    if return_shells:
        return energy
    else:
        return energy["total"]