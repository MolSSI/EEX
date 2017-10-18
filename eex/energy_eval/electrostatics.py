"""
Computes the electrostatics of a expression
"""

import itertools

import numpy as np


# def _coulomb_energy(qij, R):

#     return np.sum(qij / R)


# def lattice_sum(coords, charges, nboxes, boxlength, func="coulomb"):

#     chargeij = np.prod(list(itertools.product(charges, charges)), axis=1)
#     dRij = np.repeat(coords, natoms, axis=0) - np.tile(coords, (natoms, 1))
#     energy = 0.0
#     truncated_energy = 0.0
#     for nvecs in itertools.product(boxlist, boxlist, boxlist):
#         nvecs = np.array(nvecs)
#         shiftvec = nvecs * boxlength
#         absvecs = np.abs(nvecs)
#         if nvecs.any():
#             # This is an image box
#             shiftedRij = np.linalg.norm(dRij + shiftvec, axis=1)
#             eterm = 0.5 * _coulomb_energy(qij=chargeij, R=shiftedRij)
#         else:
#             # Handle the 'home' box
#             eterm = 0.0
#             for i in range(natoms):
#                 for j in range(i):
#                     dR = np.linalg.norm(coords[i] - coords[j])
#                     eterm += _coulomb_energy(qij=charges[i] * charges[j], R=dR)
#         eterm *= coulomb
#         if np.linalg.norm(absvecs) <= nboxes:
#             truncated_energy += eterm
#         energy += eterm
#     return energy, truncated_energy