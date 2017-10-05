"""
Evaluates the EEX energy expresions.
"""

import numpy as np
import numexpr as ne
from . import metadata


def _norm(points):
    """
    Return the Frobenius norm across axis=-1, NumPy's internal norm is crazy slow (~4x)
    """

    tmp = np.atleast_2d(points)
    return np.sqrt(np.einsum("ij,ij->i", tmp, tmp))


def compute_distance(points1, points2):
    """
    Computes the pairwise distance between all points in points1 and points2.

    Parameters
    ----------
    points1 : np.ndarray
        The first list of points, can be 1D or 2D
    points2 : np.ndarray
        The second list of points, can be 1D or 2D

    Returns
    -------
    distances : np.ndarray
        The array of distances between points1 and points2

    Notes
    -----
    Units are not considered inside these expressions, please preconvert to the same units before using.
    """

    points1 = np.atleast_2d(points1)
    points2 = np.atleast_2d(points2)

    return _norm(points1 - points2)


def compute_angle(points1, points2, points3, degrees=False):
    """
    Computes the angle (p1, p2 [vertex], p3) between the provided points on a per-row basis.

    Parameters
    ----------
    points1 : np.ndarray
        The first list of points, can be 1D or 2D
    points2 : np.ndarray
        The second list of points, can be 1D or 2D
    points3 : np.ndarray
        The third list of points, can be 1D or 2D
    degrees : bool, options
        Returns the angle in degress rather than radians if True

    Returns
    -------
    angles : np.ndarray
        The angle between the three points in radians

    Notes
    -----
    Units are not considered inside these expressions, please preconvert to the same units before using.
    """

    points1 = np.atleast_2d(points1)
    points2 = np.atleast_2d(points2)
    points3 = np.atleast_2d(points3)

    v12 = points1 - points2
    v23 = points2 - points3

    denom = _norm(v12) * _norm(v23)
    cosine_angle = np.einsum("ij,ij->i", v12, v23) / denom

    angle = np.pi - np.arccos(cosine_angle)

    if degrees:
        return np.degrees(angle)
    else:
        return angle


def compute_dihedral(points1, points2, points3, points4, degrees=False):
    """
    Computes the dihedral angle (p1, p2, p3, p4) between the provided points on a per-row basis.

    Parameters
    ----------
    points1 : np.ndarray
        The first list of points, can be 1D or 2D
    points2 : np.ndarray
        The second list of points, can be 1D or 2D
    points3 : np.ndarray
        The third list of points, can be 1D or 2D
    points4 : np.ndarray
        The third list of points, can be 1D or 2D
    degrees : bool, options
        Returns the dihedral angle in degress rather than radians if True

    Returns
    -------
    dihedrals : np.ndarray
        The dihedral angle between the four points in radians

    Notes
    -----
    Units are not considered inside these expressions, please preconvert to the same units before using.
    """

    points1 = np.atleast_2d(points1)
    points2 = np.atleast_2d(points2)
    points3 = np.atleast_2d(points3)
    points4 = np.atleast_2d(points4)

    # Build the three vectors
    v12 = points1 - points2
    v23 = points2 - points3
    v34 = points3 - points4

    # Build vectors normal to the two planes
    n123 = np.cross(v12, v23)
    n234 = np.cross(v23, v34)
    n1234 = np.cross(n123, n234)

    left = np.einsum("ij,ij->i", n1234, v23) / _norm(v23)
    right = np.einsum("ij,ij->i", n123, n234)

    angle = np.arctan2(left, right)

    if degrees:
        return np.degrees(angle)
    else:
        return angle


def evaluate_form(form, parameters, global_dict=None, out=None, evaluate=True):
    """
    Evaluates a functional form from a string.

    """

    known_vals = {"PI": np.pi}

    if global_dict is None:
        global_dict = known_vals
    else:
        global_dict = global_dict.copy()
        global_dict.update(known_vals)

    if evaluate:
        return ne.evaluate(form, local_dict=parameters, global_dict=global_dict, out=out)
    else:
        return ne.NumExpr(form)


def evaluate_energy_expression(dl):

    energy = {2: 0.0, 3: 0.0, 4: 0.0}

    # Two-body
    xyz = dl.get_atoms("xyz")
    bonds = dl.get_bonds()
    print(xyz)
    print(bonds)

    # Loop over unique parameters/terms

    for idx, df in bonds.groupby("term_index"):
        two_body_dict = {}
        two_body_dict["r"] = compute_distance(xyz.loc[df["atom1"]].values, xyz.loc[df["atom2"]].values)

        form_type, parameters = dl.get_parameters(2, idx)

        form = metadata.get_term_metadata(2, "forms", form_type)["form"]

        energy[2] += evaluate_form(form, parameters, two_body_dict)

    print(energy)
