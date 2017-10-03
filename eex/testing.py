"""
A utility file for testing helpers
"""

import pandas as pd
import numpy as np
import copy


def df_compare(left, right, columns=None, atol=1.e-8, rtol=1.e-5, equal_nan=True):
    """
    Compares two dataframe in an approximate manner.

    Checks:
    - Columns
    - Indices
    - Float columns (through tolerances)
    - Integer columns
    - Other columns
    """

    if columns is not None:
        if not isinstance(columns, (list, tuple)):
            columns = [columns]
        col_set = set(columns)
        if not (set(left.columns) >= col_set):
            raise KeyError("Left DataFrame did not contain all tested columns")

        if not (set(right.columns) >= col_set):
            raise KeyError("Right DataFrame did not contain all tested columns")

        left = left[columns]
        right = right[columns]
    else:
        if set(left.columns) != set(right.columns):
            raise KeyError("Right and Left DataFrames do not have the same columns")

        # Order the dataframe
        left = left[right.columns]

    # Check index and sort
    assert right.index.equals(left.index)
    left = left.loc[right.index]

    # Check floats
    fcols = [name for name, tp in zip(left.columns, left.dtypes) if tp.kind == "f"]
    fclose = np.allclose(left[fcols], right[fcols], atol=atol, rtol=rtol, equal_nan=equal_nan)
    if not fclose:
        raise AssertionError("DF_compare: Mismatch in float columns.")

    # Check ints
    icols = [name for name, tp in zip(left.columns, left.dtypes) if tp.kind == "i"]
    iclose = np.allclose(left[icols], right[icols])
    if not fclose:
        raise AssertionError("DF_compare: Mismatch in integer columns.")

    # Check everything else
    remaining_cols = list(set(left.columns) - set(fcols) - set(icols))
    rclose = left[remaining_cols].equals(right[remaining_cols])
    if not rclose:
        raise AssertionError("DF_compare: Mismatch in non-numeric columns.")

    return True


def dict_compare(left, right, atol=1.e-9, rtol=1.e-5):
    """
    A testing function that attempts to compare two different complex dictionaries.
    This function can currently handle the following data types:
    - int
    - str
    - float
    - set, list, tuple
    - np.ndarray
    - pd.DataFrame
    """

    if set(left) != set(right):
        raise KeyError("Right and Left dicts do not contain the same keys.")

    for key in list(left):
        lv = left[key]
        rv = right[key]

        match = True
        if isinstance(lv, (int, str)):
            match = lv == rv
        elif isinstance(lv, set):
            match = lv == set(rv)
        elif isinstance(lv, (float, np.ndarray)):
            match = np.allclose(lv, rv, atol=atol, rtol=rtol)
        elif isinstance(lv, (pd.DataFrame)):
            match = df_compare(lv, rv, atol=atol, rtol=rtol)
        else:
            raise TypeError("dict_compare: Misunderstood compare type '%s'." % str(type(lv)))

        if match is False:
            raise AssertionError("dict_compare: Mismatch for key %s, comparing %s to %s" % (key, str(lv), str(rv)))

    return True


def dl_compare(left, right, atom_checks=["xyz"]):
    """
    Attempts to compare two dataframes

    Checks:
    - Term and parameter lengths
    - Similar term and parameters
    - one/two/three/four body lengths
    - one/two/three/four body values

    Currently assumes both DL's can be loaded into memory.
    """

    ### Compare all terms within the DL

    left_uids = left.list_parameter_uids()
    right_uids = left.list_parameter_uids()

    # First make sure the number of terms is the same
    for k in list(left_uids):
        if not set(left_uids[k]) == set(right_uids[k]):
            raise KeyError("dl_compare: Mismatch in the number of parameters between left (%d) and right (%d)." %
                           (len(left_uids[k]), len(right_uids[k])))

    # Make sure the terms in left matches the terms in right.
    conversion_dict = {}

    # Loop over orders
    for k in list(left_uids):

        # Loop over left uids
        conversion_dict[k] = {}
        ruid_tmps = right_uids[k][:]
        for luid in left_uids[k]:
            pl = left.get_parameters(k, luid)

            # Loop over right uid's popping ones we used
            for ruid in ruid_tmps:
                pr = left.get_parameters(k, ruid)

                # Check match, pop right uid, and break this loop back to luid iterator
                if (pr[0] == pl[0]) and dict_compare(pr[1], pl[1]):
                    conversion_dict[k][ruid] = luid
                    ruid_tmps.remove(ruid)
                    break

        # After all luid's have been used, our ruid list should be empty
        if len(ruid_tmps):
            raise KeyError("dl_compare: Did not find a match for all parameter terms")

    ### Find matching atoms and compare atom properties
    left_atom_missing = set(atom_checks) - set(left.list_atom_properties())
    if len(left_atom_missing):
        raise KeyError("dl_compare: left dataframe was missing %s atom properies" % str(left_atom_missing))

    right_atom_missing = set(atom_checks) - set(right.list_atom_properties())
    if len(right_atom_missing):
        raise KeyError("dl_compare: right dataframe was missing %s atom properies" % str(right_atom_missing))

    left_atom = left.get_atoms(atom_checks, by_value=True)
    right_atom = right.get_atoms(atom_checks, by_value=True)

    if left_atom.shape != right_atom.shape:
        raise IndexError("dl_compare: The number of atoms in the left and right DL's does not match.")

    # Assume order is the same for now
    assert df_compare(left_atom.reset_index(), right_atom.reset_index())

    # Reorder based on coordinates
    # left_coords = left_atom[["X", "Y", "Z"]].values
    # right_coords = left_atom[["X", "Y", "Z"]].values

    # tmp_mat = left_coords[:, None, :] - right_coords
    # distance_matrix = np.sqrt(np.einsum('ijk,ijk->ij', tmp_mat, tmp_mat))

    # if np.sum(distance_matrix.min(axis=0) > 1.e-6):
    #     raise IndexError("dl_compare: Not all coordintes match")

    # reorder = distance_matrix.argmin(axis=0)
    # if reorder.shape[0] != np.unqiue(reorder).shape[0]:
    #     raise IndexError("dl_compare: ")

    ### Find matching terms

    # Build uid to uid dict
    for order in list(conversion_dict):
        for k, v in conversion_dict[order].items():
            if k != v:
                raise KeyError("dl_compare: Does not yet support non-identical parameter key dictionaries")

    for order in list(conversion_dict):
        if len(conversion_dict[order]) < 1:
            continue

        assert df_compare(left.get_terms(order).reset_index(), right.get_terms(order).reset_index())

    return True
