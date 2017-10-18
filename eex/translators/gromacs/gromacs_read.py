"""
GROMACS EEX I/O
"""

import pandas as pd
import os

import numpy as np
import eex

import logging
logger = logging.getLogger(__name__)


def read_gromacs_gro_file(dl, gro_folder, ffdir=None):

    if ffdir is None:
        if "GROMACS_DIR" in os.environ:
            ffdir = os.environ["GROMACS_DIR"]
        else:
            raise KeyError("GROMACS read: Must provide `ffdir` if 'GROMACS_DIR' not in environmental variables.")

    ### Read in conf.gro file first

    conf_fname = os.path.join(gro_folder, "conf.gro")
    if not os.path.exists(conf_fname):
        raise OSError("GROMACS read: Could not find conf.gro file, expected at '%s'." % conf_fname)

    with open(conf_fname, 'r') as conf_file:
        conf_title = next(conf_file)
        conf_read_size = int(next(conf_file))

    reader = pd.read_table(
        conf_fname,
        header=None,
        iterator=True,
        names=["atom_name", "atomic_symbol", "atom_index", "X", "Y", "Z"],
        engine="c",
        comment=";",
        delim_whitespace=True,
        skiprows=2)

    data = reader.get_chunk(conf_read_size).dropna(axis=1, how="all")
    data["atomic_number"] = 0
    data.set_index("atom_index", inplace=True, drop=True)
    for val, idx in data.groupby("atomic_symbol"):
        data.loc[idx.index, "atomic_number"] = eex.metadata.atom_symbol_to_number[val]

    dl.add_atoms(data)

    # Get the box size
    size = reader.get_chunk(1).dropna(axis=1, how="all")
    half_box_length = size.values[0] / 2
    box_size = {k: (-v, v) for k, v in zip(["x", "y", "z"], half_box_length)}
    dl.set_box_size(box_size)

    ### Read in topol.top file next

    top_fname = os.path.join(gro_folder, "topol.top")
    if not os.path.exists(conf_fname):
        raise OSError("GROMACS read: Could not find topol.top file, expected at '%s'." % conf_fname)

    data = pd.read_table(
        top_fname,
        header=None,
        # iterator=True,
        names=range(9),
        engine="c",
        comment=";",
        delim_whitespace=True)

    print("")
    data = data.loc[~data.iloc[:, 0].str.contains("#")]
    indices = np.where(data.iloc[:, 0].str.contains("\["))[0]
    for indx in range(indices.shape[0] - 1):
        start = indices[indx]
        end = indices[indx + 1]

        label = data.iloc[start, 1]
        tmp = data.iloc[(start + 1):end].convert_objects(convert_numeric=True).dropna(axis=1, how="all")
        print(label)
        print(tmp)
        print('----')

    # print(indices)
    # print("")
    # print(data)

    print(os.environ)
