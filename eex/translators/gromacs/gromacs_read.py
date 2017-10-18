"""
GROMACS EEX I/O
"""

import pandas as pd
import os

import eex

import logging
logger = logging.getLogger(__name__)


def read_gromacs_gro_file(dl, gro_folder):

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
            names=["atom_name", "atomic_symbol", "atom_index", "x", "y", "z"],
            engine="c",
            comment=";",
            delim_whitespace=True,
            skiprows=2,
        )

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
