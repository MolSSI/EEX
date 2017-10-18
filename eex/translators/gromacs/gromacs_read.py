"""
GROMACS EEX I/O
"""

import pandas as pd
import math

import eex

import logging
logger = logging.getLogger(__name__)


def read_gromacs_gro_file(dl, filename):

    box_size = {}

    with open(filename, 'r') as gro_file:
        title     = next(gro_file)
        read_size = int(next(gro_file))
        
    reader = pd.read_table(
            filename,
            header=None,
            iterator=True,
            names=range(6),
            engine="c",
            comment="#",
            delim_whitespace=True,
            skiprows=2,
            chunksize=read_size
        )
    data = reader.get_chunk(read_size).dropna(axis=1, how="all")

    size = reader.get_chunk(1).dropna(axis=1, how="all")
    half_box_length = size.values[0] / 2
    #Can use .to_dict() instead?
    box_size = dict(zip(['x', 'y', 'z'], zip(-half_box_length, half_box_length)))
    dl.set_box_size(box_size)
