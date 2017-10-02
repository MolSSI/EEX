"""
This file contains all of the code required to format EEX metadata
"""

from . import metadata as md
from .units import convert_contexts

# Build the utype for atom data
for k, v in md.atom_metadata.items():
    if md.atom_metadata[k]["units"] is None:
        continue
    md.atom_metadata[k]["utype"] = convert_contexts(v["units"])

# Build the utype for the term data
for md_body in [md.two_body_metadata, md.three_body_metadata, md.four_body_metadata]:
    for k, v in md_body["forms"].items():
        tmp = {}
        for uk, uv in v["units"].items():
            tmp[uk] = convert_contexts(uv)
        v["utype"] = tmp
