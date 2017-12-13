"""
This file contains all of the code required to format EEX metadata
"""

from . import metadata as md
from .units import convert_contexts

# Build the utype for atom data
for k, v in md.atom_metadata.items():
    if md.atom_metadata[k]["units"] is None:
        continue
    md.atom_metadata[k]["utype"] = {k: convert_contexts(v) for k, v in md.atom_metadata[k]["units"].items()}
    # md.atom_metadata[k]["utype"] = list(md.atom_metadata[k]["utype"].values())[0]

# Build the utype for the term data
for md_body in [md.two_body_metadata, md.three_body_metadata, md.four_body_metadata]:
    for k, v in md_body["forms"].items():
        v["utype"] = {uk: convert_contexts(uv) for uk, uv in v["units"].items()}

# Build the utype for the nb data
for k, v in md.nb_metadata["forms"].items():
    for uk, uv in v.items():
        uv["utype"] = {uuk: convert_contexts(uuv) for uuk, uuv in uv["units"].items()}
    #print(md_nb)
    """
    for form_type in md_nb:
        print(form_type)
        print(md_nb[form_type])
        md_nb[form_type]["utype"] = {uk: convert_contexts(uv) for uk, uv in md_nb[form_type]["units"].items()}
        """

