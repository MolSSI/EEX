"""
Contains all metadata and wrappers used within EEX
"""

# Bring in atom metadata
from .atom_fields import atom_property_to_column, _valid_atom_properties 

# Bring in term metadata
from .two_body_terms import two_body_metadata
from .three_body_terms import three_body_metadata

# Bring in validation function
from .validator import * 
