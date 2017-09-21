"""
A wrapper for the pint ureg data
"""

# Import some nice unicode features
from __future__ import unicode_literals

import pint
import os

# We only want the ureg exposed
__all__ = ["ureg"]

# Find the local path
file_dir = os.path.dirname(os.path.abspath(__file__))
pint_default_units_path = os.path.join(file_dir, "default_units.txt")

# Build the ureg and set context
ureg = pint.UnitRegistry(pint_default_units_path)
ureg.enable_contexts('chemistry')
