"""
The primary EEX init file
"""

# Import individual layers
from . import datalayer
from . import energy_eval as energy_eval
from . import filelayer
from . import testing
from . import translators
from . import utility
from . import form_converters
from . import metadata
from .units import ureg


# Initializes the eex_init metadata
from . import eex_init
del eex_init

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
