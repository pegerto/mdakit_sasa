"""
mdakit-sasa
This kit allows the calculation of a solvent-accessible-surface area of a trajectory
"""

# Add imports here

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
