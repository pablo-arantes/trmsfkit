"""
tRMSF
A package to perform time-dependent RMSF analysis on molecular dynamics data.
"""

# Add imports here
from .trmsfkit import trmsfkit

from importlib.metadata import version

__version__ = version("trmsfkit")
