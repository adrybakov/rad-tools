r"""
The module provides input-output routines. 
It does not mean to absorb all interfaces to the external data formats,
but designed to be the place for the constructor of the internal 
data structures from the input data of the external programs, 
as well as from the internal-specified formats.
"""

from .internal import read_template
from .tb2j import read_tb2j_model

__all__ = ["read_template", "read_tb2j_model"]
