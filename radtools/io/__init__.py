r"""
The module provides input-output routines. 
It does not mean to absorb all interfaces to the external data formats,
but designed to be the place for the constructor of the internal 
data structures from the input data of the external programs, 
as well as from the internal-specified formats.
"""

from .internal import *
from .tb2j import *
from .vampire import *

__all__ = []
__all__.extend(internal.__all__)
__all__.extend(tb2j.__all__)
__all__.extend(vampire.__all__)
