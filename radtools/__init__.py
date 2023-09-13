r"""
RAD-tools
"""

__version__ = "0.8.3"
__doclink__ = "rad-tools.org"
__git_hash__ = "5b4dc9b04aecb7efebb9cd710c02e0ad7fb68e44"
__release_date__ = "14 September 2023"

from . import (
    constants,
    crystal,
    decorate,
    dos,
    exceptions,
    geometry,
    io,
    magnons,
    numerical,
    score,
    spinham,
)
from .constants import *
from .crystal import *
from .decorate import *
from .dos import *
from .geometry import *
from .io import *
from .magnons import *
from .numerical import *
from .score import *
from .spinham import *

__all__ = ["__version__", "__doclink__", "__git_hash__", "__release_date__"]
__all__.extend(crystal.__all__)
__all__.extend(decorate.__all__)
__all__.extend(dos.__all__)
__all__.extend(spinham.__all__)
__all__.extend(io.__all__)
__all__.extend(magnons.__all__)
__all__.extend(score.__all__)
__all__.extend(constants.__all__)
__all__.extend(geometry.__all__)
__all__.extend(numerical.__all__)
