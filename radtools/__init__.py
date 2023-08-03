r"""
RAD-tools
"""

__version__ = "0.8.0"
__doclink__ = "rad-tools.org"
__git_hash__ = "9c9b087fa02be0cafdadaef6ec1c7926fe36e3d6"
__release_date__ = "2 August 2023"

from . import (
    crystal,
    decorate,
    dos,
    exchange,
    io,
    magnons,
    score,
    constants,
    geometry,
    numerical,
    exceptions,
)
from .crystal import *
from .decorate import *
from .dos import *
from .exchange import *
from .io import *
from .magnons import *
from .score import *
from .constants import *
from .geometry import *
from .numerical import *

__all__ = ["__version__", "__doclink__", "__git_hash__", "__release_date__"]
__all__.extend(crystal.__all__)
__all__.extend(decorate.__all__)
__all__.extend(dos.__all__)
__all__.extend(exchange.__all__)
__all__.extend(io.__all__)
__all__.extend(magnons.__all__)
__all__.extend(score.__all__)
__all__.extend(constants.__all__)
__all__.extend(geometry.__all__)
__all__.extend(numerical.__all__)
