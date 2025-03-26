# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

r"""
RAD-tools
"""


__version__ = "0.9.2"
__doclink__ = "rad-tools.org"
__git_hash__ = "2707dea716bf655b0d62172b853f5002bca808df"
__release_date__ = "26 March 2025"

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
