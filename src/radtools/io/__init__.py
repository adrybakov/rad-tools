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
The module provides input-output routines.
It does not mean to absorb all interfaces to the external data formats,
but designed to be the place for the constructor of the internal
data structures from the input data of the external programs,
as well as from the internal-specified formats.
"""

from .internal import *
from .tb2j import *
from .vampire import *
from .vasp import *

__all__ = []
__all__.extend(internal.__all__)
__all__.extend(tb2j.__all__)
__all__.extend(vampire.__all__)
__all__.extend(vasp.__all__)
