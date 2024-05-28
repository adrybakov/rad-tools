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
Exchange Module describes the spin Hamiltonian, defined on some :py:class:`.Crystal`.
"""


from .constants import *
from .hamiltonian import *
from .parameter import *
from .template import *

__all__ = []
__all__.extend(constants.__all__)
__all__.extend(hamiltonian.__all__)
__all__.extend(parameter.__all__)
__all__.extend(template.__all__)
