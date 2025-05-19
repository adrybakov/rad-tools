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

from radtools._redirect import *

__all__ = ["load_poscar", "dump_poscar"]


def load_poscar(file_object=None, return_crystal=True, return_comment=False):
    redirect_to_wulfric()


def dump_poscar(
    crystal_like, file_object="POSCAR", comment=None, decimals=8, mode: str = "Direct"
):
    redirect_to_wulfric()
