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


__all__ = ["load_template", "read_template", "dump_spinham_txt"]

from radtools._redirect import *


def load_template(filename):
    redirect_to_magnopy()


# For backward compatibility
read_template = load_template


def dump_spinham_txt(
    spinham,
    filename=None,
    anisotropic=True,
    matrix=True,
    dmi=True,
    template=None,
    decimals=4,
    additional_stats=None,
):
    redirect_to_magnopy()
