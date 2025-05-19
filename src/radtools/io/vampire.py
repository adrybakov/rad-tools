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


__all__ = ["dump_vampire", "dump_mat", "dump_ucf"]
from radtools._redirect import *


def dump_vampire(
    spinham,
    seedname="vampire",
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    nologo=False,
):
    redirect_to_magnopy()


def dump_mat(spinham, filename=None, materials=None, nologo=False):

    redirect_to_magnopy()


def dump_ucf(
    spinham,
    filename=None,
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    nologo=False,
):

    redirect_to_magnopy()
