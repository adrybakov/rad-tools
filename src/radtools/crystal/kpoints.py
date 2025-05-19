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

__all__ = ["Kpoints"]


class Kpoints:

    def __init__(
        self, b1, b2, b3, coordinates=None, names=None, labels=None, path=None, n=100
    ) -> None:
        redirect_to_wulfric()

    def add_hs_point(self, name, coordinates, label, relative=True):
        redirect_to_wulfric()

    @property
    def path(self):
        redirect_to_wulfric()

    @path.setter
    def path(self, new_path):
        redirect_to_wulfric()

    @property
    def path_string(self):
        redirect_to_wulfric()

    @property
    def n(self):
        redirect_to_wulfric()

    @n.setter
    def n(self, new_n):
        redirect_to_wulfric()

    @property
    def labels(self):
        redirect_to_wulfric()

    def coordinates(self, relative=False):
        redirect_to_wulfric()

    def points(self, relative=False):
        redirect_to_wulfric()

    def flatten_points(self, relative=False):
        redirect_to_wulfric()
