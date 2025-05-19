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


class Atom:

    def __init__(
        self,
        name="X",
        position=None,
        spin=None,
        magmom=None,
        charge=None,
        index=None,
    ) -> None:
        redirect_to_wulfric()

    def __str__(self):
        redirect_to_wulfric()

    def __format__(self, format_spec):
        redirect_to_wulfric()

    def __eq__(self, other):
        redirect_to_wulfric()

    def __hash__(self):
        redirect_to_wulfric()

    def __neq__(self, other):
        redirect_to_wulfric()

    @property
    def position(self):
        redirect_to_wulfric()

    @position.setter
    def position(self, new_position):
        redirect_to_wulfric()

    @property
    def name(self):
        redirect_to_wulfric()

    @name.setter
    def name(self, new_name):
        redirect_to_wulfric()

    @property
    def type(self):
        redirect_to_wulfric()

    @property
    def index(self):
        redirect_to_wulfric()

    @index.setter
    def index(self, new_index):
        redirect_to_wulfric()

    @property
    def spin(self):
        redirect_to_wulfric()

    @spin.setter
    def spin(self, new_spin):
        redirect_to_wulfric()

    @property
    def spin_direction(self):
        redirect_to_wulfric()

    @spin_direction.setter
    def spin_direction(self, new_spin_direction):
        redirect_to_wulfric()

    @property
    def spin_vector(self):
        redirect_to_wulfric()

    @spin_vector.setter
    def spin_vector(self, new_spin_vector):
        redirect_to_wulfric()

    @property
    def magmom(self):
        redirect_to_wulfric()

    @magmom.setter
    def magmom(self, new_magmom):
        redirect_to_wulfric()

    @property
    def charge(self):
        redirect_to_wulfric()

    @charge.setter
    def charge(self, new_charge):
        redirect_to_wulfric()

    @property
    def fullname(self):
        redirect_to_wulfric()
