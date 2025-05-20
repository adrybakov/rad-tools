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


class Crystal:

    def __init__(
        self,
        lattice=None,
        atoms=None,
        relative=True,
        standardize=True,
        **kwargs,
    ) -> None:
        redirect_to_wulfric()

    def __iter__(self):
        redirect_to_wulfric()

    def __contains__(self, atom):
        redirect_to_wulfric()

    def __len__(self):
        redirect_to_wulfric()

    def __getitem__(self, name):
        redirect_to_wulfric()

    def __getattr__(self, name):
        redirect_to_wulfric()

    @property
    def lattice(self):
        redirect_to_wulfric()

    def add_atom(self, new_atom=None, relative=True, **kwargs):
        redirect_to_wulfric()

    def remove_atom(self, atom, index=None):
        redirect_to_wulfric()

    def get_atom(self, name, index=None, return_all=False):
        redirect_to_wulfric()

    def get_atom_coordinates(self, atom, R=(0, 0, 0), index=None, relative=True):
        redirect_to_wulfric()

    def get_vector(
        self,
        atom1,
        atom2,
        R=(0, 0, 0),
        index1=None,
        index2=None,
        relative=False,
    ):
        redirect_to_wulfric()

    def get_distance(
        self,
        atom1,
        atom2,
        R=(0, 0, 0),
        index1=None,
        index2=None,
        relative=False,
    ):
        redirect_to_wulfric()

    def find_primitive_cell(self):
        redirect_to_wulfric()

    def mag_dipdip_energy(self, na, nb, nc, progress_bar=True):
        redirect_to_wulfric()

    def converge_mag_dipdip_energy(
        self,
        start=(10, 10, 10),
        step=(10, 10, 10),
        eps=10e-3,
        progress_bar=True,
        verbose=False,
    ):
        redirect_to_wulfric()
