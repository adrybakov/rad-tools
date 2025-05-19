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

__all__ = ["Lattice"]


class Lattice:

    def __init__(self, *args, standardize=True, **kwargs) -> None:
        redirect_to_wulfric()

    @property
    def cell(self):
        redirect_to_wulfric()

    def _set_cell(self, new_cell, standardize=True):
        redirect_to_wulfric()

    @cell.setter
    def cell(self, new_cell):
        redirect_to_wulfric()

    @property
    def a1(self):
        redirect_to_wulfric()

    @property
    def a2(self):
        redirect_to_wulfric()

    @property
    def a3(self):
        redirect_to_wulfric()

    @property
    def a(self):
        redirect_to_wulfric()

    @property
    def b(self):
        redirect_to_wulfric()

    @property
    def c(self):
        redirect_to_wulfric()

    @property
    def alpha(self):
        redirect_to_wulfric()

    @property
    def beta(self):
        redirect_to_wulfric()

    @property
    def gamma(self):
        redirect_to_wulfric()

    @property
    def unit_cell_volume(self):
        redirect_to_wulfric()

    @property
    def parameters(self):
        redirect_to_wulfric()

    @property
    def conv_cell(self):
        redirect_to_wulfric()

    @property
    def conv_a1(self):
        redirect_to_wulfric()

    @property
    def conv_a2(self):
        redirect_to_wulfric()

    @property
    def conv_a3(self):
        redirect_to_wulfric()

    @property
    def conv_a(self):
        redirect_to_wulfric()

    @property
    def conv_b(self):
        redirect_to_wulfric()

    @property
    def conv_c(self):
        redirect_to_wulfric()

    @property
    def conv_alpha(self):
        redirect_to_wulfric()

    @property
    def conv_beta(self):
        redirect_to_wulfric()

    @property
    def conv_gamma(self):
        redirect_to_wulfric()

    @property
    def conv_unit_cell_volume(self):
        redirect_to_wulfric()

    @property
    def conv_parameters(self):
        redirect_to_wulfric()

    @property
    def reciprocal_cell(self):
        redirect_to_wulfric()

    @property
    def b1(self):
        redirect_to_wulfric()

    @property
    def b2(self):
        redirect_to_wulfric()

    @property
    def b3(self):
        redirect_to_wulfric()

    @property
    def k_a(self):
        redirect_to_wulfric()

    @property
    def k_b(self):
        redirect_to_wulfric()

    @property
    def k_c(self):
        redirect_to_wulfric()

    @property
    def k_alpha(self):
        redirect_to_wulfric()

    @property
    def k_beta(self):
        redirect_to_wulfric()

    @property
    def k_gamma(self):
        redirect_to_wulfric()

    @property
    def reciprocal_cell_volume(self):
        redirect_to_wulfric()

    @property
    def reciprocal_parameters(self):
        redirect_to_wulfric()

    @property
    def eps(self):
        redirect_to_wulfric()

    @property
    def eps_rel(self):
        redirect_to_wulfric()

    @eps_rel.setter
    def eps_rel(self, new_value):
        redirect_to_wulfric()

    def type(self, eps_rel=None):
        redirect_to_wulfric()

    @property
    def variation(self):
        redirect_to_wulfric()

    @property
    def name(self):
        redirect_to_wulfric()

    @property
    def pearson_symbol(self):
        redirect_to_wulfric()

    @property
    def crystal_family(self):
        redirect_to_wulfric()

    @property
    def centring_type(self):
        redirect_to_wulfric()

    def lattice_points(self, relative=False, reciprocal=False, normalize=False):
        redirect_to_wulfric()

    def voronoi_cell(self, reciprocal=False, normalize=False):
        redirect_to_wulfric()

    @property
    def kpoints(self):
        redirect_to_wulfric()

    @kpoints.setter
    def kpoints(self, new_kpoints):
        redirect_to_wulfric()
