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


__all__ = ["SpinHamiltonian", "ExchangeHamiltonian"]

from radtools._redirect import *


class SpinHamiltonian:

    def __init__(self, crystal=None, notation=None, **kwargs):
        redirect_to_magnopy()

    def __len__(self):
        redirect_to_magnopy()

    # Notation attributes
    @property
    def notation(self):
        redirect_to_magnopy()

    @notation.setter
    def notation(self, new_notation):
        redirect_to_magnopy()

    @property
    def notation_string(self):
        redirect_to_magnopy()

    def set_interpretation(
        self, double_counting=None, spin_normalized=None, factor=None
    ):
        redirect_to_magnopy()

    @property
    def double_counting(self):
        redirect_to_magnopy()

    def _ensure_double_counting(self):
        redirect_to_magnopy()

    def _ensure_no_double_counting(self):
        redirect_to_magnopy()

    @double_counting.setter
    def double_counting(self, new_value):
        redirect_to_magnopy()

    @property
    def spin_normalized(self):
        redirect_to_magnopy()

    @spin_normalized.setter
    def spin_normalized(self, new_value):
        redirect_to_magnopy()

    @property
    def factor(self):
        redirect_to_magnopy()

    @factor.setter
    def factor(self, new_factor):
        redirect_to_magnopy()

    def __iter__(self):
        redirect_to_magnopy()

    def __contains__(self, key):
        redirect_to_magnopy()

    def __getitem__(self, key):
        redirect_to_magnopy()

    def __getattr__(self, name):
        redirect_to_magnopy()

    @property
    def crystal(self):
        redirect_to_magnopy()

    @property
    def cell_list(self):
        redirect_to_magnopy()

    @property
    def magnetic_atoms(self):
        redirect_to_magnopy()

    @property
    def number_spins_in_unit_cell(self):
        redirect_to_magnopy()

    @property
    def space_dimensions(self):
        redirect_to_magnopy()

    def __setitem__(self, key, value):
        redirect_to_magnopy()

    def add_bond(self, atom1, atom2, R, J=None, **kwargs):
        redirect_to_magnopy()

    def __delitem__(self, key):
        redirect_to_magnopy()

    def remove_bond(self, atom1, atom2, R):
        redirect_to_magnopy()

    def remove_atom(self, atom):
        redirect_to_magnopy()

    def filter(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        redirect_to_magnopy()

    def filtered(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        redirect_to_magnopy()

    def form_model(self, template):
        redirect_to_magnopy()

    def formed_model(self, template):
        redirect_to_magnopy()

    def ferromagnetic_energy(self, theta=0, phi=0):
        redirect_to_magnopy()

    def input_for_magnons(self, nodmi=False, noaniso=False, custom_mask=None):
        redirect_to_magnopy()


class ExchangeHamiltonian(SpinHamiltonian):
    pass
