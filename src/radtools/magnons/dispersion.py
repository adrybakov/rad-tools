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

__all__ = ["MagnonDispersion"]


class MagnonDispersion:

    def __init__(
        self,
        model,
        Q=None,
        n=None,
        nodmi=False,
        noaniso=False,
        custom_mask=None,
    ):
        redirect_to_magnopy()

    def J(self, k):
        redirect_to_magnopy()

    def A(self, k):
        redirect_to_magnopy()

    def B(self, k):
        redirect_to_magnopy()

    def C(self):
        redirect_to_magnopy()

    def h(self, k):
        redirect_to_magnopy()

    def omega(self, k, zeros_to_none=False, return_G=False, return_imaginary=False):
        redirect_to_magnopy()

    def omegas(self, kpoints, zeros_to_none=False):
        redirect_to_magnopy()

    def __call__(self, *args, **kwargs):
        redirect_to_magnopy()
