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

__all__ = ["MatplotlibBackend", "PlotlyBackend"]


class AbstractBackend:
    def __init__(self) -> None:
        redirect_to_wulfric()

    def plot(self, *args, kind, **kwargs):
        redirect_to_wulfric()

    def remove(self, *args, **kwargs):
        raise NotImplementedError

    def show(self, *args, **kwargs):
        redirect_to_wulfric()

    def save(self, *args, **kwargs):
        redirect_to_wulfric()

    def clear(self, *args, **kwargs):
        redirect_to_wulfric()

    def legend(self, *args, **kwargs):
        redirect_to_wulfric()

    # Backend-independent functions
    def plot_brillouin(self, *args, color="#FF4D67", **kwargs):
        redirect_to_wulfric()

    def plot_brillouin_kpath(
        self, *args, zone_color="#FF4D67", path_color="black", **kwargs
    ):
        redirect_to_wulfric()

    def plot_primitive(self, *args, **kwargs):
        redirect_to_wulfric()

    def plot_conventional(self, *args, **kwargs):
        redirect_to_wulfric()

    def plot_unit_cell(self, *args, **kwargs):
        redirect_to_wulfric()

    def plot_wigner_seitz(self, *args, **kwargs):
        redirect_to_wulfric()

    def plot_kpath(self, *args, **kwargs):
        redirect_to_wulfric()


class MatplotlibBackend(AbstractBackend):

    def __init__(self, fig=None, ax=None, background=True, focal_length=0.2):
        redirect_to_wulfric()

    def remove(self, kind="primitive"):
        redirect_to_wulfric()

    def plot(self, lattice, kind="primitive", **kwargs):
        redirect_to_wulfric()

    def show(self, elev=30, azim=-60):
        redirect_to_wulfric()

    def save(self, output_name="lattice_graph.png", elev=30, azim=-60, **kwargs):
        redirect_to_wulfric()

    def clear(self):
        redirect_to_wulfric()

    def legend(self, **kwargs):
        redirect_to_wulfric()

    def plot_unit_cell(
        self,
        lattice,
        vectors=True,
        color="#274DD1",
        label=None,
        vector_pad=1.1,
        conventional=False,
        reciprocal=False,
        normalize=False,
    ):
        redirect_to_wulfric()

    def plot_wigner_seitz(
        self,
        lattice,
        vectors=True,
        color="black",
        label=None,
        vector_pad=1.1,
        reciprocal=False,
        normalize=False,
    ):
        redirect_to_wulfric()

    def plot_kpath(self, lattice, color="black", label=None, normalize=False):
        redirect_to_wulfric()


class PlotlyBackend(AbstractBackend):

    def __init__(self, fig=None):
        redirect_to_wulfric()

    def show(self, **kwargs):
        redirect_to_wulfric()

    def save(
        self,
        output_name="lattice_graph.png",
        kwargs_update_layout=None,
        kwargs_write_html=None,
    ):
        redirect_to_wulfric()

    def plot_unit_cell(
        self,
        lattice,
        vectors=True,
        color="#274DD1",
        label=None,
        conventional=False,
        reciprocal=False,
        normalize=False,
    ):
        redirect_to_wulfric()

    def plot_wigner_seitz(
        self,
        lattice,
        vectors=True,
        label=None,
        color="black",
        reciprocal=False,
        normalize=False,
    ):
        redirect_to_wulfric()

    def plot_kpath(
        self, lattice, color="#000000", label=None, normalize=False, **kwargs
    ):
        redirect_to_wulfric()
