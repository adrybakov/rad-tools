# RAD-tools - program for spin Hamiltonian and magnons.
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: adrybakov.com
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

from typing import Iterable

__all__ = ["plot_hlines", "plot_vlines"]


def plot_hlines(ax, positions, **kwargs):
    r"""
    Plot horizontal lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of y positions.
    **kwargs
        Keyword arguments to be passed to :meth:`matplotlib.axes.Axes.hlines`.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    if "transfrom" in kwargs:
        kwargs.pop("transform")
    if "color" not in kwargs:
        kwargs["color"] = "grey"
    if "linewidths" not in kwargs:
        kwargs["linewidths"] = 0.5
    if "linestyles" not in kwargs:
        kwargs["linestyles"] = "dashed"

    ax.hlines(positions, 0, 1, transform=ax.get_yaxis_transform(), **kwargs)


def plot_vlines(ax, positions, **kwargs):
    r"""
    Plot vertical lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of x positions.
    **kwargs
        Keyword arguments to be passed to :meth:`matplotlib.axes.Axes.vlines`.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    if "transfrom" in kwargs:
        kwargs.pop("transform")
    if "color" not in kwargs:
        kwargs["color"] = "grey"
    if "linewidths" not in kwargs:
        kwargs["linewidths"] = 0.5
    if "linestyles" not in kwargs:
        kwargs["linestyles"] = "dashed"

    ax.vlines(positions, 0, 1, transform=ax.get_xaxis_transform(), **kwargs)
