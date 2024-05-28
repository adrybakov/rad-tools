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

from matplotlib.colors import LinearSegmentedColormap, to_rgb

__all__ = ["custom_cmap"]


def custom_cmap(start_color, finish_color):
    r"""
    Prepare custom colormap. From one color to the other.

    Parameters
    ----------
    start_color : Color
        Start color.
    finish_color : Color
        Finish color.
    """

    r1, g1, b1 = to_rgb(start_color)

    r2, g2, b2 = to_rgb(finish_color)

    cdict = {
        "red": ((0, r1, r1), (1, r2, r2)),
        "green": ((0, g1, g1), (1, g2, g2)),
        "blue": ((0, b1, b1), (1, b2, b2)),
    }

    return LinearSegmentedColormap("custom_cmap", cdict, N=256)
