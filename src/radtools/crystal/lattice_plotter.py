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

from random import choices
from string import ascii_lowercase
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d

from radtools.crystal.constants import HS_PLOT_NAMES
from radtools.geometry import volume

try:
    import plotly.graph_objects as go

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

__all__ = ["MatplotlibBackend", "PlotlyBackend"]


class AbstractBackend:
    def __init__(self) -> None:
        self.kinds = {
            "conventional": self.plot_conventional,
            "primitive": self.plot_primitive,
            "brillouin": self.plot_brillouin,
            "kpath": self.plot_kpath,
            "brillouin-kpath": self.plot_brillouin_kpath,
            "brillouin_kpath": self.plot_brillouin_kpath,
            "wigner-seitz": self.plot_wigner_seitz,
            "wigner_seitz": self.plot_wigner_seitz,
            "unit-cell": self.plot_unit_cell,
            "unit_cell": self.plot_unit_cell,
        }

    # Backend-independent functions
    def plot(self, *args, kind, **kwargs):
        r"""
        Main plotting entry point.

        Actual list of supported kinds can be check with:

        .. doctest::

            >>> self.kinds.keys() # doctest: +SKIP

        Parameters
        ----------
        kind : str or list of str
            Type of the plot to be plotted. Supported plots:

            * "conventional"
            * "primitive"
            * "brillouin"
            * "kpath"
            * "brillouin-kpath"
            * "wigner-seitz"
            * "unit-cell"
        *args
            Passed directly to the plotting functions.
        **kwargs
            Passed directly to the plotting functions.
        """
        if isinstance(kind, str):
            kinds = [kind]
        else:
            kinds = kind
        for kind in kinds:
            if kind in self.kinds:
                self.kinds[kind](*args, **kwargs)
            else:
                raise ValueError(f"Plot kind '{kind}' does not exist!")

    # Backend-dependent functions
    def remove(self, *args, **kwargs):
        raise NotImplementedError

    def show(self, *args, **kwargs):
        raise NotImplementedError

    def save(self, *args, **kwargs):
        raise NotImplementedError

    def clear(self, *args, **kwargs):
        raise NotImplementedError

    def legend(self, *args, **kwargs):
        raise NotImplementedError

    # Backend-independent functions
    def plot_brillouin(self, *args, color="#FF4D67", **kwargs):
        r"""
        Plot brillouin zone.

        Parameters
        ----------
        *args
            Passed to the :py:meth:`.plot_wigner_seitz` function.
        color : str, default "#FF4D67"
            Colour for the brillouin zone. Any format supported by the used backend.
        **kwargs
            Passed to the :py:meth:`.plot_wigner_seitz` function.

        See Also
        --------
        plot_unit_cell : for the list of parameters
        """

        self.plot_wigner_seitz(*args, reciprocal=True, color="#FF4D67", **kwargs)

    def plot_brillouin_kpath(
        self, *args, zone_color="#FF4D67", path_color="black", **kwargs
    ):
        r"""
        Plot brillouin zone and kpath.

        Parameters
        ----------
        *args
            Passed to the :py:meth:`.plot_brillouin` and :py:meth:`.plot_kpath` functions.
        zone_color : str, default "#FF4D67"
            Colour for the brillouin zone. Any format supported by the used backend.
        zone_color : str, default "black"
            Colour for the k path. Any format supported by the used backend.
        **kwargs
            Passed to the :py:meth:`.plot_brillouin` and :py:meth:`.plot_kpath` functions.

        See Also
        --------
        plot_brillouin : plot brillouin zone
        plot_kpath : plot kpath
        """

        self.plot_brillouin(*args, color=zone_color, **kwargs)
        self.plot_kpath(*args, color=path_color, **kwargs)

    def plot_primitive(self, *args, **kwargs):
        r"""
        Plot primitive unit cell.

        Parameters
        ----------
        **kwargs
            Passed to the :py:meth:`.plot_unit_cell` function.

        See Also
        --------
        plot_unit_cell : for the list of parameters
        """

        self.plot_unit_cell(*args, conventional=False, **kwargs)

    def plot_conventional(self, *args, **kwargs):
        r"""
        Plot conventional unit cell.

        See Also
        --------
        plot_unit_cell : for the list of parameters
        """

        self.plot_unit_cell(*args, conventional=True, **kwargs)

    # Backend-dependent functions
    def plot_unit_cell(self, *args, **kwargs):
        raise NotImplementedError

    def plot_wigner_seitz(self, *args, **kwargs):
        raise NotImplementedError

    def plot_kpath(self, *args, **kwargs):
        raise NotImplementedError


# Better 3D arrows, see: https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-a-3d-plot
class Arrow3D(FancyArrowPatch):
    def __init__(self, ax, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
        self.ax = ax

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.ax.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def do_3d_projection(self, *_, **__):
        return 0


class MatplotlibBackend(AbstractBackend):
    r"""
    Plotting engine for the lattice with |matplotlib|_.

    .. versionadded:: 0.8.5

    Parameters
    ----------
    fig : matplotlib figure, optional
        Figure to plot on. If not provided, a new figure and ``ax`` is created.
    ax : matplotlib axis, optional
        Axis to plot on. If not provided, a new axis is created.
    background : bool, default True
        Whether to keep the axis on the plot.
    focal_length : float, default 0.2
        See: |matplotlibFocalLength|_

    Attributes
    ----------
    fig : matplotlib figure
        Figure to plot on.
    ax : matplotlib axis
        Axis to plot on.
    artists : dict
        Dictionary of the artists. Keys are the plot kinds, values are the lists of artists.
    """

    def __init__(self, fig=None, ax=None, background=True, focal_length=0.2):
        super().__init__()
        if fig is None:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(projection="3d")
        elif ax is None:
            ax = fig.add_subplot(projection="3d")

        rcParams["axes.linewidth"] = 0
        rcParams["xtick.color"] = "#B3B3B3"
        ax.set_proj_type("persp", focal_length=focal_length)
        if background:
            ax.axes.linewidth = 0
            ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            ax.set_xlabel("x", fontsize=15, alpha=0.5)
            ax.set_ylabel("y", fontsize=15, alpha=0.5)
            ax.set_zlabel("z", fontsize=15, alpha=0.5)
            ax.tick_params(axis="both", zorder=0, color="#B3B3B3")
        else:
            ax.axis("off")
        self.fig = fig
        self.ax = ax
        self.artists = {}

    def remove(self, kind="primitive"):
        r"""
        Remove a set of artists from the plot.

        Parameters
        ----------
        kind : str or list of str
            Type of the plot to be removed. Supported kinds:

            * "conventional"
            * "primitive"
            * "brillouin"
            * "kpath"
            * "brillouin_kpath"
            * "wigner_seitz"

        ax : axis, optional
            3D matplotlib axis for the plot.
        """

        if kind == "brillouin_kpath":
            kinds = ["brillouin", "kpath"]
        else:
            kinds = [kind]

        for kind in kinds:
            if kind not in self.artists:
                raise ValueError(f"No artists for the {kind} kind.")
            for artist in self.artists[kind]:
                if isinstance(artist, list):
                    for i in artist:
                        i.remove()
                else:
                    artist.remove()
            del self.artists[kind]
        self.ax.relim(visible_only=True)
        self.ax.set_aspect("equal")

    def plot(self, lattice, kind="primitive", **kwargs):
        r"""
        Main plotting function of the Lattice.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        kind : str or list od str
            Type of the plot to be plotted. Supported plots:

            * "conventional"
            * "primitive"
            * "brillouin"
            * "kpath"
            * "brillouin-kpath"
            * "wigner-seitz"
            * "unit-cell"

        **kwargs
            Parameters to be passed to the plotting function.
            See each function for the list of supported parameters.

        Raises
        ------
        ValueError
            If the plot kind is not supported.


        See Also
        --------
        plot_conventional : "conventional" plot.
        plot_primitive : "primitive" plot.
        plot_brillouin : "brillouin" plot.
        plot_kpath : "kpath" plot.
        plot_brillouin_kpath : "brillouin_kpath" plot.
        plot_wigner_seitz : "wigner-seitz" plot.
        plot_unit_cell : "unit-cell" plot.
        show : Shows the plot.
        save : Save the figure in the file.
        """

        super().plot(lattice, kind=kind, **kwargs)
        self.ax.relim()
        self.ax.set_aspect("equal")

    def show(self, elev=30, azim=-60):
        r"""
        Show the figure in the interactive mode.

        Parameters
        ----------
        elev : float, default 30
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        azim : float, default -60
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        """
        self.ax.set_aspect("equal")
        self.ax.view_init(elev=elev, azim=azim)
        plt.show()
        self.fig = None
        self.ax = None
        plt.close()

    def save(self, output_name="lattice_graph.png", elev=30, azim=-60, **kwargs):
        r"""
        Save the figure in the file.

        .. versionchanged:: 0.8.5 Renamed from ``savefig``

        Parameters
        ----------
        output_name : str, default "lattice_graph.png"
            Name of the file to be saved. With extension.
        elev : float, default 30
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        azim : float, default -60
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        **kwargs
            Parameters to be passed to the |matplotlibSavefig|_.
        """

        self.ax.set_aspect("equal")
        self.ax.view_init(elev=elev, azim=azim)
        self.fig.savefig(output_name, **kwargs)

    def clear(self):
        r"""
        Clear the axis.
        """

        if self.ax is not None:
            self.ax.cla()

    def legend(self, **kwargs):
        r"""
        Add legend to the figure.

        Parameters
        ----------
        **kwargs :
            Directly passed to the |matplotlibLegend|_.
        """

        self.ax.legend(**kwargs)

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
        r"""
        Plot real or reciprocal space unit cell.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        vectors : bool, default True
            Whether to plot lattice vectors.
        color : str, default "#274DD1"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        label : str, optional
            Label for the plot.
        vector_pad : float, default 1.1
            Multiplier for the position of the vectors labels. 1 = position of the vector.
        conventional : bool, default False
            Whether to plot conventional cell. Affects result only for the
            Bravais lattice classes. Ignored for the general :py:class:`.Lattice`.
            Only primitive unit cell is supported for reciprocal space.
        reciprocal : bool, default False
            Whether to plot reciprocal or real unit cell.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """
        if reciprocal and conventional:
            raise ValueError("Conventional cell is not supported in reciprocal space.")
        if conventional:
            artist_group = "conventional"
        else:
            artist_group = "primitive"

        if reciprocal:
            artist_group += "_reciprocal"
            vector_label = "b"
        else:
            artist_group += "_real"
            vector_label = "a"

        self.artists[artist_group] = []

        if conventional:
            cell = lattice.conv_cell
        elif reciprocal:
            cell = lattice.reciprocal_cell
        else:
            cell = lattice.cell

        if normalize:
            cell /= abs(volume(cell) ** (1 / 3.0))

        if label is not None:
            self.artists[artist_group].append(
                self.ax.scatter(0, 0, 0, color=color, label=label)
            )
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            for i in range(3):
                self.artists[artist_group].append(
                    self.ax.text(
                        cell[i][0] * vector_pad[i],
                        cell[i][1] * vector_pad[i],
                        cell[i][2] * vector_pad[i],
                        f"${vector_label}_{i+1}$",
                        fontsize=20,
                        color=color,
                        ha="center",
                        va="center",
                    )
                )
                # Try beautiful arrows
                try:
                    self.artists[artist_group].append(
                        self.ax.add_artist(
                            Arrow3D(
                                self.ax,
                                [0, cell[i][0]],
                                [0, cell[i][1]],
                                [0, cell[i][2]],
                                mutation_scale=20,
                                arrowstyle="-|>",
                                color=color,
                                lw=2,
                                alpha=0.7,
                            )
                        )
                    )
                # Go to default
                except:
                    self.artists[artist_group].append(
                        self.ax.quiver(
                            0,
                            0,
                            0,
                            *tuple(cell[i]),
                            arrow_length_ratio=0.2,
                            color=color,
                            alpha=0.7,
                            linewidth=2,
                        )
                    )
                # Ghost point to account for the plot range
                self.artists[artist_group].append(self.ax.scatter(*tuple(cell[i]), s=0))

        def plot_line(line, shift):
            self.artists[artist_group].append(
                self.ax.plot(
                    [shift[0], shift[0] + line[0]],
                    [shift[1], shift[1] + line[1]],
                    [shift[2], shift[2] + line[2]],
                    color=color,
                )
            )

        for i in range(0, 3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            plot_line(cell[i], np.zeros(3))
            plot_line(cell[i], cell[j])
            plot_line(cell[i], cell[k])
            plot_line(cell[i], cell[j] + cell[k])

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
        r"""
        Plot Wigner-Seitz unit cell.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        vectors : bool, default True
            Whether to plot lattice vectors.
        color : str, default "black"
            Colour for the plot. Any format supported by the used backend.
        label : str, optional
            Label for the plot.
        vector_pad : float, default 1.1
            Multiplier for the position of the vectors labels. 1 = position of the vector.
        reciprocal : bool, default False
            Whether to plot reciprocal or real Wigner-Seitz cell.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        if reciprocal:
            artist_group = "brillouin"
        else:
            artist_group = "wigner_seitz"

        self.artists[artist_group] = []

        if reciprocal:
            v1, v2, v3 = lattice.b1, lattice.b2, lattice.b3
            vector_label = "b"
        else:
            v1, v2, v3 = lattice.a1, lattice.a2, lattice.a3
            vector_label = "a"

        if color is None:
            color = "black"

        if normalize:
            factor = volume(v1, v2, v3) ** (1 / 3.0)
            v1 /= factor
            v2 /= factor
            v3 /= factor

        vs = [v1, v2, v3]

        if label is not None:
            self.artists[artist_group].append(
                self.ax.scatter(0, 0, 0, color=color, label=label)
            )
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            for i in range(3):
                self.artists[artist_group].append(
                    self.ax.text(
                        vs[i][0] * vector_pad[i],
                        vs[i][1] * vector_pad[i],
                        vs[i][2] * vector_pad[i],
                        f"${vector_label}_{i+1}$",
                        fontsize=20,
                        color=color,
                        ha="center",
                        va="center",
                    )
                )
                # Try beautiful arrows
                try:
                    self.artists[artist_group].append(
                        self.ax.add_artist(
                            Arrow3D(
                                self.ax,
                                [0, vs[i][0]],
                                [0, vs[i][1]],
                                [0, vs[i][2]],
                                mutation_scale=20,
                                arrowstyle="-|>",
                                color=color,
                                lw=2,
                                alpha=0.8,
                            )
                        )
                    )
                # Go to default
                except:
                    self.artists[artist_group].append(
                        self.ax.quiver(
                            0,
                            0,
                            0,
                            *tuple(vs[i]),
                            arrow_length_ratio=0.2,
                            color=color,
                            alpha=0.5,
                        )
                    )
                # Ghost point to account for the plot range
                self.artists[artist_group].append(self.ax.scatter(*tuple(vs[i]), s=0))

        edges, _ = lattice.voronoi_cell(reciprocal=reciprocal, normalize=normalize)
        for p1, p2 in edges:
            self.artists[artist_group].append(
                self.ax.plot(
                    [p1[0], p2[0]],
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    color=color,
                )
            )

    def plot_kpath(self, lattice, color="black", label=None, normalize=False):
        r"""
        Plot k path in the reciprocal space.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        color : str, default "black"
            Colour for the plot. Any format supported by the used backend.
        label : str, optional
            Label for the plot.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        artist_group = "kpath"

        self.artists[artist_group] = []

        cell = lattice.reciprocal_cell

        kp = lattice.kpoints

        if normalize:
            cell /= volume(cell) ** (1 / 3.0)

        for point in kp.hs_names:
            self.artists[artist_group].append(
                self.ax.scatter(
                    *tuple(kp.hs_coordinates[point] @ cell),
                    s=36,
                    color=color,
                )
            )
            if point == "S" and lattice.type() == "BCT":
                label = "$\\Sigma$"
            elif point == "S1" and lattice.type() == "BCT":
                label = "$\\Sigma_1$"
            else:
                label = HS_PLOT_NAMES[point]
            self.artists[artist_group].append(
                self.ax.text(
                    *tuple(
                        kp.hs_coordinates[point] @ cell
                        + 0.025 * cell[0]
                        + +0.025 * cell[1]
                        + 0.025 * cell[2]
                    ),
                    label,
                    fontsize=20,
                    color=color,
                )
            )
        if label is not None:
            self.artists[artist_group].append(
                self.ax.scatter(
                    0,
                    0,
                    0,
                    s=36,
                    color=color,
                    label=label,
                )
            )

        for subpath in kp.path:
            for i in range(len(subpath) - 1):
                self.artists[artist_group].append(
                    self.ax.plot(
                        *tuple(
                            np.concatenate(
                                (
                                    kp.hs_coordinates[subpath[i]] @ cell,
                                    kp.hs_coordinates[subpath[i + 1]] @ cell,
                                )
                            )
                            .reshape(2, 3)
                            .T
                        ),
                        color=color,
                        alpha=0.5,
                        linewidth=3,
                    )
                )


class PlotlyBackend(AbstractBackend):
    r"""
    Plotting engine for the lattice with |plotly|_.

    .. versionadded:: 0.8.5

    Parameters
    ----------
    fig : plotly graph object
        Figure to plot on. If not provided, a new figure is created.

    Attributes
    ----------
    fig : plotly graph object
        Figure to plot on.

    Notes
    -----
    Plotly is not in the dependencies of the package.
    To use this backend, install it with ``pip``:

    .. code-block:: bash

        pip install plotly

    or ``pip3``

    .. code-block:: bash

        pip3 install plotly

    """

    def __init__(self, fig=None):
        if not PLOTLY_AVAILABLE:
            raise ImportError(
                'Plotly is not available. Install it with "pip install plotly"'
            )
        super().__init__()
        if fig is None:
            fig = go.Figure()
        self.fig = fig

    def show(self, **kwargs):
        r"""
        Show the figure in the interactive mode.

        Parameters
        ----------
        **kwargs
            Passed directly to the |plotly-update-layout|_.
        """

        # Set up defaults
        if "width" not in kwargs:
            kwargs["width"] = 800
        if "height" not in kwargs:
            kwargs["height"] = 700
        if "yaxis_scaleanchor" not in kwargs:
            kwargs["yaxis_scaleanchor"] = "x"
        if "showlegend" not in kwargs:
            kwargs["showlegend"] = False
        if "autosize" not in kwargs:
            kwargs["autosize"] = False

        self.fig.update_layout(**kwargs)
        self.fig.show()

    def save(
        self,
        output_name="lattice_graph.png",
        kwargs_update_layout=None,
        kwargs_write_html=None,
    ):
        r"""
        Save the figure in the html file.

        Parameters
        ----------
        output_name : str, default "lattice_graph.png"
            Name of the file to be saved. With extension.
        kwargs_update_layout : dict, optional
            Passed directly to the |plotly-update-layout|_.
        kwargs_write_html : dict, optional
            Passed directly to the |plotly-write-html|_.
        """

        if kwargs_update_layout is None:
            kwargs_update_layout = {}
        if kwargs_write_html is None:
            kwargs_write_html = {}

        self.fig.update_scenes(aspectmode="data")

        self.fig.update_layout(**kwargs_update_layout)
        self.fig.write_html(output_name, **kwargs_write_html)

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
        r"""
        Plot real or reciprocal space unit cell.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        vectors : bool, default True
            Whether to plot lattice vectors.
        color : str, default "#274DD1"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        label : str, optional
            Label for the plot.
        conventional : bool, default False
            Whether to plot conventional cell. Affects result only for the
            Bravais lattice classes. Ignored for the general :py:class:`.Lattice`.
            Only primitive unit cell is supported for reciprocal space.
        reciprocal : bool, default False
            Whether to plot reciprocal or real unit cell.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """
        if reciprocal and conventional:
            raise ValueError("Conventional cell is not supported in reciprocal space.")
        if conventional:
            artist_group = "conventional"
        else:
            artist_group = "primitive"

        if reciprocal:
            artist_group += "_reciprocal"
            vector_label = "b"
        else:
            artist_group += "_real"
            vector_label = "a"

        if conventional:
            cell = lattice.conv_cell
        elif reciprocal:
            cell = lattice.reciprocal_cell
        else:
            cell = lattice.cell

        if normalize:
            cell /= abs(volume(cell) ** (1 / 3.0))

        legendgroup = "".join(choices(ascii_lowercase, k=10))

        if vectors:
            labels = [f"{vector_label}{i+1}" for i in range(3)]
            for i in range(3):
                x = [0, cell[i][0]]
                y = [0, cell[i][1]]
                z = [0, cell[i][2]]
                self.fig.add_traces(
                    data=[
                        {
                            "x": x,
                            "y": y,
                            "z": z,
                            "mode": "lines",
                            "type": "scatter3d",
                            "hoverinfo": "none",
                            "line": {"color": color, "width": 3},
                            "showlegend": False,
                            "legendgroup": legendgroup,
                        },
                        {
                            "type": "cone",
                            "x": [x[1]],
                            "y": [y[1]],
                            "z": [z[1]],
                            "u": [0.2 * (x[1] - x[0])],
                            "v": [0.2 * (y[1] - y[0])],
                            "w": [0.2 * (z[1] - z[0])],
                            "anchor": "tip",
                            "hoverinfo": "none",
                            "colorscale": [[0, color], [1, color]],
                            "showscale": False,
                            "showlegend": False,
                            "legendgroup": legendgroup,
                        },
                    ]
                )
                self.fig.add_traces(
                    data=go.Scatter3d(
                        mode="text",
                        x=[1.2 * x[1]],
                        y=[1.2 * y[1]],
                        z=[1.2 * z[1]],
                        marker=dict(size=0, color=color),
                        text=labels[i],
                        hoverinfo="none",
                        textposition="top center",
                        textfont=dict(size=12),
                        showlegend=False,
                        legendgroup=legendgroup,
                    )
                )

        def plot_line(line, shift, showlegend=False):
            self.fig.add_traces(
                data=go.Scatter3d(
                    mode="lines",
                    x=[shift[0], shift[0] + line[0]],
                    y=[shift[1], shift[1] + line[1]],
                    z=[shift[2], shift[2] + line[2]],
                    line=dict(color=color),
                    hoverinfo="none",
                    legendgroup=legendgroup,
                    name=label,
                    showlegend=showlegend,
                ),
            )

        showlegend = label is not None
        for i in range(0, 3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            plot_line(cell[i], np.zeros(3), showlegend=showlegend)
            if showlegend:
                showlegend = False
            plot_line(cell[i], cell[j])
            plot_line(cell[i], cell[k])
            plot_line(cell[i], cell[j] + cell[k])

    def plot_wigner_seitz(
        self,
        lattice,
        vectors=True,
        label=None,
        color="black",
        reciprocal=False,
        normalize=False,
    ):
        r"""
        Plot Wigner-Seitz unit cell.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        vectors : bool, default True
            Whether to plot lattice vectors.
        label : str, optional
            Label for the plot.
        color : str, default "black" or "#FF4D67"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        reciprocal : bool, default False
            Whether to plot reciprocal or real Wigner-Seitz cell.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        if reciprocal:
            v1, v2, v3 = lattice.b1, lattice.b2, lattice.b3
            vector_label = "b"
        else:
            v1, v2, v3 = lattice.a1, lattice.a2, lattice.a3
            vector_label = "a"

        if color is None:
            color = "black"

        if normalize:
            factor = volume(v1, v2, v3) ** (1 / 3.0)
            v1 /= factor
            v2 /= factor
            v3 /= factor

        vs = [v1, v2, v3]

        legendgroup = "".join(choices(ascii_lowercase, k=10))

        if vectors:
            labels = [f"{vector_label}{i+1}" for i in range(3)]
            for i in range(3):
                x = [0, vs[i][0]]
                y = [0, vs[i][1]]
                z = [0, vs[i][2]]
                self.fig.add_traces(
                    data=[
                        {
                            "x": x,
                            "y": y,
                            "z": z,
                            "mode": "lines",
                            "type": "scatter3d",
                            "hoverinfo": "none",
                            "line": {"color": color, "width": 3},
                            "showlegend": False,
                            "legendgroup": legendgroup,
                        },
                        {
                            "type": "cone",
                            "x": [x[1]],
                            "y": [y[1]],
                            "z": [z[1]],
                            "u": [0.2 * (x[1] - x[0])],
                            "v": [0.2 * (y[1] - y[0])],
                            "w": [0.2 * (z[1] - z[0])],
                            "anchor": "tip",
                            "hoverinfo": "none",
                            "colorscale": [[0, color], [1, color]],
                            "showscale": False,
                            "showlegend": False,
                            "legendgroup": legendgroup,
                        },
                    ]
                )
                self.fig.add_traces(
                    data=go.Scatter3d(
                        mode="text",
                        x=[1.2 * x[1]],
                        y=[1.2 * y[1]],
                        z=[1.2 * z[1]],
                        marker=dict(size=0, color=color),
                        text=labels[i],
                        hoverinfo="none",
                        textposition="top center",
                        textfont=dict(size=12),
                        showlegend=False,
                        legendgroup=legendgroup,
                    )
                )

        edges, _ = lattice.voronoi_cell(reciprocal=reciprocal, normalize=normalize)
        showlegend = label is not None
        for p1, p2 in edges:
            xyz = np.array([p1, p2]).T
            self.fig.add_traces(
                data=go.Scatter3d(
                    mode="lines",
                    x=xyz[0],
                    y=xyz[1],
                    z=xyz[2],
                    line=dict(color=color),
                    hoverinfo="none",
                    showlegend=showlegend,
                    legendgroup=legendgroup,
                    name=label,
                ),
            )
            if showlegend:
                showlegend = False

    def plot_kpath(
        self, lattice, color="#000000", label=None, normalize=False, **kwargs
    ):
        r"""
        Plot k path in the reciprocal space.

        Parameters
        ----------
        lattice : :py:class:`.Lattice`
            Lattice to be plotted.
        color : str, default "#000000"
            Colour for the plot. Any format supported by the used backend.
        label : str, optional
            Label for the plot.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """
        cell = lattice.reciprocal_cell

        kp = lattice.kpoints

        if normalize:
            cell /= volume(cell) ** (1 / 3.0)

        p_abs = []
        p_rel = []
        labels = []
        for point in kp.hs_names:
            p_abs.append(tuple(kp.hs_coordinates[point] @ cell))
            p_rel.append(kp.hs_coordinates[point])
            if point == "S" and lattice.type() == "BCT":
                p_label = R"$\\Sigma$"
            elif point == "S1" and lattice.type() == "BCT":
                p_label = R"$\\Sigma_1$"
            else:
                p_label = HS_PLOT_NAMES[point]
            labels.append(point)

        p_abs = np.array(p_abs).T

        self.fig.add_traces(
            data=go.Scatter3d(
                mode="markers+text",
                x=p_abs[0],
                y=p_abs[1],
                z=p_abs[2],
                marker=dict(size=6, color=color),
                text=labels,
                hoverinfo="text",
                hovertext=p_rel,
                textposition="top center",
                textfont=dict(size=16),
                showlegend=False,
            )
        )

        for subpath in kp.path:
            xyz = []
            for i in range(len(subpath)):
                xyz.append(kp.hs_coordinates[subpath[i]] @ cell)

            xyz = np.array(xyz).T
            self.fig.add_traces(
                data=go.Scatter3d(
                    mode="lines",
                    x=xyz[0],
                    y=xyz[1],
                    z=xyz[2],
                    line=dict(color=color),
                    hoverinfo="none",
                    showlegend=False,
                ),
            )
