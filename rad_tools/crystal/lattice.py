r"""
General 3D lattice.
"""

from math import acos, cos, floor, log10, pi, sqrt, sin
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d
from scipy.spatial import Voronoi
from termcolor import cprint

from rad_tools.routines import _todegrees, _toradians, angle, volume, reciprocal_cell

__all__ = ["Lattice", "get_niggli"]


# Better 3D arrows, see: https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-a-3d-plot
class Arrow3D(FancyArrowPatch):
    def __init__(self, ax, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
        self._ax = ax

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self._ax.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def do_3d_projection(self, *_, **__):
        return 0


class Lattice:
    r"""
    General 3D lattice.

    In absence of the atoms (which is always the case for the lattice)
    and additional lattice points. Every cell is primitive, since lattice points are
    constructed form the translations. Therefore, general Lattice class does not
    dos not distinguishes between primitive and conventional lattice as in the
    Bravais Lattice. Therefore only cell attribute is present and it is always
    interpreted as the primitive unit cell. In the case of the Bravais lattice additional
    attribute conv_cell appears.

    Lattice can be created in a three alternative ways:

    .. doctest::

        >>> import rad_tools as rad
        >>> l = rad.Lattice(cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> l = rad.Lattice(a1=[1,0,0], a2=[0,1,0], a3=[0,0,1]])
        >>> l = rad.Lattice(1, 1, 1, 90, 90, 90)

    Parameters
    ----------
    cell : (3,3) |array_like|_
        Unit cell, rows are vectors, columns are coordinates.
    a1 : (3,) |array_like|_
        First vector of unit cell (cell[0]).
    a2 : (3,) |array_like|_
        SEcond vector of unit cell (cell[1]).
    a3 : (3,) |array_like|_
        Third vector of unit cell (cell[2]).
    a : float, default=1
        Length of the :math:`a_1` vector.
    b : float, default=1
        Length of the :math:`a_2` vector.
    c : float, default=1
        Length of the :math:`a_3` vector.
    alpha : float, default=90
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float, default=90
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float, default=90
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.

    Attributes
    ----------
    path : list
        K-path.
    kpoints : dist
        Dictionary of the high symmetry points.
        Coordinates are given in relative coordinates.

        .. code-block:: python

            kpoints = {"Name" : [k_x, k_y, k_z], ...}

    pearson_symbol : str
    crystal_family : str
    centring_type : str
    """

    _pearson_symbol = None

    def __init__(self, *args) -> None:
        self._cell = None
        if len(args) == 1:
            self.cell = np.array(args[0])
        elif len(args) == 3:
            self.cell = np.array(args)
        elif len(args) == 6:
            a, b, c, alpha, beta, gamma = args
            alpha = alpha * _toradians
            beta = beta * _toradians
            gamma = gamma * _toradians
            self.cell = np.array(
                [
                    [a, 0, 0],
                    [b * cos(gamma), b * sin(gamma), 0],
                    [
                        c * cos(beta),
                        c / sin(gamma) * (cos(alpha) - cos(beta) * cos(gamma)),
                        c
                        / sin(gamma)
                        * sqrt(
                            1
                            + 2 * cos(alpha) * cos(beta) * cos(gamma)
                            - cos(alpha) ** 2
                            - cos(beta) ** 2
                            - cos(gamma) ** 2
                        ),
                    ],
                ]
            )
        else:
            raise ValueError(
                "Unable to identify input parameters. "
                + "Supported: one (3,3) array_like, or three (3,) array_like, or 6 floats."
            )
        self.points = {}
        self._path = None
        self._default_path = None
        self._fig = None
        self._ax = None
        self._artists = {}
        self._PLOT_NAMES = {
            "G": "$\\Gamma$",
            "M": "$M$",
            "R": "$R$",
            "X": "$X$",
            "K": "$K$",
            "L": "$L$",
            "U": "$U$",
            "W": "$W$",
            "H": "$H$",
            "P": "$P$",
            "N": "$N$",
            "A": "$A$",
            "Z": "$Z$",
            "Z1": "$Z_1$",
            "Y": "$Y$",
            "Y1": "$Y_1$",
            "S": "$S$",  # it is overwritten to sigma if needed.
            "S1": "$S_1$",  # it is overwritten to sigma if needed.
            "T": "$T$",
            "A1": "$A_1$",
            "X1": "$X_1$",
            "C": "$C$",
            "C1": "$C_1$",
            "D": "$D$",
            "D1": "$D_1$",
            "H1": "$H_1$",
            "L1": "$L_1$",
            "L2": "$L_2$",
            "B": "$B$",
            "B1": "$B_1$",
            "F": "$F$",
            "P1": "$P_1$",
            "P2": "$P_2$",
            "Q": "$Q$",
            "Q1": "$Q_1$",
            "E": "$E$",
            "H2": "$H_2$",
            "M1": "$M_1$",
            "M2": "$M_2$",
            "N1": "$N_1$",
            "F1": "$F_1$",
            "F2": "$F_2$",
            "F3": "$F_3$",
            "I": "$I$",
            "I1": "$I_1$",
            "X2": "$X_2$",
            "Y2": "$Y_2$",
            "Y3": "$Y_3$",
        }

    @property
    def cell(self):
        if self._cell is None:
            raise AttributeError(f"Cell is not defined for lattice {self}")
        return self._cell

    @cell.setter
    def cell(self, new_cell):
        try:
            new_cell = np.array(new_cell)
        except:
            raise ValueError(f"New cell is not array_like: {new_cell}")
        if new_cell.shape != (3, 3):
            raise ValueError(f"New cell is not 3 x 3 matrix.")
        self._cell = new_cell

    @property
    def path(self):
        r"""
        K-point path.

        Manage default path for the predefined lattices and custom user-defined path.
        """
        if self._path is None and self._default_path is None:
            return []
        if self._path is None:
            return self._default_path
        return self._path

    @property
    def pearson_symbol(self):
        r"""
        Pearson symbol.

        Returns
        -------
        pearson_symbol : str
            Pearson symbol of the lattice.

        Raises
        ------
        RuntimeError
            If the type of the lattice is not defined.

        Notes
        -----
        See: |PearsonSymbol|_
        """

        if self._pearson_symbol is not None:
            return self._pearson_symbol
        raise RuntimeError("Type of the lattice is not defined.")

    @property
    def crystal_family(self):
        r"""
        Crystal family.

        Returns
        -------
        crystal_family : str
            Crystal family of the lattice.

        Raises
        ------
        ValueError
            If the type of the lattice is not defined.

        Notes
        -----
        See: |PearsonSymbol|_
        """

        return self.pearson_symbol[0]

    @property
    def centring_type(self):
        r"""
        Centring type.

        Returns
        -------
        centring_type : str
            Centring type of the lattice.

        Raises
        ------
        ValueError
            If the type of the lattice is not defined.

        Notes
        -----
        See: |PearsonSymbol|_
        """

        return self.pearson_symbol[1]

    @property
    def a1(self):
        r"""
        First lattice vector :math:`\vec{a}_1`.
        """
        return self.cell[0]

    @property
    def a2(self):
        r"""
        Second lattice vector :math:`\vec{a}_2`.
        """
        return self.cell[1]

    @property
    def a3(self):
        r"""
        Third lattice vector :math:`\vec{a}_3`.
        """
        return self.cell[2]

    @property
    def a(self):
        r"""
        Length of the first lattice vector :math:`\vert\vec{a}_1\vert`.
        """

        return np.linalg.norm(self.cell[0])

    @property
    def b(self):
        r"""
        Length of the second lattice vector :math:`\vert\vec{a}_2\vert`.
        """

        return np.linalg.norm(self.cell[1])

    @property
    def c(self):
        r"""
        Length of the third lattice vector :math:`\vert\vec{a}_3\vert`.
        """

        return np.linalg.norm(self.cell[2])

    @property
    def alpha(self):
        r"""
        Angle between second and third lattice vector.

        Returns
        -------
        angle : float
            In degrees
        """

        return angle(self.a2, self.a3)

    @property
    def beta(self):
        r"""
        Angle between first and third lattice vector.

        Returns
        -------
        angle : float
            In degrees
        """

        return angle(self.a1, self.a3)

    @property
    def gamma(self):
        r"""
        Angle between first and second lattice vector.

        Returns
        -------
        angle : float
            In degrees
        """

        return angle(self.a1, self.a2)

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.
        """

        return volume(self.a1, self.a2, self.a3)

    @property
    def reciprocal_cell(self):
        r"""
        Reciprocal cell. Always primitive.
        """

        return reciprocal_cell(self.cell)

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector.

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{a}_2\times\vec{a}_3

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`
        """

        return self.reciprocal_cell[0]

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector.

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{a}_3\times\vec{a}_1

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`
        """

        return self.reciprocal_cell[1]

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector.

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}_1\times\vec{a}_2

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`
        """

        return self.reciprocal_cell[2]

    @property
    def k_a(self):
        r"""
        Length of the first reciprocal lattice vector :math:`\vert\vec{b}_1\vert`.
        """

        return np.linalg.norm(self.b1)

    @property
    def k_b(self):
        r"""
        Length of the second reciprocal lattice vector :math:`\vert\vec{b}_2\vert`.
        """

        return np.linalg.norm(self.b2)

    @property
    def k_c(self):
        r"""
        Length of the third reciprocal lattice vector :math:`\vert\vec{b}_3\vert`.
        """

        return np.linalg.norm(self.b3)

    @property
    def k_alpha(self):
        r"""
        Angle between second and third reciprocal lattice vector.
        """

        return angle(self.b2, self.b3)

    @property
    def k_beta(self):
        r"""
        Angle between first and third reciprocal lattice vector.
        """

        return angle(self.b1, self.b3)

    @property
    def k_gamma(self):
        r"""
        Angle between first and second reciprocal lattice vector.
        """

        return angle(self.b1, self.b2)

    @property
    def reciprocal_cell_volume(self):
        r"""
        Volume of the reciprocal cell.

        .. math::

            V = \vec{b}_1\cdot(\vec{b}_2\times\vec{b}_3)
        """

        return volume(self.b1, self.b2, self.b3)

    @property
    def variation(self):
        r"""There is no variations of the lattice"""
        return self.__class__.__name__

    def lattice_points(self, relative=False, reciprocal=False, normalize=False):
        r"""
        Compute lattice points

        Parameters
        ----------
        relative : bool, default False
            Whether to return relative or absolute coordinates.
        reciprocal : bool, default False
            Whether to use reciprocal or real cell.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        if reciprocal:
            cell = self.reciprocal_cell
        else:
            cell = self.cell

        if normalize:
            cell /= volume(cell) ** (1 / 3.0)

        lattice_points = np.zeros((27, 3), dtype=float)
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    point = np.array([i, j, k])
                    if not relative:
                        point = point @ cell
                    lattice_points[9 * (i + 1) + 3 * (j + 1) + (k + 1)] = point
        return lattice_points

    def voronoi_cell(self, reciprocal=False, normalize=False):
        r"""
        Computes Voronoy edges around (0,0,0) point.

        Parameters
        ----------
        reciprocal : bool, default False
            Whether to use reciprocal or real cell.

        Returns
        -------
        edges : (N, 2, 3) :numpy:`ndarray`
            N edges of the Voronoi cell around (0,0,0) point.
            Each elements contains two vectors of the points
            of the voronoi vertices forming an edge.
        vertices : (M, 3) :numpy:`ndarray`
            M vertices of the Voronoi cell around (0,0,0) point.
            Each element is a vector :math:`v = (v_x, v_y, v_z)`.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        voronoi = Voronoi(
            self.lattice_points(
                relative=False, reciprocal=reciprocal, normalize=normalize
            )
        )
        edges_index = set()
        # Thanks ase for the idea. 13 - is the index of (0,0,0) point.
        for rv, rp in zip(voronoi.ridge_vertices, voronoi.ridge_points):
            if -1 not in rv and 13 in rp:
                for j in range(0, len(rv)):
                    if (rv[j - 1], rv[j]) not in edges_index and (
                        rv[j],
                        rv[j - 1],
                    ) not in edges_index:
                        edges_index.add((rv[j - 1], rv[j]))
        edges_index = np.array(list(edges_index))
        edges = np.zeros((edges_index.shape[0], 2, 3), dtype=voronoi.vertices.dtype)
        for i in range(edges_index.shape[0]):
            edges[i][0] = voronoi.vertices[edges_index[i][0]]
            edges[i][1] = voronoi.vertices[edges_index[i][1]]
        return edges, voronoi.vertices[np.unique(edges_index.flatten())]

    def prepare_figure(self, background=True, focal_length=0.2) -> None:
        r"""
        Prepare style of the figure for the plot.

        Parameters
        ----------
        background : bool, default True
            Whether to keep the axis on the plot.
        focal_length : float, default 0.2
            See: |matplotlibFocalLength|_
        """

        if self._fig is None:
            self._fig = plt.figure(figsize=(6, 6))
            rcParams["axes.linewidth"] = 0
            rcParams["xtick.color"] = "#B3B3B3"
            self._ax = self._fig.add_subplot(projection="3d")
            self._ax.set_proj_type("persp", focal_length=focal_length)
            if background:
                self._ax.axes.linewidth = 0
                self._ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
                self._ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
                self._ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
                self._ax.set_xlabel("x", fontsize=15, alpha=0.5)
                self._ax.set_ylabel("y", fontsize=15, alpha=0.5)
                self._ax.set_zlabel("z", fontsize=15, alpha=0.5)
                self._ax.tick_params(axis="both", zorder=0, color="#B3B3B3")
            else:
                self._ax.axis("off")
            self._ax.set_aspect("equal")

    def plot(self, kind="primitive", ax=None, **kwargs):
        r"""
        Main plotting function of the Lattice.

        Parameters
        ----------
        kind : str or list od str
            Type of the plot to be plotted. Supported plots:

            * "conventional"
            * "primitive"
            * "brillouin"
            * "kpath"
            * "brillouin_kpath"
            * "wigner_seitz"

        ax : axis, optional
            3D matplotlib axis for the plot.
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
        plot_wigner_seitz : "wigner_seitz" plot.
        show : Shows the plot.
        savefig : Save the figure in the file.
        """

        if ax is None:
            self.prepare_figure()
            ax = self._ax
        if isinstance(kind, str):
            kinds = [kind]
        else:
            kinds = kind
        try:
            for kind in kinds:
                getattr(self, f"plot_{kind}")(ax=ax, **kwargs)
        except AttributeError:
            raise ValueError(f"Plot kind '{kind}' does not exist!")

    def remove(self, kind="primitive", ax=None):
        r"""
        Remove a set of artists from the plot.

        Parameters
        ----------
        kind : str or list of str
            Type of the plot to be removed. Supported plots:

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
            if kind not in self._artists:
                raise ValueError(f"No artists for the {kind} kind.")
            for artist in self._artists[kind]:
                if isinstance(artist, list):
                    for i in artist:
                        i.remove()
                else:
                    artist.remove()
            del self._artists[kind]

    def show(self):
        r"""
        Show the figure in the interactive matplotlib window.
        """
        self._ax.set_aspect("equal")
        plt.show()
        del self._fig
        del self._ax
        plt.close()

    def savefig(self, output_name="lattice_graph.png", elev=30, azim=-60, **kwargs):
        r"""
        Save the figure in the file

        Parameters
        ----------
        output_name : str, default "lattice_graph.png"
            Name of the file to be saved.
        elev : float, default 30
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        azim : float, default -60
            Passed directly to matplotlib. See |matplotlibViewInit|_.
        **kwargs
            Parameters to be passed to the |matplotlibSavefig|_.
        """

        self._ax.set_aspect("equal")
        self._ax.view_init(elev=elev, azim=azim)
        self._fig.savefig(output_name, **kwargs)

    def legend(self, **kwargs):
        r"""
        Add legend to the figure.
        Directly passed to the |matplotlibLegend|_.
        """

        self._ax.legend(**kwargs)

    def plot_real_space(
        self,
        ax=None,
        vectors=True,
        colour="#274DD1",
        label=None,
        vector_pad=1.1,
        conventional=False,
        normalize=False,
    ):
        r"""
        Plot real space unit cell.

        ax : axis, optional
            3D matplotlib axis for the plot.
        vectors : bool, default True
            Whether to plot lattice vectors.
        colour : str, default "#274DD1"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        label : str, optional
            Label for the plot.
        vector_pad : float, default 1.1
            Multiplier for the position of the vectors labels. 1 = position of the vector.
        conventional : bool, default False
            Whether to plot conventional cell. Affects result only for the
            Bravais lattice classes. Ignored for the general :py:class:`.Lattice`.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        if conventional:
            artist_group = "conventional"
        else:
            artist_group = "primitive"

        self._artists[artist_group] = []

        if ax is None:
            self.prepare_figure()
            ax = self._ax

        if conventional:
            try:
                cell = self.conv_cell
            except AttributeError:
                cell = self.cell
        else:
            cell = self.cell

        if normalize:
            cell /= volume(cell) ** (1 / 3.0)

        if label is not None:
            self._artists[artist_group].append(
                ax.scatter(0, 0, 0, color=colour, label=label)
            )
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            self._artists[artist_group].append(
                ax.text(
                    cell[0][0] * vector_pad[0],
                    cell[0][1] * vector_pad[0],
                    cell[0][2] * vector_pad[0],
                    "$a_1$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            self._artists[artist_group].append(
                ax.text(
                    cell[1][0] * vector_pad[2],
                    cell[1][1] * vector_pad[2],
                    cell[1][2] * vector_pad[2],
                    "$a_2$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            self._artists[artist_group].append(
                ax.text(
                    cell[2][0] * vector_pad[2],
                    cell[2][1] * vector_pad[2],
                    cell[2][2] * vector_pad[2],
                    "$a_3$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            for i in cell:
                # Try beautiful arrows
                try:
                    self._artists[artist_group].append(
                        ax.add_artist(
                            Arrow3D(
                                ax,
                                [0, i[0]],
                                [0, i[1]],
                                [0, i[2]],
                                mutation_scale=20,
                                arrowstyle="-|>",
                                color=colour,
                                lw=2,
                                alpha=0.7,
                            )
                        )
                    )
                # Go to default
                except:
                    self._artists[artist_group].append(
                        ax.quiver(
                            0,
                            0,
                            0,
                            *tuple(i),
                            arrow_length_ratio=0.2,
                            color=colour,
                            alpha=0.7,
                            linewidth=2,
                        )
                    )
                # Ghost point to account for the plot range
                self._artists[artist_group].append(ax.scatter(*tuple(i), s=0))

        def plot_line(line, shift):
            self._artists[artist_group].append(
                ax.plot(
                    [shift[0], shift[0] + line[0]],
                    [shift[1], shift[1] + line[1]],
                    [shift[2], shift[2] + line[2]],
                    color=colour,
                )
            )

        for i in range(0, 3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            plot_line(cell[i], np.zeros(3))
            plot_line(cell[i], cell[j])
            plot_line(cell[i], cell[k])
            plot_line(cell[i], cell[j] + cell[k])

    def plot_brillouin(
        self,
        ax=None,
        vectors=True,
        colour="#FF4D67",
        label=None,
        vector_pad=1.1,
        normalize=False,
    ):
        r"""
        Plot brillouin zone.

        ax : axis, optional
            3D matplotlib axis for the plot.
        vectors : bool, default True
            Whether to plot reciprocal lattice vectors.
        colour : str, default "#FF4D67"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        label : str, optional
            Label for the plot.
        vector_pad : float, default 1.1
            Multiplier for the position of the vectors labels. 1 = position of the vector.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        self.plot_wigner_seitz(
            ax,
            vectors=vectors,
            colour=colour,
            label=label,
            vector_pad=vector_pad,
            reciprocal=True,
            normalize=normalize,
        )

    def plot_wigner_seitz(
        self,
        ax=None,
        vectors=True,
        colour="black",
        label=None,
        vector_pad=1.1,
        reciprocal=False,
        normalize=False,
    ):
        r"""
        Plot Wigner-Seitz unit cell.

        ax : axis, optional
            3D matplotlib axis for the plot.
        vectors : bool, default True
            Whether to plot lattice vectors.
        colour : str, default "black" or "#FF4D67"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
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

        self._artists[artist_group] = []

        if ax is None:
            self.prepare_figure()
            ax = self._ax

        if reciprocal:
            v1, v2, v3 = self.b1, self.b2, self.b3
            v_literal = "b"
            if colour is None:
                colour = "#FF4D67"
        else:
            v1, v2, v3 = self.a1, self.a2, self.a3
            v_literal = "a"
            if colour is None:
                colour = "black"

        if normalize:
            factor = volume(v1, v2, v3) ** (1 / 3.0)
            v1 /= factor
            v2 /= factor
            v3 /= factor

        if label is not None:
            self._artists[artist_group].append(
                ax.scatter(0, 0, 0, color=colour, label=label)
            )
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            self._artists[artist_group].append(
                ax.text(
                    v1[0] * vector_pad[0],
                    v1[1] * vector_pad[0],
                    v1[2] * vector_pad[0],
                    f"${v_literal}_1$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            self._artists[artist_group].append(
                ax.text(
                    v2[0] * vector_pad[1],
                    v2[1] * vector_pad[1],
                    v2[2] * vector_pad[1],
                    f"${v_literal}_2$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            self._artists[artist_group].append(
                ax.text(
                    v3[0] * vector_pad[2],
                    v3[1] * vector_pad[2],
                    v3[2] * vector_pad[2],
                    f"${v_literal}_3$",
                    fontsize=20,
                    color=colour,
                    ha="center",
                    va="center",
                )
            )
            for v in [v1, v2, v3]:
                # Try beautiful arrows
                try:
                    self._artists[artist_group].append(
                        ax.add_artist(
                            Arrow3D(
                                ax,
                                [0, v[0]],
                                [0, v[1]],
                                [0, v[2]],
                                mutation_scale=20,
                                arrowstyle="-|>",
                                color=colour,
                                lw=2,
                                alpha=0.8,
                            )
                        )
                    )
                # Go to default
                except:
                    self._artists[artist_group].append(
                        ax.quiver(
                            0,
                            0,
                            0,
                            *tuple(v),
                            arrow_length_ratio=0.2,
                            color=colour,
                            alpha=0.5,
                        )
                    )
                # Ghost point to account for the plot range
                self._artists[artist_group].append(ax.scatter(*tuple(v), s=0))

        edges, vertices = self.voronoi_cell(reciprocal=reciprocal, normalize=normalize)
        for p1, p2 in edges:
            self._artists[artist_group].append(
                ax.plot(
                    [p1[0], p2[0]],
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    color=colour,
                )
            )

    def plot_conventional(self, **kwargs):
        r"""
        Plot conventional unit cell.

        See Also
        --------
        plot_real_space : for the list of parameters
        """

        self.plot_real_space(conventional=True, **kwargs)

    def plot_primitive(self, **kwargs):
        r"""
        Plot primitive unit cell.

        See Also
        --------
        plot_real_space : for the list of parameters
        """

        self.plot_real_space(**kwargs)

    def plot_kpath(self, ax=None, colour="black", label=None, normalize=False):
        r"""
        Plot k path in the reciprocal space.

        ax : axes
            Axes for the plot. 3D.
        colour : str, default "black"
            Colour for the plot. Any format supported by matplotlib. See |matplotlibColor|_.
        label : str, optional
            Label for the plot.
        normalize : bool, default False
            Whether to normalize corresponding vectors to have the volume equal to one.
        """

        artist_group = "kpath"

        self._artists[artist_group] = []

        if ax is None:
            self.prepare_figure()
            ax = self._ax

        cell = self.reciprocal_cell

        if normalize:
            cell /= volume(cell) ** (1 / 3.0)

        for point in self.points:
            self._artists[artist_group].append(
                ax.scatter(
                    *tuple(self.points[point] @ cell),
                    s=36,
                    color=colour,
                )
            )
            self._artists[artist_group].append(
                ax.text(
                    *tuple(
                        self.points[point] @ cell
                        + 0.025 * cell[0]
                        + +0.025 * cell[1]
                        + 0.025 * cell[2]
                    ),
                    self._PLOT_NAMES[point],
                    fontsize=20,
                    color=colour,
                )
            )
        if label is not None:
            self._artists[artist_group].append(
                ax.scatter(
                    0,
                    0,
                    0,
                    s=36,
                    color=colour,
                    label=label,
                )
            )

        for subpath in self.path:
            for i in range(len(subpath) - 1):
                self._artists[artist_group].append(
                    ax.plot(
                        *tuple(
                            np.concatenate(
                                (
                                    self.points[subpath[i]] @ cell,
                                    self.points[subpath[i + 1]] @ cell,
                                )
                            )
                            .reshape(2, 3)
                            .T
                        ),
                        color=colour,
                        alpha=0.5,
                        linewidth=3,
                    )
                )

    def plot_brillouin_kpath(
        self, zone_colour="#FF4D67", path_colour="black", **kwargs
    ):
        r"""
        Plot brillouin zone and kpath.

        Parameters
        ----------
        zone_colour : str, default "#FF4D67"
            Colour for the brillouin zone. Any format supported by matplotlib. See |matplotlibColor|_.
        zone_colour : str, default "black"
            Colour for the k path. Any format supported by matplotlib. See |matplotlibColor|_.

        See Also
        --------
        plot_brillouin : plot brillouin zone
        plot_kpath : plot kpath
        """

        self.plot_brillouin(colour=zone_colour, **kwargs)
        self.plot_kpath(colour=path_colour, **kwargs)


def get_niggli(
    a=1,
    b=1,
    c=1,
    alpha=90,
    beta=90,
    gamma=90,
    eps_rel=1e-5,
    verbose=False,
    return_cell=False,
):
    r"""
    Computes Niggli matrix form.

    Parameters
    ----------
    a : float, default 1
        Length of the :math:`a_1` vector.
    b : float, default 1
        Length of the :math:`a_2` vector.
    c : float, default 1
        Length of the :math:`a_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [2]_.
    verbose : bool, default False
        Whether to print the steps of an algorithm.
    return_cell : bool, default False
        Whether to return cell parameters instead of Niggli matrix form.

    Returns
    -------
    result : (3,2) :numpy:`ndarray`
        Niggli matrix form as defined in [1]_:

        .. math::

            \begin{pmatrix}
                A & B & C \\
                \xi/2 & \eta/2 & \zeta/2
            \end{pmatrix}
        
        If return_cell == True, then return Niggli cell: (a, b, c, alpha, beta, gamma).


    References
    ----------
    .. [1] Křivý, I. and Gruber, B., 1976.
        A unified algorithm for determining the reduced (Niggli) cell.
        Acta Crystallographica Section A: Crystal Physics, Diffraction,
        Theoretical and General Crystallography,
        32(2), pp.297-298.
    .. [2] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    Examples
    --------
    Example from [1]_
    (parameters are reproducing :math:`A=9`, :math:`B=27`, :math:`C=4`, 
    :math:`\xi` = -5, :math:`\eta` = -4, :math:`\zeta = -22`):

    .. doctest::

        >>> import rad_tools as rad
        >>> from rad_tools.routines import _todegrees
        >>> from math import acos, sqrt
        >>> a = 3
        >>> b = sqrt(27)
        >>> c = 2
        >>> print(f"{a} {b:.3f} {c}")
        3 5.196 2
        >>> alpha = acos(-5 / 2 / b / c) * _todegrees
        >>> beta = acos(-4 / 2 / a / c) * _todegrees
        >>> gamma = acos(-22 / 2 / a / b) * _todegrees
        >>> print(f"{alpha:.2f} {beta:.2f} {gamma:.2f}")
        103.92 109.47 134.88
        >>> niggli_matrix_form = rad.get_niggli(a, b, c, alpha, beta, gamma, verbose=True)
                       A         B         C        xi        eta      zeta   
        start:       9.00000  27.00000   4.00000  -5.00000  -4.00000 -22.00000
        2 appl. to   9.00000  27.00000   4.00000  -5.00000  -4.00000 -22.00000
        1 appl. to   9.00000   4.00000  27.00000  -5.00000 -22.00000  -4.00000
        4 appl. to   4.00000   9.00000  27.00000 -22.00000  -5.00000  -4.00000
        5 appl. to   4.00000   9.00000  27.00000 -22.00000  -5.00000  -4.00000
        4 appl. to   4.00000   9.00000  14.00000  -4.00000  -9.00000  -4.00000
        6 appl. to   4.00000   9.00000  14.00000  -4.00000  -9.00000  -4.00000
        4 appl. to   4.00000   9.00000   9.00000  -8.00000  -1.00000  -4.00000
        7 appl. to   4.00000   9.00000   9.00000  -8.00000  -1.00000  -4.00000
        3 appl. to   4.00000   9.00000   9.00000  -9.00000  -1.00000   4.00000
        5 appl. to   4.00000   9.00000   9.00000   9.00000   1.00000   4.00000
        3 appl. to   4.00000   9.00000   9.00000  -9.00000  -3.00000   4.00000
        result:      4.00000   9.00000   9.00000   9.00000   3.00000   4.00000
        >>> niggli_matrix_form
        array([[4. , 9. , 9. ],
               [4.5, 1.5, 2. ]])

    """

    eps = eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)
    n = abs(floor(log10(abs(eps))))

    # 0
    A = a**2
    B = b**2
    C = c**2
    xi = 2 * b * c * cos(alpha * _toradians)
    eta = 2 * a * c * cos(beta * _toradians)
    zeta = 2 * a * b * cos(gamma * _toradians)
    N = (
        max(
            len(str(A).split(".")[0]),
            len(str(B).split(".")[0]),
            len(str(C).split(".")[0]),
            len(str(xi).split(".")[0]),
            len(str(eta).split(".")[0]),
            len(str(zeta).split(".")[0]),
        )
        + 1
        + n
    )

    def compare(x, condition, y):
        if condition == "<":
            return x < y - eps
        if condition == ">":
            return y < x - eps
        if condition == "<=":
            return not y < x - eps
        if condition == ">=":
            return not x < y - eps
        if condition == "==":
            return not (x < y - eps or y < x - eps)

    def summary_line():
        return (
            f"{A:{N}.{n}f} {B:{N}.{n}f} {C:{N}.{n}f} "
            + f"{xi:{N}.{n}f} {eta:{N}.{n}f} {zeta:{N}.{n}f}"
        )

    if verbose:
        print(
            f"           {'A':^{N}} {'B':^{N}} {'C':^{N}} "
            + f"{'xi':^{N}} {'eta':^{N}} {'zeta':^{N}}"
        )
        cprint(
            f"start:     {summary_line()}",
            color="yellow",
        )
        phrase = "appl. to"
    while True:
        # 1
        if compare(A, ">", B) or (
            compare(A, "==", B) and compare(abs(xi), ">", abs(eta))
        ):
            if verbose:
                print(f"1 {phrase} {summary_line()}")
            A, xi, B, eta = B, eta, A, xi
        # 2
        if compare(B, ">", C) or (
            compare(B, "==", C) and compare(abs(eta), ">", abs(zeta))
        ):
            if verbose:
                print(f"2 {phrase} {summary_line()}")
            B, eta, C, zeta = C, zeta, B, eta
            continue
        # 3
        if compare(xi * eta * zeta, ">", 0):
            if verbose:
                print(f"3 {phrase} {summary_line()}")
            xi, eta, zeta = abs(xi), abs(eta), abs(zeta)
        # 4
        if compare(xi * eta * zeta, "<=", 0):
            if verbose:
                print(f"4 {phrase} {summary_line()}")
            xi, eta, zeta = -abs(xi), -abs(eta), -abs(zeta)
        # 5
        if (
            compare(abs(xi), ">", B)
            or (compare(xi, "==", B) and compare(2 * eta, "<", zeta))
            or (compare(xi, "==", -B) and compare(zeta, "<", 0))
        ):
            if verbose:
                print(f"5 {phrase} {summary_line()}")
            C = B + C - xi * np.sign(xi)
            eta = eta - zeta * np.sign(xi)
            xi = xi - 2 * B * np.sign(xi)
            continue
        # 6
        if (
            compare(abs(eta), ">", A)
            or (compare(eta, "==", A) and compare(2 * xi, "<", zeta))
            or (compare(eta, "==", -A) and compare(zeta, "<", 0))
        ):
            if verbose:
                print(f"6 {phrase} {summary_line()}")
            C = A + C - eta * np.sign(eta)
            xi = xi - zeta * np.sign(eta)
            eta = eta - 2 * A * np.sign(eta)
            continue
        # 7
        if (
            compare(abs(zeta), ">", A)
            or (compare(zeta, "==", A) and compare(2 * xi, "<", eta))
            or (compare(zeta, "==", -A) and compare(eta, "<", 0))
        ):
            if verbose:
                print(f"7 {phrase} {summary_line()}")
            B = A + B - zeta * np.sign(zeta)
            xi = xi - eta * np.sign(zeta)
            zeta = zeta - 2 * A * np.sign(zeta)
            continue
        # 8
        if compare(xi + eta + zeta + A + B, "<", 0) or (
            compare(xi + eta + zeta + A + B, "==", 0)
            and compare(2 * (A + eta) + zeta, ">", 0)
        ):
            if verbose:
                print(f"8 {phrase} {summary_line()}")
            C = A + B + C + xi + eta + zeta
            xi = 2 * B + xi + zeta
            eta = 2 * A + eta + zeta
            continue
        break
    if verbose:
        cprint(
            f"result:    {summary_line()}",
            color="green",
        )
    if return_cell:
        a = round(sqrt(A), n)
        b = round(sqrt(B), n)
        c = round(sqrt(C), n)
        alpha = round(acos(xi / 2 / b / c) * _todegrees, n)
        beta = round(acos(eta / 2 / a / c) * _todegrees, n)
        gamma = round(acos(zeta / 2 / a / b) * _todegrees, n)
        return a, b, c, alpha, beta, gamma
    return np.around(np.array([[A, B, C], [xi / 2, eta / 2, zeta / 2]]), decimals=n)


def lepage(
    a=1,
    b=1,
    c=1,
    alpha=90,
    beta=90,
    gamma=90,
    eps_rel=1e-4,
    verbose=False,
    delta_max=None,
):
    r"""
    Le Page algorithm [1_].

    Parameters
    ----------
    a : float, default 1
        Length of the :math:`a_1` vector.
    b : float, default 1
        Length of the :math:`a_2` vector.
    c : float, default 1
        Length of the :math:`a_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.
    eps_rel : float, default 1e-4
        Relative epsilon as defined in [2]_.
    verbose : bool, default False
        Whether to print the steps of an algorithm.

    References
    ----------
    .. [1] Le Page, Y., 1982.
        The derivation of the axes of the conventional unit cell from
        the dimensions of the Buerger-reduced cell.
        Journal of Applied Crystallography, 15(3), pp.255-259.
    """

    eps = eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)
    if delta_max is None:
        delta_max = eps
    target_axes = {
        "CUB": np.concatenate(
            (
                np.zeros(24),
                abs(cos(60 * _toradians)) * np.ones(24),
                abs(cos(45 * _toradians)) * np.ones(24),
                np.ones(9),
            )
        ),
        "HEX": np.concatenate(
            (
                np.zeros(18),
                abs(cos(60 * _toradians)) * np.ones(12),
                abs(cos(30 * _toradians)) * np.ones(12),
                np.ones(7),
            )
        ),
        "TET": np.concatenate(
            (np.zeros(12), abs(cos(45 * _toradians)) * np.ones(8), np.ones(5))
        ),
        "RHL": np.concatenate((abs(cos(60 * _toradians)) * np.ones(6), np.ones(3))),
        "ORC": np.concatenate((np.zeros(6), np.ones(3))),
        "MCL": np.concatenate((np.zeros(6), np.ones(3))),
    }

    conventional_axis = {
        "CUB": np.concatenate(
            (
                np.zeros(4),
                abs(cos(45 * _toradians)) * np.ones(4),
                np.ones(1),
            )
        ),
        "TET": {
            "z": np.concatenate((np.zeros(4), np.ones(1))),
            "xy": np.concatenate(
                (np.zeros(2), abs(cos(45 * _toradians)) * np.ones(2), np.ones(1))
            ),
        },
    }

    limit = 1.5

    # Niggli reduction
    a, b, c, alpha, beta, gamma = get_niggli(
        a, b, c, alpha, beta, gamma, return_cell=True
    )
    lattice = Lattice(a, b, c, alpha, beta, gamma)

    # find all axes
    miller_indices = (np.indices((5, 5, 5)) - 2).transpose((1, 2, 3, 0)).reshape(125, 3)
    axes = []
    for d_i in miller_indices:
        for r_i in miller_indices:
            if abs(d_i @ r_i) == 2:
                t = d_i @ lattice.cell
                tau = r_i @ lattice.reciprocal_cell
                delta = (
                    np.arctan(np.linalg.norm(np.cross(t, tau)) / abs(t @ tau))
                    * _todegrees
                )
                if delta < limit:
                    axes.append([d_i, t / np.linalg.norm(t), abs(d_i @ r_i), delta])

    # sort and filter
    axes.sort(key=lambda x: x[-1])
    keep_index = np.ones(len(axes))
    for i in range(len(axes)):
        if keep_index[i]:
            for j in range(i + 1, len(axes)):
                if (
                    (axes[i][0] == axes[j][0]).all()
                    or (axes[i][0] == -axes[j][0]).all()
                    or (axes[i][0] == 2 * axes[j][0]).all()
                    or (axes[i][0] == 0.5 * axes[j][0]).all()
                ):
                    keep_index[i] = 0
                    break
    new_axes = []
    for i in range(len(axes)):
        if keep_index[i]:
            if set(axes[i][0]) == {0, 2}:
                axes[i][0] = axes[i][0] / 2
            new_axes.append(axes[i])
    axes = new_axes

    n = len(axes)
    angle_set = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            angle_set[i][j] = np.array(axes[i][1]) @ np.array(axes[j][1])

    delta = None
    while delta is None or delta >= delta_max:
        found_candidate = False
        n = len(axes)
        result = None
        try:
            delta = max(axes, key=lambda x: x[-1])[-1]
        except ValueError:
            delta = 0
        print(f"delta     = {delta:11.8f}")
        print(f"delta_max = {delta_max:11.8f}")

        # CUB
        if (
            n**2 == target_axes["CUB"].shape[0]
            and (
                np.abs(np.sort(np.abs(angle_set.flatten())) - target_axes["CUB"]) < eps
            ).all()
        ):
            xyz = []
            for i in range(n):
                if (
                    np.abs(np.sort(np.abs(angle_set[i])) - conventional_axis["CUB"])
                    < eps
                ).all():
                    xyz.append(axes[i])
            det = np.abs(
                np.linalg.det(
                    [
                        xyz[0][0],
                        xyz[1][0],
                        xyz[2][0],
                    ]
                )
            )
            if det == 1:
                result = "CUB"
                found_candidate = True
            elif det == 4:
                result = "FCC"
                found_candidate = True
            elif det == 2:
                result = "BCC"
                found_candidate = True

        # HEX
        if not found_candidate and (
            n**2 == target_axes["HEX"].shape[0]
            and (
                np.abs(np.sort(np.abs(angle_set.flatten())) - target_axes["HEX"]) < eps
            ).all()
        ):
            result = "HEX"
            found_candidate = True

        # TET
        if not found_candidate and (
            n**2 == target_axes["TET"].shape[0]
            and (
                np.abs(np.sort(np.abs(angle_set.flatten())) - target_axes["TET"]) < eps
            ).all()
        ):
            x, y, z = None, None, None
            x_alt, y_alt = None, None
            for i in range(n):
                if (
                    np.abs(
                        np.sort(np.abs(angle_set[i])) - conventional_axis["TET"]["z"]
                    )
                    < eps
                ).all():
                    z = axes[i]
                if (
                    np.abs(
                        np.sort(np.abs(angle_set[i])) - conventional_axis["TET"]["xy"]
                    )
                    < eps
                ).all():
                    if x is None:
                        x = axes[i]
                    elif y is None and abs(axes[i][1] @ x[1]) < eps:
                        y = axes[i]
                    elif x_alt is None:
                        x_alt = axes[i]
                    else:
                        y_alt = axes[i]

            if np.linalg.norm(x_alt[0] @ lattice.cell) < np.linalg.norm(
                x[0] @ lattice.cell
            ):
                x, y = x_alt, y_alt

            det = np.abs(
                np.linalg.det(
                    [
                        x[0],
                        y[0],
                        z[0],
                    ]
                )
            )
            if det == 1:
                result = "TET"
                found_candidate = True
            elif det == 2:
                result = "BCT"
                found_candidate = True

        # RHL
        if not found_candidate and (
            n**2 == target_axes["RHL"].shape[0]
            and (
                np.abs(np.sort(np.abs(angle_set.flatten())) - target_axes["RHL"]) < eps
            ).all()
        ):
            result = "RHL"
            found_candidate = True

        # ORC
        if not found_candidate and (
            n**2 == target_axes["ORC"].shape[0]
            and (
                np.abs(np.sort(np.abs(angle_set.flatten())) - target_axes["ORC"]) < eps
            ).all()
        ):
            C = np.array(
                [
                    axes[0][0],
                    axes[1][0],
                    axes[2][0],
                ],
                dtype=float,
            ).T
            det = np.abs(np.linalg.det(C))
            if det == 1:
                result = "ORC"
                found_candidate = True
            if det == 4:
                result = "ORCF"
                found_candidate = True
            if det == 2:
                v1, v2, v3, v4 = (
                    C @ [0, 1, 1],
                    C @ [1, 0, 1],
                    C @ [1, 1, 0],
                    C @ [1, 1, 1],
                )

                def gcd(p, q):
                    while q != 0:
                        p, q = q, p % q
                    return p

                if (
                    gcd(abs(v4[0]), abs(v4[1])) > 1
                    and gcd(abs(v4[0]), abs(v4[2])) > 1
                    and gcd(abs(v4[1]), abs(v4[2])) > 1
                ):
                    result = "ORCI"
                    found_candidate = True
                else:
                    result = "ORCC"
                    found_candidate = True

        # MCL
        if not found_candidate and (n == 1):
            v = axes[0][0] @ lattice.cell
            a, b, c = lattice.cell

            ax = []
            test_ax = []

            if abs(a @ v) < eps:
                ax.append(np.array([1, 0, 0]))
            else:
                test_ax.append(np.array([1, 0, 0]))
            if abs(b @ v) < eps:
                ax.append(np.array([0, 1, 0]))
            else:
                test_ax.append(np.array([0, 1, 0]))
            if abs(c @ v) < eps:
                ax.append(np.array([0, 0, 1]))
            else:
                test_ax.append(np.array([0, 0, 1]))

            indices = [[1, 0], [0, 1], [1, 1], [1, -1]]
            if len(ax) == 2:
                a, b = ax
                ax = [a, b, a + b, a - b]
                ax.sort(key=lambda x: np.linalg.norm(x @ lattice.cell))
                a = ax[0]
                b = ax[1]
                c = axes[0][0]
            elif len(test_ax) == 2:
                tmp = ax
                a, b = test_ax
                ax = [a, b, a + b, a - b]
                new_ax = []
                for i in ax:
                    if abs((i @ lattice.cell) @ v) < eps:
                        new_ax.append(i)

                a, b = tmp[0], new_ax[0]
                ax = [a, b, a + b, a - b]
                ax.sort(key=lambda x: np.linalg.norm(x @ lattice.cell))
                a = ax[0]
                b = ax[1]
                c = axes[0][0]

            C = np.array(
                [
                    a,
                    b,
                    c,
                ],
                dtype=float,
            ).T
            det = np.abs(np.linalg.det(C))
            if det == 1:
                result = "MCL"
                found_candidate = True
            if det == 2:
                result = "MCLC"
                found_candidate = True

        # TRI
        if not found_candidate:
            result = "TRI"

        if len(axes) > 0:
            # remove worst axes
            while len(axes) >= 2 and axes[-1][-1] == axes[-2][-1]:
                axes = axes[:-1]
                angle_set = np.delete(angle_set, -1, -1)[:-1]
            axes = axes[:-1]
            angle_set = np.delete(angle_set, -1, -1)[:-1]
        print(f"System: {result}")
        print("=====END WHILE=====")

    return result


if __name__ == "__main__":
    print(
        lepage(
            4,
            4.472,
            4.583,
            79.030,
            64.130,
            64.150,
            verbose=True,
            eps_rel=0.001,
            delta_max=0.006,
        )
    )

    from rad_tools.crystal.bravais_lattice import lattice_example

    delta_max = [
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        1e-5,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        1e-5,
        None,
        None,
        None,
        None,
    ]
    for i, name in enumerate(
        [
            "CUB",
            "FCC",
            "BCC",
            "HEX",
            "TET",
            "BCT1",
            "BCT2",
            "RHL1",
            "RHL2",
            "ORC",
            "ORCF1",
            "ORCF2",
            "ORCF3",
            "ORCI",
            "ORCC",
            "MCL",
            "MCLC1",
            "MCLC2",
            "MCLC3",
            "MCLC4",
            "MCLC5",
            "TRI1a",
            "TRI2a",
            "TRI1b",
            "TRI2b",
        ]
    ):
        lattice = lattice_example(name)
        print("\n" + "=" * 80)
        print(
            name,
            lepage(
                lattice.a,
                lattice.b,
                lattice.c,
                lattice.alpha,
                lattice.beta,
                lattice.gamma,
                verbose=True,
                delta_max=delta_max[i],
            ),
        )

    print("\n" + "=" * 80)
