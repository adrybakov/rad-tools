r"""
General 3D lattice.
"""

from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d
from scipy.spatial import Voronoi

from radtools.crystal.identify import lepage
from radtools.routines import angle, cell_from_param, reciprocal_cell, volume

__all__ = ["Lattice"]


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


def get_3D_axes(background=True, focal_length=0.2):
    r"""
    Prepare style of the figure for the plot.

    Parameters
    ----------
    background : bool, default True
        Whether to keep the axis on the plot.
    focal_length : float, default 0.2
        See: |matplotlibFocalLength|_

    Returns
    -------
    fig
    ax
    """

    fig = plt.figure(figsize=(6, 6))
    rcParams["axes.linewidth"] = 0
    rcParams["xtick.color"] = "#B3B3B3"
    ax = fig.add_subplot(projection="3d")
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
    return fig, ax


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

        >>> import radtools as rad
        >>> l = rad.Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> l = rad.Lattice([1,0,0], [0,1,0], [0,0,1])
        >>> l = rad.Lattice(1, 1, 1, 90, 90, 90)

    Parameters
    ----------
    cell : (3, 3) |array_like|_
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
    kpoints : dist
        Dictionary of the high symmetry points.
        Coordinates are given in relative coordinates.

        .. code-block:: python

            kpoints = {"Name" : [k_x, k_y, k_z], ...}
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
            self.cell = cell_from_param(a, b, c, alpha, beta, gamma)
        else:
            raise ValueError(
                "Unable to identify input parameters. "
                + "Supported: one (3,3) array_like, or three (3,) array_like, or 6 floats."
            )
        self.points = {}
        self._path = None
        self._default_path = None
        self.fig = None
        self.ax = None
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
        r"""
        Unit cell of the lattice.

        Returns
        -------
        cell : (3, 3) :numpy:`ndarray`
            Unit cell, rows are vectors, columns are coordinates.
        """
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

        Returns
        -------
        path : list of str
            List of the high symmetry points.
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

        Returns
        -------
        a1 : (3,) :numpy:`ndarray`
            First lattice vector :math:`\vec{a}_1`.
        """
        return self.cell[0]

    @property
    def a2(self):
        r"""
        Second lattice vector :math:`\vec{a}_2`.

        Returns
        -------
        a2 : (3,) :numpy:`ndarray`
            Second lattice vector :math:`\vec{a}_2`.
        """
        return self.cell[1]

    @property
    def a3(self):
        r"""
        Third lattice vector :math:`\vec{a}_3`.

        Returns
        -------
        a3 : (3,) :numpy:`ndarray`
            Third lattice vector :math:`\vec{a}_3`.
        """
        return self.cell[2]

    @property
    def a(self):
        r"""
        Length of the first lattice vector :math:`\vert\vec{a}_1\vert`.

        Returns
        -------
        a : float
        """

        return np.linalg.norm(self.cell[0])

    @property
    def b(self):
        r"""
        Length of the second lattice vector :math:`\vert\vec{a}_2\vert`.

        Returns
        -------
        b : float
        """

        return np.linalg.norm(self.cell[1])

    @property
    def c(self):
        r"""
        Length of the third lattice vector :math:`\vert\vec{a}_3\vert`.

        Returns
        -------
        c : float
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
    def parameters(self):
        r"""
        Return cell parameters.

        :math:`(a, b, c, \alpha, \beta, \gamma)`

        Returns
        -------
        a : float
        b : float
        c : float
        alpha : float
        beta : float
        gamma : float
        """
        return self.a, self.b, self.c, self.alpha, self.beta, self.gamma

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        Returns
        -------
        volume : float
            Unit cell volume.
        """

        return volume(self.a1, self.a2, self.a3)

    @property
    def reciprocal_cell(self):
        r"""
        Reciprocal cell. Always primitive.

        Returns
        -------
        reciprocal_cell : (3, 3) :numpy:`ndarray`
            Reciprocal cell, rows are vectors, columns are coordinates.
        """

        return reciprocal_cell(self.cell)

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector.

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{a}_2\times\vec{a}_3

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`

        Returns
        -------
        b1 : (3,) :numpy:`ndarray`
            First reciprocal lattice vector :math:`\vec{b}_1`.
        """

        return self.reciprocal_cell[0]

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector.

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{a}_3\times\vec{a}_1

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`

        Returns
        -------
        b2 : (3,) :numpy:`ndarray`
            Second reciprocal lattice vector :math:`\vec{b}_2`.
        """

        return self.reciprocal_cell[1]

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector.

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}_1\times\vec{a}_2

        where :math:`V = \vec{a}_1\cdot(\vec{a}_2\times\vec{a}_3)`

        Returns
        -------
        b3 : (3,) :numpy:`ndarray`
            Third reciprocal lattice vector :math:`\vec{b}_3`.
        """

        return self.reciprocal_cell[2]

    @property
    def k_a(self):
        r"""
        Length of the first reciprocal lattice vector :math:`\vert\vec{b}_1\vert`.

        Returns
        -------
        k_a : float
        """

        return np.linalg.norm(self.b1)

    @property
    def k_b(self):
        r"""
        Length of the second reciprocal lattice vector :math:`\vert\vec{b}_2\vert`.

        Returns
        -------
        k_b : float
        """

        return np.linalg.norm(self.b2)

    @property
    def k_c(self):
        r"""
        Length of the third reciprocal lattice vector :math:`\vert\vec{b}_3\vert`.

        Returns
        -------
        k_c : float
        """

        return np.linalg.norm(self.b3)

    @property
    def k_alpha(self):
        r"""
        Angle between second and third reciprocal lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.b2, self.b3)

    @property
    def k_beta(self):
        r"""
        Angle between first and third reciprocal lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.b1, self.b3)

    @property
    def k_gamma(self):
        r"""
        Angle between first and second reciprocal lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.b1, self.b2)

    @property
    def reciprocal_cell_volume(self):
        r"""
        Volume of the reciprocal cell.

        .. math::

            V = \vec{b}_1\cdot(\vec{b}_2\times\vec{b}_3)

        Returns
        -------
        volume : float
            Volume of the reciprocal cell.
        """

        return volume(self.b1, self.b2, self.b3)

    @property
    def variation(self):
        r"""
        Variation of the lattice, if any.

        For the :py:class:`.Lattice` return "Lattice".

        For the Bravais lattice with only one variation the name of the class is returned.

        For the variation of each Bravais lattice type see
        corresponding ``variation`` attribute.

        Returns
        -------
        variation : str
            Variation of the lattice.

        Examples
        --------

        .. doctest::

            >>> import radtools as rad
            >>> l = rad.lattice_example("cub")
            >>> l.variation
            'CUB'

        .. doctest::

            >>> import radtools as rad
            >>> l = rad.lattice_example("BCT1")
            >>> l.variation
            'BCT1'

        .. doctest::

            >>> import radtools as rad
            >>> l = rad.lattice_example("MCLC4")
            >>> l.variation
            'MCLC4'

        See Also
        --------
        :py:attr:`.CUB.variation`
        :py:attr:`.FCC.variation`
        :py:attr:`.BCC.variation`
        :py:attr:`.TET.variation`
        :py:attr:`.BCT.variation`
        :py:attr:`.ORC.variation`
        :py:attr:`.ORCF.variation`
        :py:attr:`.ORCI.variation`
        :py:attr:`.ORCC.variation`
        :py:attr:`.HEX.variation`
        :py:attr:`.RHL.variation`
        :py:attr:`.MCL.variation`
        :py:attr:`.MCLC.variation`
        :py:attr:`.TRI.variation`
        """
        return self.__class__.__name__

    def identify(self):
        r"""
        Identify the Bravais lattice type.
        """
        return lepage(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

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

        Returns
        -------
        lattice_points : (N, 3) :numpy:`ndarray`
            N lattice points. Each element is a vector :math:`v = (v_x, v_y, v_z)`.
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

        if self.fig is None:
            self.fig, self.ax = get_3D_axes(
                background=background, focal_length=focal_length
            )

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
            * "brillouin-kpath"
            * "wigner-seitz"

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
        plot_wigner_seitz : "wigner-seitz" plot.
        show : Shows the plot.
        savefig : Save the figure in the file.
        """

        if ax is None:
            self.prepare_figure()
            ax = self.ax
        if isinstance(kind, str):
            kinds = [kind]
        else:
            kinds = kind
        try:
            for kind in kinds:
                kind = kind.replace("-", "_")
                getattr(self, f"plot_{kind}")(ax=ax, **kwargs)
        except AttributeError:
            raise ValueError(f"Plot kind '{kind}' does not exist!")
        ax.relim()
        ax.set_aspect("equal")

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
        if ax is None:
            ax = self.ax
        ax.relim(visible_only=True)
        ax.set_aspect("equal")

    def show(self, elev=30, azim=-60):
        r"""
        Show the figure in the interactive matplotlib window.

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
        Directly passed to the |matplotlibLegend|_.
        """

        self.ax.legend(**kwargs)

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
            ax = self.ax

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
            ax = self.ax

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
            ax = self.ax

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
