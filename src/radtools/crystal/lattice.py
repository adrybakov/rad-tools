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

r"""
General 3D lattice.
"""

import numpy as np
from scipy.spatial import Voronoi

import radtools.crystal.cell as Cell
from radtools.crystal.bravais_lattice.hs_points import (
    BCC_hs_points,
    BCT_hs_points,
    CUB_hs_points,
    FCC_hs_points,
    HEX_hs_points,
    MCL_hs_points,
    MCLC_hs_points,
    ORC_hs_points,
    ORCC_hs_points,
    ORCF_hs_points,
    ORCI_hs_points,
    RHL_hs_points,
    TET_hs_points,
    TRI_hs_points,
)
from radtools.crystal.bravais_lattice.standardize import standardize_cell
from radtools.crystal.bravais_lattice.variations import (
    BCT_variation,
    MCLC_variation,
    ORCF_variation,
    RHL_variation,
    TRI_variation,
)
from radtools.crystal.constants import (
    BRAVAIS_LATTICE_NAMES,
    DEFAULT_K_PATHS,
    HS_PLOT_NAMES,
    PEARSON_SYMBOLS,
    REL_TOL,
    TRANSFORM_TO_CONVENTIONAL,
)
from radtools.crystal.identify import lepage
from radtools.crystal.kpoints import Kpoints
from radtools.geometry import angle, volume

__all__ = ["Lattice"]


class Lattice:
    r"""
    General 3D lattice.

    When created from the cell orientation of the cell is respected,
    however the lattice vectors may be renamed with respect to [1]_.

    Creation may change the angles and the lengths of the cell vectors.
    It preserve the volume, right- or left- handedness, lattice type and variation
    of the cell.

    The lattice vector`s lengths are preserved as a set.

    The angles between the lattice vectors are preserved as a set with possible
    changes of the form: :math:`angle \rightarrow 180 - angle`.

    The returned cell may not be the same as the input one, but it is translationally
    equivalent.

    Lattice can be created in a three alternative ways:

    .. doctest::

        >>> import radtools as rad
        >>> l = rad.Lattice(cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> l = rad.Lattice(a1 = [1,0,0], a2 = [0,1,0], a3 = [0,0,1])
        >>> l = rad.Lattice(a=1, b=1, c=1, alpha=90, beta=90, gamma=90)

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
    standardize : bool, default True
        Whether to standardize the cell.
        The consistence of the predefined k paths is not guaranteed in the cell is not unified.

    Attributes
    ----------
    eps_rel : float, default 1e-4
        Relative error for the :ref:`library_lepage` algorithm.

    References
    ----------
    .. [1] Setyawan, W. and Curtarolo, S., 2010.
        High-throughput electronic band structure calculations: Challenges and tools.
        Computational materials science, 49(2), pp.299-312.
    """

    def __init__(self, *args, standardize=True, **kwargs) -> None:
        self._eps_rel = REL_TOL
        self._cell = None
        self._type = None
        self._kpoints = None
        if "cell" in kwargs:
            cell = kwargs["cell"]
        elif "a1" in kwargs and "a2" in kwargs and "a3" in kwargs:
            cell = np.array([kwargs["a1"], kwargs["a2"], kwargs["a3"]])
        elif (
            "a" in kwargs
            and "b" in kwargs
            and "c" in kwargs
            and "alpha" in kwargs
            and "beta" in kwargs
            and "gamma" in kwargs
        ):
            cell = Cell.from_params(
                kwargs["a"],
                kwargs["b"],
                kwargs["c"],
                kwargs["alpha"],
                kwargs["beta"],
                kwargs["gamma"],
            )
        elif len(args) == 1:
            cell = np.array(args[0])
        elif len(args) == 3:
            cell = np.array(args)
        elif len(args) == 6:
            a, b, c, alpha, beta, gamma = args
            cell = Cell.from_params(a, b, c, alpha, beta, gamma)
        elif len(args) == 0 and len(kwargs) == 0:
            cell = np.eye(3)
        else:
            raise ValueError(
                "Unable to identify input parameters. "
                + "Supported: cell ((3,3) array_like), "
                + "or a1, a2, a3 (each is an (3,) array_like), "
                + "or a, b, c, alpha, beta, gamma (floats)."
            )

        self.fig = None
        self.ax = None

        self._set_cell(cell, standardize=standardize)

    # Primitive cell parameters
    @property
    def cell(self):
        r"""
        Unit cell of the lattice.

        Notes
        -----
        In order to rotate the cell with an arbitrary rotation matrix :math:`R` use the syntax:

        .. code-block:: python

            rotated_cell = cell @ R.T

        Transpose is required, since the vectors are stored as rows.

        Returns
        -------
        cell : (3, 3) :numpy:`ndarray`
            Unit cell, rows are vectors, columns are coordinates.
        """
        if self._cell is None:
            raise AttributeError(f"Cell is not defined for lattice {self}")
        return np.array(self._cell)

    # For the child`s overriding
    def _set_cell(self, new_cell, standardize=True):
        try:
            new_cell = np.array(new_cell)
        except:
            raise ValueError(f"New cell is not array_like: {new_cell}")
        if new_cell.shape != (3, 3):
            raise ValueError(f"New cell is not 3 x 3 matrix.")
        self._cell = new_cell
        # Reset type
        self._type = None
        # Standardize cell
        if standardize:
            self._cell = standardize_cell(
                self._cell, self.type(), rtol=self.eps_rel, atol=self.eps
            )

    @cell.setter
    def cell(self, new_cell):
        self._set_cell(new_cell)

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

    # Conventional cell parameters
    @property
    def conv_cell(self):
        r"""
        Conventional cell.

        Returns
        -------
        conv_cell : (3, 3) :numpy:`ndarray`
            Conventional cell, rows are vectors, columns are coordinates.
        """

        return TRANSFORM_TO_CONVENTIONAL[self.type()] @ self.cell

    @property
    def conv_a1(self):
        r"""
        First vector of the conventional cell.

        Returns
        -------
        conv_a1 : (3,) :numpy:`ndarray`
            First vector of the conventional cell.
        """

        return self.conv_cell[0]

    @property
    def conv_a2(self):
        r"""
        Second vector of the conventional cell.

        Returns
        -------
        conv_a2 : (3,) :numpy:`ndarray`
            Second vector of the conventional cell.
        """

        return self.conv_cell[1]

    @property
    def conv_a3(self):
        r"""
        Third vector of the conventional cell.

        Returns
        -------
        conv_a3 : (3,) :numpy:`ndarray`
            Third vector of the conventional cell.
        """

        return self.conv_cell[2]

    @property
    def conv_a(self):
        r"""
        Length of the first vector of the conventional cell.

        Returns
        -------
        conv_a : float
            Length of the first vector of the conventional cell.
        """

        return np.linalg.norm(self.conv_a1)

    @property
    def conv_b(self):
        r"""
        Length of the second vector of the conventional cell.

        Returns
        -------
        conv_b : float
            Length of the second vector of the conventional cell.
        """

        return np.linalg.norm(self.conv_a2)

    @property
    def conv_c(self):
        r"""
        Length of the third vector of the conventional cell.

        Returns
        -------
        conv_c : float
            Length of the third vector of the conventional cell.
        """

        return np.linalg.norm(self.conv_a3)

    @property
    def conv_alpha(self):
        r"""
        Angle between second and third conventional lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.conv_a2, self.conv_a3)

    @property
    def conv_beta(self):
        r"""
        Angle between first and third conventional lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.conv_a1, self.conv_a3)

    @property
    def conv_gamma(self):
        r"""
        Angle between first and second conventional lattice vector.

        Returns
        -------
        angle : float
            In degrees.
        """

        return angle(self.conv_a1, self.conv_a2)

    @property
    def conv_unit_cell_volume(self):
        r"""
        Volume of the conventional unit cell.

        Returns
        -------
        volume : float
            Unit cell volume.
        """

        return volume(self.conv_a1, self.conv_a2, self.conv_a3)

    @property
    def conv_parameters(self):
        r"""
        Return conventional cell parameters.

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
        return (
            self.conv_a,
            self.conv_b,
            self.conv_c,
            self.conv_alpha,
            self.conv_beta,
            self.conv_gamma,
        )

    # Reciprocal parameters
    @property
    def reciprocal_cell(self):
        r"""
        Reciprocal cell. Always primitive.

        Returns
        -------
        reciprocal_cell : (3, 3) :numpy:`ndarray`
            Reciprocal cell, rows are vectors, columns are coordinates.
        """

        return Cell.reciprocal(self.cell)

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
    def reciprocal_parameters(self):
        r"""
        Return reciprocal cell parameters.

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
        return self.k_a, self.k_b, self.k_c, self.k_alpha, self.k_beta, self.k_gamma

    # Lattice type routines and properties
    @property
    def eps(self):
        r"""
        Epsilon parameter.

        Derived from :py:attr:`.eps_rel` as
        .. math::

            \epsilon = \epsilon_{rel}\cdot V^{\frac{1}{3}}
        """

        return self.eps_rel * abs(self.unit_cell_volume) ** (1 / 3.0)

    @property
    def eps_rel(self):
        r"""
        Relative epsilon

        Returns
        -------
        eps_rel : float
        """

        return self._eps_rel

    @eps_rel.setter
    def eps_rel(self, new_value):
        try:
            self._eps_rel = float(new_value)
        except ValueError:
            raise ValueError("Not a float")
        self._type = None

    def type(self, eps_rel=None):
        r"""
        Identify the lattice type.

        Parameters
        ----------
        eps_rel : float, optional
            Relative error for the :ref:`library_lepage` algorithm.

        Returns
        -------
        lattice_type : str
            Bravais lattice type.

        See Also
        --------
        lepage : Algoritm for the lattice type identification
        variation : Variation of the lattice, if any.
        """

        if self._type is None or eps_rel is not None:
            if eps_rel is None:
                eps_rel = self.eps_rel

            lattice_type = lepage(
                self.a,
                self.b,
                self.c,
                self.alpha,
                self.beta,
                self.gamma,
                eps_rel=eps_rel,
            )

            self._type = lattice_type
        return self._type

    @property
    def variation(self):
        r"""
        Variation of the lattice, if any.

        For the Bravais lattice with only one variation the :py:meth:`.Lattice.type` is returned.

        Returns
        -------
        variation : str
            Variation of the lattice.

        Examples
        --------

        .. doctest::

            >>> import radtools as rad
            >>> l = rad.lattice_example("CUB")
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
        :py:meth:`.Lattice.type`
        """
        lattice_type = self.type()

        if lattice_type == "BCT":
            result = BCT_variation(self.conv_a, self.conv_c)
        elif lattice_type == "ORCF":
            result = ORCF_variation(self.conv_a, self.conv_b, self.conv_c, self.eps)
        elif lattice_type == "RHL":
            result = RHL_variation(self.conv_alpha, self.eps)
        elif lattice_type == "MCLC":
            result = MCLC_variation(
                self.conv_a,
                self.conv_b,
                self.conv_c,
                self.conv_alpha,
                self.k_gamma,
                self.eps,
            )
        elif lattice_type == "TRI":
            result = TRI_variation(self.k_alpha, self.k_beta, self.k_gamma, self.eps)
        else:
            result = lattice_type

        return result

    @property
    def name(self):
        r"""
        Human-readable name of the Bravais lattice type.

        Returns
        -------
        name : str
            Name of the Bravais lattice type.
        """

        return BRAVAIS_LATTICE_NAMES[self.type()]

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

        return PEARSON_SYMBOLS[self.type()]

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
        Computes Voronoi edges around (0,0,0) point.

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

    @property
    def kpoints(self) -> Kpoints:
        r"""
        Instance of :py:class:`.Kpoints` with the high symmetry points and path.

        Notes
        -----
        When a new instance of the :py:class:`.Kpoints` is assigned to the lattice,
        reciprocal vectors of the lattice are not updated. Reciprocal vectors of new kpoints
        are not updated as well.

        Returns
        -------
        kpoints : :py:class:`.Kpoints`
            Instance of the :py:class:`.Kpoints` class.

        See Also
        --------
        Kpoints : Class for the high symmetry points and path.
        """

        if self._kpoints is None:
            self._kpoints = Kpoints(self.b1, self.b2, self.b3)

            if self.type() == "CUB":
                hs_points = CUB_hs_points()
            elif self.type() == "FCC":
                hs_points = FCC_hs_points()
            elif self.type() == "BCC":
                hs_points = BCC_hs_points()
            elif self.type() == "TET":
                hs_points = TET_hs_points()
            elif self.type() == "BCT":
                hs_points = BCT_hs_points(self.variation, self.conv_a, self.conv_c)
            elif self.type() == "ORC":
                hs_points = ORC_hs_points()
            elif self.type() == "ORCF":
                hs_points = ORCF_hs_points(
                    self.variation, self.conv_a, self.conv_b, self.conv_c
                )
            elif self.type() == "ORCI":
                hs_points = ORCI_hs_points(self.conv_a, self.conv_b, self.conv_c)
            elif self.type() == "ORCC":
                hs_points = ORCC_hs_points(self.conv_a, self.conv_b)
            elif self.type() == "HEX":
                hs_points = HEX_hs_points()
            elif self.type() == "RHL":
                hs_points = RHL_hs_points(self.variation, self.conv_alpha)
            elif self.type() == "MCL":
                hs_points = MCL_hs_points(self.conv_b, self.conv_c, self.conv_alpha)
            elif self.type() == "MCLC":
                hs_points = MCLC_hs_points(
                    self.variation,
                    self.conv_a,
                    self.conv_b,
                    self.conv_c,
                    self.conv_alpha,
                )
            elif self.type() == "TRI":
                hs_points = TRI_hs_points(self.variation)

            for point in hs_points:
                if point == "S" and self.type() == "BCT":
                    self._kpoints.add_hs_point(
                        point, hs_points[point], label="$\\Sigma$"
                    )
                elif point == "S1" and self.type() == "BCT":
                    self._kpoints.add_hs_point(
                        point, hs_points[point], label="$\\Sigma_1$"
                    )
                else:
                    self._kpoints.add_hs_point(
                        point, hs_points[point], label=HS_PLOT_NAMES[point]
                    )

            self._kpoints.path = DEFAULT_K_PATHS[self.variation]

        return self._kpoints

    @kpoints.setter
    def kpoints(self, new_kpoints: Kpoints):
        if not isinstance(new_kpoints, Kpoints):
            raise ValueError(
                f"New kpoints should be an instance of the Kpoints class. "
                + f"Got {type(new_kpoints)} instead."
            )
        self._kpoints = new_kpoints
