r"""
Bravais lattices.
"""

from math import cos, floor, log10, pi, sin, sqrt, tan

import numpy as np

from radtools.crystal.identify import lepage
from radtools.crystal.lattice import Lattice
from radtools.routines import param_from_cell, toradians, volume

__all__ = [
    "CUB",
    "FCC",
    "BCC",
    "TET",
    "BCT",
    "ORC",
    "ORCF",
    "ORCI",
    "ORCC",
    "HEX",
    "RHL",
    "MCL",
    "MCLC",
    "TRI",
    "lattice_example",
    "bravais_lattice_from_param",
    "bravais_lattice_from_cell",
]


class NotEnoughParameters(Exception):
    r"""
    Raised if one tries to create a Bravais lattice without enough parameters.

    Gives a summary of required parameters for the Bravais lattice type.

    Parameters
    ----------
    lattice : str
        Lattice type name.
    parameters : dict
        Dictionary of parameters names and values
    """

    def __init__(self, lattice, parameters):
        self.message = (
            f"\nFor the Bravais lattice type '{lattice}' the cell "
            + f"\nor the following lattice parameters are needed:\n    "
        )
        for i, param in enumerate(parameters):
            self.message += f"{param[0]}"
            if i != len(parameters) - 1:
                self.message += ", "
        self.message += "\nGot:\n"
        for i, param in enumerate(parameters):
            self.message += f"    {param[0]} = {param[1]}"
            if i != len(parameters) - 1:
                self.message += ", "

    def __str__(self):
        return self.message


class CellTypeMismatch(Exception):
    r"""
    Raised when one tries to create a Bravais lattice with a cell (or set of parameters)
    which does not match desired Bravais lattice type.

    Parameters
    ----------
    lattice_type : str
        Target lattice type
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.
    a : float
        Length of the :math:`a_1` vector.
    b : float
        Length of the :math:`a_2` vector.
    c : float
        Length of the :math:`a_3` vector.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.
    correct_lattice_type : str, optional
        Correct lattice type

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.
    """

    def __init__(
        self,
        lattice_type,
        eps_rel,
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        correct_lattice_type=None,
    ):
        n = abs(
            floor(
                log10(abs(eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)))
            )
        )
        self.message = (
            f"\nCell is not '{lattice_type}':\n"
            + f"    Relative epsilon = {eps_rel},\n"
            + f"    a = {a:{n+6}.{n+1}f},\n"
            + f"    b = {a:{n+6}.{n+1}f},\n"
            + f"    c = {c:{n+6}.{n+1}f},\n"
            + f"    alpha = {alpha:{n+6}.{n+1}f},\n"
            + f"    beta = {beta:{n+6}.{n+1}f},\n"
            + f"    gamma = {gamma:{n+6}.{n+1}f},\n"
        )
        if correct_lattice_type is None:
            correct_lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
        self.message += (
            f"Lattice type defined from parameters: '{correct_lattice_type}'"
        )

    def __str__(self):
        return self.message


# 1
class CUB(Lattice):
    r"""
    Cubic (CUB, cP)

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, a)

    .. seealso::

        :ref:`lattice-cub` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vectors of the conventional cel.
    cell : (3,3) |array_like|_, optional
        Cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vectors of the conventional cel.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "cP"

    def __init__(self, a: float = None, cell=None, eps_rel=1e-5) -> None:
        if a is None and cell is None:
            raise NotEnoughParameters("CUB", [("a", a)])
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "CUB":
                raise CellTypeMismatch(
                    "CUB", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            super().__init__(cell)
            self.conv_a = (a + b + c) / 3
        else:
            super().__init__([a, 0, 0], [0, a, 0], [0, 0, a])
            self.conv_a = a
        self.conv_cell = self.cell
        self.points = {
            "G": np.array([0, 0, 0]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
        }
        self._default_path = [["G", "X", "M", "G", "R", "X"], ["M", "R"]]


# 2
class FCC(Lattice):
    r"""
    Face-centred cubic (FCC, cF)

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, a)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (0, a/2, a/2)

        \boldsymbol{a}_2 = (a/2, 0, a/2)

        \boldsymbol{a}_3 = (a/2, a/2, 0)

    .. seealso::

        :ref:`lattice-fcc` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional cel.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional cel.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "cF"

    def __init__(self, a: float = None, cell=None, eps_rel=1e-5) -> None:
        if a is None and cell is None:
            raise NotEnoughParameters("FCC", [("a", a)])
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "FCC":
                raise CellTypeMismatch(
                    "FCC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            super().__init__(cell)
            self.conv_cell = [
                [-1.0, 1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, -1.0],
            ] @ self.cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            self.conv_a = (a + b + c) / 3
        else:
            self.conv_a = a
            super().__init__([0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0])
            self.conv_cell = np.diag([self.conv_a, self.conv_a, self.conv_a])

        self.points = {
            "G": np.array([0, 0, 0]),
            "K": np.array([3 / 8, 3 / 8, 3 / 4]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "U": np.array([5 / 8, 1 / 4, 5 / 8]),
            "W": np.array([1 / 2, 1 / 4, 3 / 4]),
            "X": np.array([1 / 2, 0, 1 / 2]),
        }
        self._default_path = [
            ["G", "X", "W", "K", "G", "L", "U", "W", "L", "K"],
            ["U", "X"],
        ]


# 3
class BCC(Lattice):
    r"""
    Body-centered cubic (BCC, cI)

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, a)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (-a/2, a/2, a/2)

        \boldsymbol{a}_2 = (a/2, -a/2, a/2)

        \boldsymbol{a}_3 = (a/2, a/2, -a/2)

    .. seealso::

        :ref:`lattice-bcc` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "cI"

    def __init__(self, a: float = None, cell=None, eps_rel=1e-5) -> None:
        if a is None and cell is None:
            raise NotEnoughParameters("BCC", [("a", a)])
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "BCC":
                raise CellTypeMismatch(
                    "BCC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            super().__init__(cell)
            self.conv_cell = [
                [0, 1.0, 1.0],
                [1.0, 0, 1.0],
                [1.0, 1.0, 0],
            ] @ self.cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            self.conv_a = (a + b + c) / 3
        else:
            self.conv_a = a
            super().__init__(
                [-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]
            )
            self.conv_cell = np.diag([self.conv_a, self.conv_a, self.conv_a])
        self.points = {
            "G": np.array([0, 0, 0]),
            "H": np.array([1 / 2, -1 / 2, 1 / 2]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "N": np.array([0, 0, 1 / 2]),
        }
        self._default_path = [["G", "H", "N", "G", "P", "H"], ["P", "N"]]


# 4
class TET(Lattice):
    r"""
    Tetragonal (TET, tP)

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    .. seealso::

        :ref:`lattice-tet` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "tP"

    def __init__(
        self, a: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or c is None) and cell is None:
            raise NotEnoughParameters("TET", [("a", a), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "TET":
                raise CellTypeMismatch(
                    "TET", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            if not (a < c - eps or c < a - eps):
                cell = [cell[0], cell[2], cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)
            elif not (b < c - eps or b < c - eps):
                cell = [cell[1], cell[2], cell[0]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            super().__init__(cell)
            self.conv_a = (a + b) / 2
            self.conv_c = c
        else:
            lattice_type = lepage(a, a, c, 90, 90, 90, eps_rel=eps_rel)
            if lattice_type != "TET":
                raise CellTypeMismatch(
                    "TET", eps_rel, a, a, c, 90, 90, 90, lattice_type
                )
            self.conv_a = a
            self.conv_c = c
            super().__init__([a, 0, 0], [0, a, 0], [0, 0, c])
        self.conv_cell = self.cell
        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([0, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }
        self._default_path = [
            ["G", "X", "M", "G", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"],
        ]


# 5
class BCT(Lattice):
    r"""
    Body-centred tetragonal (BCT, tI)

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (-a/2, a/2, c/2)

        \boldsymbol{a}_2 = (a/2, -a/2, c/2)

        \boldsymbol{a}_3 = (a/2, a/2, -c/2)

    .. seealso::

        :ref:`lattice-bct` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of conventional lattice.
    c : float, optional
        Length of the lattice vector of conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of conventional lattice.
    conv_c : float
        Length of the lattice vector of conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "tI"

    def __init__(
        self, a: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or c is None) and cell is None:
            raise NotEnoughParameters("BCT", [("a", a), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "BCT":
                raise CellTypeMismatch(
                    "BCT", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            conv_cell = [[0, 1.0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)
            if not (a < c - eps or c < a - eps):
                cell = [cell[0], cell[2], cell[1]]
            elif not (b < c - eps or b < c - eps):
                cell = [cell[1], cell[2], cell[0]]

            self.conv_cell = [[0, 1.0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            super().__init__(cell)
            self.conv_a = (a + b) / 2
            self.conv_c = c
        else:
            self.conv_a = a
            self.conv_c = c
            super().__init__(
                [-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "BCT":
                raise CellTypeMismatch(
                    "BCT", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            self.conv_cell = np.diag([self.conv_a, self.conv_a, self.conv_c])

        self._PLOT_NAMES["S"] = "$\\Sigma$"
        self._PLOT_NAMES["S1"] = "$\\Sigma_1$"

        if self.variation == "BCT1":
            eta = (1 + self.conv_c**2 / self.conv_a**2) / 4
            self.points = {
                "G": np.array([0, 0, 0]),
                "M": np.array([-1 / 2, 1 / 2, 1 / 2]),
                "N": np.array([0, 1 / 2, 0]),
                "P": np.array([1 / 4, 1 / 4, 1 / 4]),
                "X": np.array([0, 0, 1 / 2]),
                "Z": np.array([eta, eta, -eta]),
                "Z1": np.array([-eta, 1 - eta, eta]),
            }

            self._default_path = [
                ["G", "X", "M", "G", "Z", "P", "N", "Z1", "M"],
                ["X", "P"],
            ]

        elif self.variation == "BCT2":
            eta = (1 + self.conv_a**2 / self.conv_c**2) / 4
            zeta = self.conv_a**2 / (2 * self.conv_c**2)
            self.points = {
                "G": np.array([0, 0, 0]),
                "N": np.array([0, 1 / 2, 0]),
                "P": np.array([1 / 4, 1 / 4, 1 / 4]),
                "S": np.array([-eta, eta, eta]),
                "S1": np.array([eta, 1 - eta, -eta]),
                "X": np.array([0, 0, 1 / 2]),
                "Y": np.array([-zeta, zeta, 1 / 2]),
                "Y1": np.array([1 / 2, 1 / 2, -zeta]),
                "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
            }

            self._default_path = [
                [
                    "G",
                    "X",
                    "Y",
                    "S",
                    "G",
                    "Z",
                    "S1",
                    "N",
                    "P",
                    "Y1",
                    "Z",
                ],
                ["X", "P"],
            ]

    @property
    def variation(self):
        r"""
        Two variations of the Lattice.

        :math:`\text{BCT}_1: c < a` and :math:`\text{BCT}_2: c > a`

        Returns
        -------
        variation : str
            Variation of the lattice. "BCT1" or "BCT2".
        """
        if self.conv_a > self.conv_c:
            return "BCT1"
        elif self.conv_a < self.conv_c:
            return "BCT2"


# 6
class ORC(Lattice):
    r"""
    Orthorhombic (ORC, oP)

    :math:`a < b < c`

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    .. seealso::

        :ref:`lattice-orc` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "oP"

    def __init__(
        self, a: float = None, b: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or b is None or c is None) and cell is None:
            raise NotEnoughParameters("ORC", [("a", a), ("b", b), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORC":
                raise CellTypeMismatch(
                    "ORC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            if b < a - eps:
                # minus preserves right-hand order
                cell = [cell[1], cell[0], -cell[2]]
            if c < a - eps:
                cell = [cell[2], cell[0], cell[1]]
            elif c < b - eps:
                # minus preserves right-hand order
                cell = [cell[0], cell[2], -cell[1]]

            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            super().__init__(cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
        else:
            lattice_type = lepage(a, b, c, 90, 90, 90, eps_rel=eps_rel)
            if lattice_type != "ORC":
                raise CellTypeMismatch(
                    "ORC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            a, b, c = tuple(sorted([a, b, c]))
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            super().__init__([a, 0, 0], [0, b, 0], [0, 0, c])
        self.conv_cell = self.cell
        self.points = {
            "G": np.array([0, 0, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "S": np.array([1 / 2, 1 / 2, 0]),
            "T": np.array([0, 1 / 2, 1 / 2]),
            "U": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }
        self._default_path = [
            ["G", "X", "S", "Y", "G", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"],
        ]


# 7
class ORCF(Lattice):
    r"""
    Face-centred orthorhombic (ORCF, oF)

    :math:`a < b < c`

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (0, b/2, c/2)

        \boldsymbol{a}_2 = (a/2, 0, c/2)

        \boldsymbol{a}_3 = (a/2, b/2, 0)

    .. seealso::

        :ref:`lattice-orcf` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.

    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "oF"

    def __init__(
        self, a: float = None, b: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or b is None or c is None) and cell is None:
            raise NotEnoughParameters("ORCF", [("a", a), ("b", b), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCF":
                raise CellTypeMismatch(
                    "ORCF", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            conv_cell = [[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)
            if b < a - eps:
                cell = [cell[1], cell[0], cell[2]]
            if c < a - eps:
                cell = [cell[2], cell[0], cell[1]]
            elif c < b - eps:
                cell = [cell[0], cell[2], cell[1]]

            super().__init__(cell)
            self.conv_cell = [
                [-1.0, 1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, -1.0],
            ] @ self.cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
        else:
            a, b, c = tuple(sorted([a, b, c]))
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            super().__init__([0, b / 2, c / 2], [a / 2, 0, c / 2], [a / 2, b / 2, 0])
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCF":
                raise CellTypeMismatch(
                    "ORCF", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            self.conv_cell = np.diag([self.conv_a, self.conv_b, self.conv_c])
        if self.variation == "ORCF1":
            eta = (
                1
                + self.conv_a**2 / self.conv_b**2
                + self.conv_a**2 / self.conv_c**2
            ) / 4
            zeta = (
                1
                + self.conv_a**2 / self.conv_b**2
                - self.conv_a**2 / self.conv_c**2
            ) / 4
            self.points = {
                "G": np.array([0, 0, 0]),
                "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
                "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "T": np.array([1, 1 / 2, 1 / 2]),
                "X": np.array([0, eta, eta]),
                "X1": np.array([1, 1 - eta, 1 - eta]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self._default_path = [
                ["G", "Y", "T", "Z", "G", "X", "A1", "Y"],
                ["T", "X1"],
                ["X", "A", "Z"],
                ["L", "G"],
            ]
        elif self.variation == "ORCF2":
            eta = (
                1
                + self.conv_a**2 / self.conv_b**2
                - self.conv_a**2 / self.conv_c**2
            ) / 4
            delta = (
                1
                + self.conv_b**2 / self.conv_a**2
                - self.conv_b**2 / self.conv_c**2
            ) / 4
            phi = (
                1
                + self.conv_c**2 / self.conv_b**2
                - self.conv_c**2 / self.conv_a**2
            ) / 4

            self.points = {
                "G": np.array([0, 0, 0]),
                "C": np.array([1 / 2, 1 / 2 - eta, 1 - eta]),
                "C1": np.array([1 / 2, 1 / 2 + eta, eta]),
                "D": np.array([1 / 2 - delta, 1 / 2, 1 - delta]),
                "D1": np.array([1 / 2 + delta, 1 / 2, delta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "H": np.array([1 - phi, 1 / 2 - phi, 1 / 2]),
                "H1": np.array([phi, 1 / 2 + phi, 1 / 2]),
                "X": np.array([0, 1 / 2, 1 / 2]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self._default_path = [
                ["G", "Y", "C", "D", "X", "G", "Z", "D1", "H", "C"],
                ["C1", "Z"],
                ["X", "H1"],
                ["H", "Y"],
                ["L", "G"],
            ]
        elif self.variation == "ORCF3":
            eta = (
                1
                + self.conv_a**2 / self.conv_b**2
                + self.conv_a**2 / self.conv_c**2
            ) / 4
            zeta = (
                1
                + self.conv_a**2 / self.conv_b**2
                - self.conv_a**2 / self.conv_c**2
            ) / 4

            self.points = {
                "G": np.array([0, 0, 0]),
                "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
                "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "T": np.array([1, 1 / 2, 1 / 2]),
                "X": np.array([0, eta, eta]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self._default_path = [
                ["G", "Y", "T", "Z", "G", "X", "A1", "Y"],
                ["X", "A", "Z"],
                ["L", "G"],
            ]

    @property
    def variation(self):
        r"""
        Three variations of the Lattice.

        :math:`\text{ORCF}_1: \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_2: \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_3: \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}`,

        Returns
        -------
        variation : str
            Variation of the lattice. "ORCF1", "ORCF2" or "ORCF3".
        """

        expresion = 1 / self.conv_a**2 - 1 / self.conv_b**2 - 1 / self.conv_c**2
        if np.abs(expresion) < 10 * np.finfo(float).eps:
            return "ORCF3"
        elif expresion > 0:
            return "ORCF1"
        elif expresion < 0:
            return "ORCF2"


# 8
class ORCI(Lattice):
    r"""
    Body-centred orthorhombic (ORCI, oI)

    :math:`a < b < c`

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (-a/2, b/2, c/2)

        \boldsymbol{a}_2 = (a/2, -b/2, c/2)

        \boldsymbol{a}_3 = (a/2, b/2, -c/2)

    .. seealso::

        :ref:`lattice-orci` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "oI"

    def __init__(
        self, a: float = None, b: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or b is None or c is None) and cell is None:
            raise NotEnoughParameters("ORCI", [("a", a), ("b", b), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCI":
                raise CellTypeMismatch(
                    "ORCI", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            conv_cell = [[0, 1.0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)
            if b < a - eps:
                cell = [cell[1], cell[0], cell[2]]
            if c < a - eps:
                cell = [cell[2], cell[0], cell[1]]
            elif c < b - eps:
                cell = [cell[0], cell[2], cell[1]]

            super().__init__(cell)
            self.conv_cell = [
                [0, 1.0, 1.0],
                [1.0, 0, 1.0],
                [1.0, 1.0, 0],
            ] @ self.cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
        else:
            a, b, c = tuple(sorted([a, b, c]))
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            super().__init__(
                [-a / 2, b / 2, c / 2], [a / 2, -b / 2, c / 2], [a / 2, b / 2, -c / 2]
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCI":
                raise CellTypeMismatch(
                    "ORCI", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            self.conv_cell = np.diag([self.conv_a, self.conv_b, self.conv_c])

        zeta = (1 + self.conv_a**2 / self.conv_c**2) / 4
        eta = (1 + self.conv_b**2 / self.conv_c**2) / 4
        delta = (self.conv_b**2 - self.conv_a**2) / (4 * self.conv_c**2)
        mu = (self.conv_a**2 + self.conv_b**2) / (4 * self.conv_c**2)

        self.points = {
            "G": np.array([0, 0, 0]),
            "L": np.array([-mu, mu, 1 / 2 - delta]),
            "L1": np.array([mu, -mu, 1 / 2 + delta]),
            "L2": np.array([1 / 2 - delta, 1 / 2 + delta, -mu]),
            "R": np.array([0, 1 / 2, 0]),
            "S": np.array([1 / 2, 0, 0]),
            "T": np.array([0, 0, 1 / 2]),
            "W": np.array([1 / 4, 1 / 4, 1 / 4]),
            "X": np.array([-zeta, zeta, zeta]),
            "X1": np.array([zeta, 1 - zeta, -zeta]),
            "Y": np.array([eta, -eta, eta]),
            "Y1": np.array([1 - eta, eta, -eta]),
            "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
        }

        self._default_path = [
            ["G", "X", "L", "T", "W", "R", "X1", "Z", "G", "Y", "S", "W"],
            ["L1", "Y"],
            ["Y1", "Z"],
        ]


# 9
class ORCC(Lattice):
    r"""
    C-centred orthorhombic (ORCC, oS)

    :math:`a < b`

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (a/2, -b/2, 0)

        \boldsymbol{a}_2 = (a/2, b/2, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    .. seealso::

        :ref:`lattice-orcc` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "oS"

    def __init__(
        self, a: float = None, b: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or b is None or c is None) and cell is None:
            raise NotEnoughParameters("ORCC", [("a", a), ("b", b), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCC":
                raise CellTypeMismatch(
                    "ORCC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

            # a == c
            if not (a < c - eps or c < a - eps):
                cell = [cell[0], cell[2], cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)
            # b = c
            elif not (b < c - eps or c < b - eps):
                cell = [cell[1], cell[2], cell[0]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            conv_cell = [[1.0, 1.0, 0], [-1.0, 1.0, 0], [0, 0, 1.0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)
            # b < a
            if b < a - eps:
                cell = [cell[1], cell[0], cell[2]]

            super().__init__(cell)
            self.conv_cell = [[1.0, 1.0, 0], [-1.0, 1.0, 0], [0, 0, 1.0]] @ self.cell
            a, b, c, alpha, beta, gamma = param_from_cell(self.conv_cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
        else:
            a, b = tuple(sorted([a, b]))
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            super().__init__([a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c])
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "ORCC":
                raise CellTypeMismatch(
                    "ORCC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            self.conv_cell = np.diag([self.conv_a, self.conv_b, self.conv_c])

        zeta = (1 + self.conv_a**2 / self.conv_b**2) / 4

        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([zeta, zeta, 1 / 2]),
            "A1": np.array([-zeta, 1 - zeta, 1 / 2]),
            "R": np.array([0, 1 / 2, 1 / 2]),
            "S": np.array([0, 1 / 2, 0]),
            "T": np.array([-1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([zeta, zeta, 0]),
            "X1": np.array([-zeta, 1 - zeta, 0]),
            "Y": np.array([-1 / 2, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        self._default_path = [
            ["G", "X", "S", "R", "A", "Z", "G", "Y", "X1", "A1", "T", "Y"],
            ["Z", "T"],
        ]


# 10
class HEX(Lattice):
    r"""
    Hexagonal (HEX, hP)

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a/2, -a\sqrt{3}, 0)

        \boldsymbol{a}_2 = (a/2, a\sqrt{3}, 0)

        \boldsymbol{a}_3 = (0, 0, c)

    .. seealso::

        :ref:`lattice-hex` for more information.

    Parameters
    ----------
    a : float, None
        Length of the lattice vector of the conventional lattice.
    c : float, None
        Length of the lattice vector of the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "hP"

    def __init__(
        self, a: float = None, c: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or c is None) and cell is None:
            raise NotEnoughParameters("HEX", [("a", a), ("c", c)])
        if cell is not None:
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "HEX":
                raise CellTypeMismatch(
                    "HEX", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

            if not (a < c - eps or c < a - eps):
                cell = [cell[0], cell[2], cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)
            elif not (b < c - eps or c < b - eps):
                cell = [cell[1], cell[2], cell[0]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            super().__init__(cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
        else:
            self.conv_a = a
            self.conv_c = c
            super().__init__(
                [a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "HEX":
                raise CellTypeMismatch(
                    "HEX", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

        self.conv_cell = self.cell

        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([0, 0, 1 / 2]),
            "H": np.array([1 / 3, 1 / 3, 1 / 2]),
            "K": np.array([1 / 3, 1 / 3, 0]),
            "L": np.array([1 / 2, 0, 1 / 2]),
            "M": np.array([1 / 2, 0, 0]),
        }

        self._default_path = [
            ["G", "M", "K", "G", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]


# 11
class RHL(Lattice):
    r"""
    Rhombohedral (RHL, hR)

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a\cos(\alpha / 2), -a\sin(\alpha/2), 0)

        \boldsymbol{a}_2 = (a\cos(\alpha / 2), a\sin(\alpha/2), 0)

        \boldsymbol{a}_3 = (\frac{\cos(\alpha)}{\cos(\alpha/2)}, 0, a\sqrt{1 - \frac{\cos^2(\alpha)}{\cos^2(\alpha/2)}})

    .. seealso::

        :ref:`lattice-rhl` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    alpha : float, optional
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "hR"

    def __init__(
        self, a: float = None, alpha: float = None, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or alpha is None) and cell is None:
            raise NotEnoughParameters("RHL", [("a", a), ("alpha", alpha)])
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "RHL":
                raise CellTypeMismatch(
                    "RHL", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

            super().__init__(cell)
            self.conv_a = (a + b + c) / 3
            self.conv_alpha = (alpha + beta + gamma) / 3
        else:
            if alpha >= 120:
                raise ValueError("alpha has to be < 120 degrees.")
            self.conv_a = a
            self.conv_alpha = alpha
            super().__init__(
                [a * cos(alpha / 180 * pi / 2), -a * sin(alpha / 180 * pi / 2), 0],
                [a * cos(alpha / 180 * pi / 2), a * sin(alpha / 180 * pi / 2), 0],
                [
                    a * cos(alpha / 180 * pi) / cos(alpha / 180 * pi / 2),
                    0,
                    a
                    * sqrt(
                        1 - cos(alpha / 180 * pi) ** 2 / cos(alpha / 180 * pi / 2) ** 2
                    ),
                ],
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "RHL":
                raise CellTypeMismatch(
                    "RHL", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

        self.conv_cell = self.cell
        if self.variation == "RHL1":
            eta = (1 + 4 * cos(alpha * toradians)) / (2 + 4 * cos(alpha * toradians))
            nu = 3 / 4 - eta / 2

            self.points = {
                "G": np.array([0, 0, 0]),
                "B": np.array([eta, 1 / 2, 1 - eta]),
                "B1": np.array([1 / 2, 1 - eta, eta - 1]),
                "F": np.array([1 / 2, 1 / 2, 0]),
                "L": np.array([1 / 2, 0, 0]),
                "L1": np.array([0, 0, -1 / 2]),
                "P": np.array([eta, nu, nu]),
                "P1": np.array([1 - nu, 1 - nu, 1 - eta]),
                "P2": np.array([nu, nu, eta - 1]),
                "Q": np.array([1 - nu, nu, 0]),
                "X": np.array([nu, 0, -nu]),
                "Z": np.array([1 / 2, 1 / 2, 1 / 2]),
            }

            self._default_path = [
                ["G", "L", "B1"],
                ["B", "Z", "G", "X"],
                ["Q", "F", "P1", "Z"],
                ["L", "P"],
            ]
        elif self.variation == "RHL2":
            eta = 1 / (2 * tan(alpha * toradians / 2) ** 2)
            nu = 3 / 4 - eta / 2

            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([1 / 2, -1 / 2, 0]),
                "L": np.array([1 / 2, 0, 0]),
                "P": np.array([1 - nu, -nu, 1 - nu]),
                "P1": np.array([nu, nu - 1, nu - 1]),
                "Q": np.array([eta, eta, eta]),
                "Q1": np.array([1 - eta, -eta, -eta]),
                "Z": np.array([1 / 2, -1 / 2, 1 / 2]),
            }

            self._default_path = [["G", "P", "Z", "Q", "G", "F", "P1", "Q1", "L", "Z"]]

    @property
    def variation(self):
        r"""
        Two variations of the Lattice.

        :math:`\text{RHL}_1 \alpha < 90^{\circ}`,
        :math:`\text{RHL}_2 \alpha > 90^{\circ}`

        Returns
        -------
        variation : str
            Variation of the lattice. Either "RHL1" or "RHL2".
        """
        if self.conv_alpha < 90:
            return "RHL1"
        elif self.conv_alpha > 90:
            return "RHL2"


# 12
class MCL(Lattice):
    r"""
    Monoclinic (MCL, mP)

    :math:`a, b \le c`, :math:`\alpha < 90^{\circ}`, :math:`\beta = \gamma = 90^{\circ}`.

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))

    .. seealso::

        :ref:`lattice-mcl` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    alpha : float, optional
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "mP"

    def __init__(
        self, a: float, b: float, c: float, alpha: float, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or alpha is None) and cell is None:
            raise NotEnoughParameters(
                "MCL", [("a", a), ("b", b), ("c", c), ("alpha", alpha)]
            )
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "MCL":
                raise CellTypeMismatch(
                    "MCL", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            eps = eps_rel * volume(cell) ** (1 / 3.0)

            # beta != 90
            if cos(beta * toradians) + eps < 0 or 0 < cos(beta * toradians) - eps:
                cell = [cell[1], cell[2], cell[0]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)
            # gamma != 90
            elif cos(gamma * toradians) + eps < 0 or 0 < cos(gamma * toradians) - eps:
                cell = [cell[2], cell[0], cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            # alpha > 90
            if cos(alpha * toradians) < -eps:
                cell = [cell[0], cell[2], -cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            # b > c
            if c < b - eps:
                cell = [-cell[0], cell[2], cell[1]]
                a, b, c, alpha, beta, gamma = param_from_cell(cell)

            super().__init__(cell)
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            self.conv_alpha = alpha
        else:
            eps = eps_rel * volume(a, b, c, alpha, 90, 90) ** (1 / 3.0)
            b, c = tuple(sorted([b, c]))
            if cos(alpha * toradians) < -eps:
                alpha = alpha - 90
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            self.conv_alpha = alpha
            super().__init__(
                [a, 0, 0],
                [0, b, 0],
                [0, c * cos(alpha * toradians), c * sin(alpha * toradians)],
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "MCL":
                raise CellTypeMismatch(
                    "MCL", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

        self.conv_cell = self.cell

        eta = (1 - b * cos(alpha * toradians) / c) / (2 * sin(alpha * toradians) ** 2)
        nu = 1 / 2 - eta * c * cos(alpha * toradians) / b

        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2, 0]),
            "C": np.array([0, 1 / 2, 1 / 2]),
            "D": np.array([1 / 2, 0, 1 / 2]),
            "D1": np.array([1 / 2, 0, -1 / 2]),
            "E": np.array([1 / 2, 1 / 2, 1 / 2]),
            "H": np.array([0, eta, 1 - nu]),
            "H1": np.array([0, 1 - eta, nu]),
            "H2": np.array([0, eta, -nu]),
            "M": np.array([1 / 2, eta, 1 - nu]),
            "M1": np.array([1 / 2, 1 - eta, nu]),
            "M2": np.array([1 / 2, eta, -nu]),
            "X": np.array([0, 1 / 2, 0]),
            "Y": np.array([0, 0, 1 / 2]),
            "Y1": np.array([0, 0, -1 / 2]),
            "Z": np.array([1 / 2, 0, 0]),
        }

        self._default_path = [
            ["G", "Y", "H", "C", "E", "M1", "A", "X", "H1"],
            ["M", "D", "Z"],
            ["Y", "D"],
        ]


# 13
class MCLC(Lattice):
    r"""
    C-centred monoclinic (MCLC, mS)

    :math:`a, b \le c`, :math:`\alpha < 90^{\circ}`, :math:`\beta = \gamma = 90^{\circ}`.

    Conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))

    Primitive lattice:

    .. math::

        \boldsymbol{a}_1 = (a/2, b/2, 0)

        \boldsymbol{a}_2 = (-a/2, b/2, 0)

        \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))

    .. seealso::

        :ref:`lattice-mclc` for more information.

    Parameters
    ----------
    a : float, optional
        Length of the lattice vector of the conventional lattice.
    b : float, optional
        Length of the lattice vector of the conventional lattice.
    c : float, optional
        Length of the lattice vector of the conventional lattice.
    alpha : float, optional
        Angle between b and c. In degrees. Corresponds to the conventional lattice
    cell : (3,3) |array_like|_, optional
        Primitive cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [1]_.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    """

    _pearson_symbol = "mS"

    def __init__(
        self, a: float, b: float, c: float, alpha: float, cell=None, eps_rel=1e-5
    ) -> None:
        if (a is None or alpha is None) and cell is None:
            raise NotEnoughParameters(
                "MCLC", [("a", a), ("b", b), ("c", c), ("alpha", alpha)]
            )
        if cell is not None:
            a, b, c, alpha, beta, gamma = param_from_cell(cell)
            eps = eps_rel * volume(cell) ** (1 / 3.0)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "MCLC":
                raise CellTypeMismatch(
                    "MCLC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )
            eps = eps_rel * volume(cell) ** (1 / 3.0)

            # a == c
            if not (a < c - eps or c < a - eps):
                cell = [cell[2], cell[1], cell[1]]
            # b == c
            elif not (b < c - eps or c < b - eps):
                cell = [cell[1], cell[2], cell[0]]

            conv_cell = [[1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]] @ cell
            a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)

            # alpha > 90
            if cos(alpha * toradians) < -eps:
                cell = [cell[0], cell[2], -cell[1]]
                conv_cell = [[1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]] @ cell
                a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)

            # b > c
            if c < b - eps:
                cell = [-cell[0], cell[2], cell[1]]
                conv_cell = [[1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]] @ cell
                a, b, c, alpha, beta, gamma = param_from_cell(conv_cell)

            super().__init__(cell)
            self.conv_cell = [
                [1.0, -1.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ] @ self.cell
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            self.conv_alpha = alpha
        else:
            eps = eps_rel * volume(a, b, c, alpha, 90, 90) ** (1 / 3.0)
            b, c = tuple(sorted([b, c]))
            if cos(alpha * toradians) < -eps:
                alpha = alpha - 90
            self.conv_a = a
            self.conv_b = b
            self.conv_c = c
            self.conv_alpha = alpha
            super().__init__(
                [a / 2, b / 2, 0],
                [-a / 2, b / 2, 0],
                [
                    0,
                    c * cos(alpha * toradians),
                    c * sin(alpha * toradians),
                ],
            )
            self.conv_cell = np.array(
                [
                    [self.conv_a, 0, 0],
                    [0, self.conv_b, 0],
                    [
                        0,
                        self.conv_c * cos(alpha * toradians),
                        self.conv_c * sin(alpha * toradians),
                    ],
                ]
            )
            a, b, c, alpha, beta, gamma = param_from_cell(self.cell)
            lattice_type = lepage(a, b, c, alpha, beta, gamma, eps_rel=eps_rel)
            if lattice_type != "MCLC":
                raise CellTypeMismatch(
                    "MCLC", eps_rel, a, b, c, alpha, beta, gamma, lattice_type
                )

        # Parameters
        if self.variation in ["MCLC1", "MCLC2"]:
            zeta = (
                2 - self.conv_b * cos(self.conv_alpha * toradians) / self.conv_c
            ) / (4 * sin(self.conv_alpha * toradians) ** 2)
            eta = (
                1 / 2
                + 2
                * zeta
                * self.conv_c
                * cos(self.conv_alpha * toradians)
                / self.conv_b
            )
            psi = 3 / 4 - self.conv_a**2 / (
                4 * self.conv_b**2 * sin(self.conv_alpha * toradians) ** 2
            )
            phi = (
                psi + (3 / 4 - psi) * self.conv_b * cos(self.conv_alpha * toradians) / c
            )
        elif self.variation in ["MCLC3", "MCLC4"]:
            mu = (1 + self.conv_b**2 / self.conv_a**2) / 4
            delta = (
                self.conv_b
                * c
                * cos(self.conv_alpha * toradians)
                / (2 * self.conv_a**2)
            )
            zeta = (
                mu
                - 1 / 4
                + (1 - self.conv_b * cos(self.conv_alpha * toradians) / self.conv_c)
                / (4 * sin(self.conv_alpha * toradians) ** 2)
            )
            eta = (
                1 / 2
                + 2
                * zeta
                * self.conv_c
                * cos(self.conv_alpha * toradians)
                / self.conv_b
            )
            phi = 1 + zeta - 2 * mu
            psi = eta - 2 * delta
        elif self.variation == "MCLC5":
            zeta = (
                self.conv_b**2 / self.conv_a**2
                + (1 - self.conv_b * cos(self.conv_alpha * toradians) / self.conv_c)
                / sin(self.conv_alpha * toradians) ** 2
            ) / 4
            eta = (
                1 / 2
                + 2
                * zeta
                * self.conv_c
                * cos(self.conv_alpha * toradians)
                / self.conv_b
            )
            mu = (
                eta / 2
                + self.conv_b**2 / (4 * self.conv_a**2)
                - self.conv_b
                * self.conv_c
                * cos(self.conv_alpha * toradians)
                / (2 * self.conv_a**2)
            )
            nu = 2 * mu - zeta
            rho = 1 - zeta * self.conv_a**2 / self.conv_b**2
            omega = (
                (
                    4 * nu
                    - 1
                    - self.conv_b**2
                    * sin(self.conv_alpha * toradians) ** 2
                    / self.conv_a**2
                )
                * self.conv_c
                / (2 * self.conv_b * cos(self.conv_alpha * toradians))
            )
            delta = (
                zeta * self.conv_c * cos(self.conv_alpha * toradians) / self.conv_b
                + omega / 2
                - 1 / 4
            )

        # Path
        if self.variation == "MCLC1":
            self._default_path = [
                ["G", "Y", "F", "L", "I"],
                ["I1", "Z", "F1"],
                ["Y", "X1"],
                ["X", "G", "N"],
                ["M", "G"],
            ]
            self.points = {
                "G": np.array([0, 0, 0]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
                "F1": np.array([zeta, zeta, eta]),
                "F2": np.array([-zeta, -zeta, 1 - eta]),
                "I": np.array([phi, 1 - phi, 1 / 2]),
                "I1": np.array([1 - phi, phi - 1, 1 / 2]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "M": np.array([1 / 2, 0, 1 / 2]),
                "X": np.array([1 - psi, psi - 1, 0]),
                "X1": np.array([psi, 1 - psi, 0]),
                "X2": np.array([psi - 1, -psi, 0]),
                "Y": np.array([1 / 2, 1 / 2, 0]),
                "Y1": np.array([-1 / 2, -1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }
        elif self.variation == "MCLC2":
            self._default_path = [
                ["G", "Y", "F", "L", "I"],
                ["I1", "Z", "F1"],
                ["N", "G", "M"],
            ]
            self.points = {
                "G": np.array([0, 0, 0]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
                "F1": np.array([zeta, zeta, eta]),
                "F2": np.array([-zeta, -zeta, 1 - eta]),
                "F3": np.array([1 - zeta, -zeta, 1 - eta]),
                "I": np.array([phi, 1 - phi, 1 / 2]),
                "I1": np.array([1 - phi, phi - 1, 1 / 2]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "M": np.array([1 / 2, 0, 1 / 2]),
                "X": np.array([1 - psi, psi - 1, 0]),
                "Y": np.array([1 / 2, 1 / 2, 0]),
                "Y1": np.array([-1 / 2, -1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }
        elif self.variation == "MCLC3":
            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([1 - phi, 1 - phi, 1 - psi]),
                "F1": np.array([phi, phi - 1, psi]),
                "F2": np.array([1 - phi, -phi, 1 - psi]),
                "H": np.array([zeta, zeta, eta]),
                "H1": np.array([1 - zeta, -zeta, 1 - eta]),
                "H2": np.array([-zeta, -zeta, 1 - eta]),
                "I": np.array([1 / 2, -1 / 2, 1 / 2]),
                "M": np.array([1 / 2, 0, 1 / 2]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "X": np.array([1 / 2, -1 / 2, 0]),
                "Y": np.array([mu, mu, delta]),
                "Y1": np.array([1 - mu, -mu, -delta]),
                "Y2": np.array([-mu, -mu, -delta]),
                "Y3": np.array([mu, mu - 1, delta]),
                "Z": np.array([0, 0, 1 / 2]),
            }
            self._default_path = [
                ["G", "Y", "F", "H", "Z", "I", "F1"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
        elif self.variation == "MCLC4":
            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([1 - phi, 1 - phi, 1 - psi]),
                "H": np.array([zeta, zeta, eta]),
                "H1": np.array([1 - zeta, -zeta, 1 - eta]),
                "H2": np.array([-zeta, -zeta, 1 - eta]),
                "I": np.array([1 / 2, -1 / 2, 1 / 2]),
                "M": np.array([1 / 2, 0, 1 / 2]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "X": np.array([1 / 2, -1 / 2, 0]),
                "Y": np.array([mu, mu, delta]),
                "Y1": np.array([1 - mu, -mu, -delta]),
                "Y2": np.array([-mu, -mu, -delta]),
                "Y3": np.array([mu, mu - 1, delta]),
                "Z": np.array([0, 0, 1 / 2]),
            }
            self._default_path = [
                ["G", "Y", "F", "H", "Z", "I"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
        elif self.variation == "MCLC5":
            self._default_path = [
                ["G", "Y", "F", "L", "I"],
                ["I1", "Z", "H", "F1"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([nu, nu, omega]),
                "F1": np.array([1 - nu, 1 - nu, 1 - omega]),
                "F2": np.array([nu, nu - 1, omega]),
                "H": np.array([zeta, zeta, eta]),
                "H1": np.array([1 - zeta, -zeta, 1 - eta]),
                "H2": np.array([-zeta, -zeta, 1 - eta]),
                "I": np.array([rho, 1 - rho, 1 / 2]),
                "I1": np.array([1 - rho, rho - 1, 1 / 2]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "M": np.array([1 / 2, 0, 1 / 2]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "X": np.array([1 / 2, -1 / 2, 0]),
                "Y": np.array([mu, mu, delta]),
                "Y1": np.array([1 - mu, -mu, -delta]),
                "Y2": np.array([-mu, -mu, -delta]),
                "Y3": np.array([mu, mu - 1, delta]),
                "Z": np.array([0, 0, 1 / 2]),
            }

    @property
    def variation(self):
        r"""
        Five variation of the Lattice.

        :math:`\text{MCLC}_1: k_{\gamma} > 90^{\circ}`,
        :math:`\text{MCLC}_2: k_{\gamma} = 90^{\circ}`,
        :math:`\text{MCLC}_3: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`
        :math:`\text{MCLC}_4: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`
        :math:`\text{MCLC}_5: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`

        Returns
        -------
        variation : str
            Variation of the lattice.
            Either "MCLC1", "MCLC2", "MCLC3", "MCLC4" or "MCLC5".
        """

        if np.abs(self.k_gamma - 90) <= 10 * np.finfo(float).eps:
            return "MCLC2"
        elif self.k_gamma > 90:
            return "MCLC1"
        elif self.k_gamma < 90:
            expression = (
                self.conv_b * cos(self.conv_alpha * toradians) / self.conv_c
                + self.conv_b**2
                * sin(self.conv_alpha * toradians) ** 2
                / self.conv_a**2
            )
            if np.abs(expression - 1) <= 10 * np.finfo(float).eps:
                return "MCLC4"
            elif expression < 1:
                return "MCLC3"
            elif expression > 1:
                return "MCLC5"


# 14
class TRI(Lattice):
    r"""
    Triclinic (TRI, aP)



    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (b\cos(\gamma), b\sin(\gamma), 0)

        \boldsymbol{a}_3 = (c\cos(\beta), \frac{c(\cos(\alpha) - \cos(\beta)\cos(\gamma))}{\sin{\gamma}}, \frac{c}{\sin(\gamma)}\sqrt{\sin^2(\gamma) - \cos^2(\alpha) - \cos^2(\beta) + 2\cos(\alpha)\cos(\beta)\cos(\gamma)})

    Variations of the trigonal lattice are defined through the angles of the reciprocal cell,
    therefore it is possible to define trigonal Bravais lattice with reciprocal cell parameters
    (argument ``reciprocal``).

    .. seealso::

        :ref:`lattice-tri` for more information.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.
    alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    beta : float
        Angle between a and c. In degrees. Corresponds to the conventional lattice.
    gamma : float
        Angle between a and b. In degrees. Corresponds to the conventional lattice.
    reciprocal : bool, default False
        Whether to interpret ``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma``
        as reciprocal lattice parameters.


    Attributes
    ----------
    conv_a : float
        Length of the lattice vector of the conventional lattice.
    conv_b : float
        Length of the lattice vector of the conventional lattice.
    conv_c : float
        Length of the lattice vector of the conventional lattice.
    conv_alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice.
    beta : float
        Angle between a and c. In degrees. Corresponds to the conventional lattice.
    conv_gamma : float
        Angle between a and b. In degrees. Corresponds to the conventional lattice.
    conv_cell : (3,3) :numpy:`ndarray`
        Conventional unit cell.

        .. code-block:: python

            conv_cell = [[a_x, a_y, a_z],
                         [b_x, b_y, b_z],
                         [c_x, c_y, c_z]]

    """

    _pearson_symbol = "aP"

    def __init__(
        self,
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
        reciprocal=False,
    ) -> None:
        if not reciprocal:
            tmp = sorted([(a, alpha), (b, beta), (c, gamma)], key=lambda x: x[0])
            a = tmp[0][0]
            alpha = tmp[0][1]
            b = tmp[1][0]
            beta = tmp[1][1]
            c = tmp[2][0]
            gamma = tmp[2][1]
        super().__init__(a, b, c, alpha, beta, gamma)
        if reciprocal:
            self.cell = self.reciprocal_cell
        self.conv_a = a
        self.conv_b = b
        self.conv_c = c
        self.conv_alpha = alpha
        self.conv_beta = beta
        self.conv_gamma = gamma
        if (
            self.k_alpha < self.k_gamma < self.k_beta
            or self.k_beta < self.k_gamma < self.k_alpha
        ):
            raise RuntimeError("k_gamma is not minimal nor maximal.")
        self.conv_cell = self.cell
        if self.variation in ["TRI1a", "TRI2a"]:
            self.points = {
                "G": np.array([0, 0, 0]),
                "L": np.array([1 / 2, 1 / 2, 0]),
                "M": np.array([0, 1 / 2, 1 / 2]),
                "N": np.array([1 / 2, 0, 1 / 2]),
                "R": np.array([1 / 2, 1 / 2, 1 / 2]),
                "X": np.array([1 / 2, 0, 0]),
                "Y": np.array([0, 1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }

            self._default_path = [
                ["X", "G", "Y"],
                ["L", "G", "Z"],
                ["N", "G", "M"],
                "R",
                "G",
            ]
        elif self.variation in ["TRI1b", "TRI2b"]:
            self.points = {
                "G": np.array([0, 0, 0]),
                "L": np.array([1 / 2, -1 / 2, 0]),
                "M": np.array([0, 0, 1 / 2]),
                "N": np.array([-1 / 2, -1 / 2, 1 / 2]),
                "R": np.array([0, -1 / 2, 1 / 2]),
                "X": np.array([0, -1 / 2, 0]),
                "Y": np.array([1 / 2, 0, 0]),
                "Z": np.array([-1 / 2, 0, 1 / 2]),
            }

            self._default_path = [
                ["X", "G", "Y"],
                ["L", "G", "Z"],
                ["N", "G", "M"],
                "R",
                "G",
            ]

    @property
    def variation(self):
        r"""
        Four variations of the Lattice.

        :math:`\text{TRI}_{1a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`

        :math:`\text{TRI}_{1b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`

        :math:`\text{TRI}_{2a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`

        :math:`\text{TRI}_{2b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`

        Returns
        -------
        variation : str
            Variation of the lattice.
            Either "TRI1a", "TRI1b", "TRI2a" or "TRI2b".
        """
        if self.k_gamma == 90:
            if self.k_alpha > 90 and self.k_beta > 90:
                return "TRI2a"
            elif self.k_alpha < 90 and self.k_beta < 90:
                return "TRI2b"
        elif (min(self.k_gamma, self.k_beta, self.k_alpha)) > 90:
            return "TRI1a"
        elif (max(self.k_gamma, self.k_beta, self.k_alpha)) < 90:
            return "TRI1b"
        else:
            return "TRI"


def bravais_lattice_from_param(a, b, c, alpha, beta, gamma) -> Lattice:
    r"""
    Create Bravais lattice from lattice parameters.

    Orientation is default as described in [1]_.

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
    lattice_type : str
        Lattice type.

    Returns
    -------
    bravais_lattice : Lattice
        Bravais lattice.

    References
    ----------
    .. [1] Setyawan, W. and Curtarolo, S., 2010.
        High-throughput electronic band structure calculations: Challenges and tools.
        Computational materials science, 49(2), pp.299-312.
    """

    lattice_type = lepage(a, b, c, alpha, beta, gamma)

    if lattice_type == "CUB":
        a = (a + b + c) / 3
        return CUB(a)
    if lattice_type == "FCC":
        a = (a + b + c) / 3
        return FCC(a)
    if lattice_type == "BCC":
        a = (a + b + c) / 3
        return BCC(a)
    if lattice_type == "TET":
        if a == b:
            return TET((a + b) / 2, c)
        elif a == c:
            return TET((a + c) / 2, b)
        elif b == c:
            return TET((b + c) / 2, a)
    if lattice_type == "BCT":
        if a == b:
            return BCT((a + b) / 2, c)
        elif a == c:
            return BCT((a + c) / 2, b)
        elif b == c:
            return BCT((b + c) / 2, a)
    if lattice_type == "ORC":
        return ORC(a, b, c)
    if lattice_type == "ORCF":
        return ORCF(a, b, c)
    if lattice_type == "ORCC":
        return ORCC(a, b, c)
    if lattice_type == "ORCI":
        return ORCI(a, b, c)
    if lattice_type == "HEX":
        if abs(a - b) < abs(a - c) and abs(a - b) < abs(b - c):
            return HEX((a + b) / 2, c)
        elif abs(a - c) < abs(a - b) and abs(a - c) < abs(b - c):
            return HEX((a + c) / 2, b)
        elif abs(b - c) < abs(a - c) and abs(b - c) < abs(a - b):
            return HEX((b + c) / 2, a)
    if lattice_type == "RHL":
        return RHL((a + b + c) / 3, (alpha + beta + gamma) / 2)
    if lattice_type == "MCL":
        alpha, beta, gamma = [alpha, beta, gamma].sort()
        return MCL(a, b, c, alpha)
    if lattice_type == "MCLC":
        alpha, beta, gamma = [alpha, beta, gamma].sort()
        return MCLC(a, b, c, alpha)
    return TRI(a, b, c, alpha, beta, gamma)


def bravais_lattice_from_cell(cell) -> Lattice:
    r"""
    Create Bravais lattice from cell matrix.

    Orientation of the cell is respected, however the lattice vectors are renamed
    with respect to [1]_.

    Parameters
    ----------
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]

    lattice_type : str
        Lattice type.

    Returns
    -------
    bravais_lattice : Lattice
        Bravais lattice.

    References
    ----------
    .. [1] Setyawan, W. and Curtarolo, S., 2010.
        High-throughput electronic band structure calculations: Challenges and tools.
        Computational materials science, 49(2), pp.299-312.
    """

    lattice_type = lepage(*param_from_cell(cell))

    if lattice_type == "CUB":
        return CUB(cell=cell)
    if lattice_type == "FCC":
        return FCC(cell=cell)
    if lattice_type == "BCC":
        return BCC(cell=cell)
    if lattice_type == "TET":
        return TET(cell=cell)
    if lattice_type == "BCT":
        return BCT(cell=cell)
    if lattice_type == "ORC":
        return ORC(cell=cell)
    if lattice_type == "ORCF":
        return ORCF(cell=cell)
    if lattice_type == "ORCC":
        return ORCC(cell=cell)
    if lattice_type == "ORCI":
        return ORCI(cell=cell)
    if lattice_type == "HEX":
        return HEX(cell=cell)
    if lattice_type == "RHL":
        return RHL(cell=cell)
    if lattice_type == "MCL":
        return MCL(cell=cell)
    if lattice_type == "MCLC":
        return MCLC(cell=cell)
    return TRI(cell=cell)
    return Lattice(cell)


def lattice_example(
    lattice=None,
):
    r"""
    Return an example of the lattice.

    Parameters
    ----------
    lattice : str, optional
        Name of the lattice to be returned.
        For available names see documentation of each Bravais lattice class.
        Lowercased before usage.

    Returns
    -------
    lattice : Lattice or list   
        Child of the :py:class:`.Lattice` class is returned.
        If no math found a list with available examples is returned.
    """

    all_examples = [
        "CUB",
        "FCC",
        "BCC",
        "TET",
        "BCT1",
        "BCT2",
        "ORC",
        "ORCF1",
        "ORCF2",
        "ORCF3",
        "ORCI",
        "ORCC",
        "HEX",
        "RHL1",
        "RHL2",
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
    if not isinstance(lattice, str):
        return all_examples

    lattice = lattice.lower()

    if lattice == "cub":
        return CUB(pi)
    elif lattice == "fcc":
        return FCC(pi)
    elif lattice == "bcc":
        return BCC(pi)
    elif lattice == "tet":
        return TET(pi, 1.5 * pi)
    elif lattice in ["bct1", "bct"]:
        return BCT(1.5 * pi, pi)
    elif lattice == "bct2":
        return BCT(pi, 1.5 * pi)
    elif lattice == "orc":
        return ORC(pi, 1.5 * pi, 2 * pi)
    elif lattice in ["orcf1", "orcf"]:
        return ORCF(0.7 * pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orcf2":
        return ORCF(1.2 * pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orcf3":
        return ORCF(pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orci":
        return ORCI(pi, 1.3 * pi, 1.7 * pi)
    elif lattice == "orcc":
        return ORCC(pi, 1.3 * pi, 1.7 * pi)
    elif lattice == "hex":
        return HEX(pi, 2 * pi)
    elif lattice in ["rhl1", "rhl"]:
        # If alpha = 60 it is effectively FCC!
        return RHL(pi, 70)
    elif lattice == "rhl2":
        return RHL(pi, 110)
    elif lattice == "mcl":
        return MCL(pi, 1.3 * pi, 1.6 * pi, alpha=75)
    elif lattice in ["mclc1", "mclc"]:
        return MCLC(pi, 1.4 * pi, 1.7 * pi, 80)
    elif lattice == "mclc2":
        return MCLC(1.4 * pi * sin(75 * toradians), 1.4 * pi, 1.7 * pi, 75)
    elif lattice == "mclc3":
        b = pi
        x = 1.1
        alpha = 78
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha) * 1.8
        a = x * b * sin(ralpha)
        return MCLC(a, b, c, alpha)
    elif lattice == "mclc4":
        b = pi
        x = 1.2
        alpha = 65
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha)
        a = x * b * sin(ralpha)
        return MCLC(a, b, c, alpha)
    elif lattice == "mclc5":
        b = pi
        x = 1.4
        alpha = 53
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha) * 0.9
        a = x * b * sin(ralpha)
        return MCLC(a, b, c, alpha)
    elif lattice in ["tri1a", "tri1", "tri", "tria"]:
        return TRI(1, 1.5, 2, 120, 110, 100, reciprocal=True)
    elif lattice in ["tri2a", "tri2"]:
        return TRI(1, 1.5, 2, 120, 110, 90, reciprocal=True)
    elif lattice in ["tri1b", "trib"]:
        return TRI(1, 1.5, 2, 60, 70, 80, reciprocal=True)
    elif lattice == "tri2b":
        return TRI(1, 1.5, 2, 60, 70, 90, reciprocal=True)
    else:
        return all_examples


if __name__ == "__main__":
    for e in lattice_example()[-4:]:
        print(e)
        l = lattice_example(e)
        l.prepare_figure()
        l.plot("brillouin_kpath", label=l.variation)
        l.legend()
        l.show()
