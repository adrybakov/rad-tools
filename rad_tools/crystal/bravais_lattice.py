r"""
For each type of Bravais lattice a class defined, for some classes there are 
several variations of lattice, each of which are treated under same class 
(see :py:attr:`variation` for each class).

For each type and variation a predefined example of the lattice is available. 
It could be accessed in a following way:

.. doctest::

    >>> import rad_tools as rad_tools
    >>> cubic_example = rad.lattice_example("cub")

Each Bravais Lattice is created by the parameters 
:math:`a`, :math:`b`, :math:`c`, :math:`\alpha`, :math:`\beta`, :math:`\gamma`,
which corresponds to the conventional lattice. However, attributes of the class 
``self.a``, ``self.b``, ``self.c``, ``self.alpha``, ``self.beta``, ``self.gamma`` 
return parameters of the primitive lattice. Conventional lattice may be accessed through
the attributes ``self.conv_cell``, ``self.conv_a``, ``self.conv_b``, ``self.conv_c``, 
``self.conv_alpha``, ``self.conv_beta``, ``self.conv_gamma`` .
"""

from math import cos, pi, sin, sqrt, tan

import numpy as np

from rad_tools.crystal.lattice import Lattice
from rad_tools.routines import _toradians

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
]


# 1
class CUB(Lattice):
    r"""
    Cubic (CUB, cP)

    Primitive and conventional lattice:

    .. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, a, 0)

        \boldsymbol{a}_3 = (0, 0, a)

    Parameters
    ----------
    a : float
        Length of the lattice vectors of the conventional cel.

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

    """

    _pearson_symbol = "cP"

    def __init__(self, a: float) -> None:
        self.conv_a = a
        super().__init__([a, 0, 0], [0, a, 0], [0, 0, a])
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cel.

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

    """

    _pearson_symbol = "cF"

    def __init__(self, a: float) -> None:
        self.conv_a = a
        super().__init__([0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0])
        self.conv_cell = np.diag([a, a, a])
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.

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

    """

    _pearson_symbol = "cI"

    def __init__(self, a: float) -> None:
        self.conv_a = a
        super().__init__(
            [-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]
        )
        self.conv_cell = np.diag([a, a, a])
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.

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

    """

    _pearson_symbol = "tP"

    def __init__(self, a: float, c: float) -> None:
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of conventional lattice.
    c : float
        Length of the lattice vector of conventional lattice.

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

    """

    _pearson_symbol = "tI"

    def __init__(self, a: float, c: float) -> None:
        if a == c:
            raise ValueError("Are you trying to create BCC Lattice (a == c)?")
        self.conv_a = a
        self.conv_c = c
        super().__init__(
            [-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]
        )
        self._PLOT_NAMES["S"] = "$\\Sigma$"
        self._PLOT_NAMES["S1"] = "$\\Sigma_1$"
        self.conv_cell = np.diag([a, a, c])
        if self.variation == "BCT1":
            eta = (1 + c**2 / a**2) / 4
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
            eta = (1 + a**2 / c**2) / 4
            zeta = a**2 / (2 * c**2)
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.

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

    """

    _pearson_symbol = "oP"

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct CUB Lattice (a = b = c)?")
        if a == b or b == c:
            raise ValueError(
                "Are you trying to construct TET Lattice (a = b != c or a != b == c)?"
            )
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.

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

    """

    _pearson_symbol = "oF"

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct FCC Lattice (a = b = c)?")
        if a == b:
            raise ValueError("FIXME, dont know which lattice it will be.")
        if b == c:
            raise ValueError("FIXME, dont know which lattice it will be.")
        self.conv_a = a
        self.conv_b = b
        self.conv_c = c
        super().__init__([0, b / 2, c / 2], [a / 2, 0, c / 2], [a / 2, b / 2, 0])
        self.conv_cell = np.diag([a, b, c])
        if self.variation == "ORCF1":
            eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
            zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4
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
            eta = (1 + a**2 / b**2 - a**2 / c**2) / 4
            delta = (1 + b**2 / a**2 - b**2 / c**2) / 4
            phi = (1 + c**2 / b**2 - c**2 / a**2) / 4

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
            eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
            zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4

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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.


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

    """

    _pearson_symbol = "oI"

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct BCC Lattice (a = b = c)?")
        if a == b:
            raise ValueError("Are you trying to construct BCT2 Lattice (a = b < c)?")
        if b == c:
            raise ValueError("Are you trying to construct BCT1 Lattice (a < b = c)?")
        self.conv_a = a
        self.conv_b = b
        self.conv_c = c
        super().__init__(
            [-a / 2, b / 2, c / 2], [a / 2, -b / 2, c / 2], [a / 2, b / 2, -c / 2]
        )
        self.conv_cell = np.diag([a, b, c])
        zeta = (1 + a**2 / c**2) / 4
        eta = (1 + b**2 / c**2) / 4
        delta = (b**2 - a**2) / (4 * c**2)
        mu = (a**2 + b**2) / (4 * c**2)

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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.


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

    """

    _pearson_symbol = "oS"

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct TET Lattice (a = b = c)?")
        self.conv_a = a
        self.conv_b = c
        self.conv_c = b
        super().__init__([a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c])
        self.conv_cell = np.diag([a, b, c])
        zeta = (1 + a**2 / b**2) / 4

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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.


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

    """

    _pearson_symbol = "hP"

    def __init__(self, a: float, c: float) -> None:
        self.conv_a = a
        self.conv_c = c
        super().__init__(
            [a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]
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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice.


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

    """

    _pearson_symbol = "hR"

    def __init__(self, a: float, alpha: float) -> None:
        if alpha == 90:
            raise ValueError("Are you trying to construct CUB Lattice (alpha == 90)?")
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
                * sqrt(1 - cos(alpha / 180 * pi) ** 2 / cos(alpha / 180 * pi / 2) ** 2),
            ],
        )
        self.conv_cell = self.cell
        if self.variation == "RHL1":
            eta = (1 + 4 * cos(alpha / 180 * pi)) / (2 + 4 * cos(alpha / 180 * pi))
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
            eta = 1 / (2 * tan(alpha / 180 * pi / 2) ** 2)
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

    """

    _pearson_symbol = "mP"

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.conv_a = a
        self.conv_b = b
        self.conv_c = c
        self.conv_alpha = alpha
        super().__init__(
            [a, 0, 0],
            [0, b, 0],
            [0, c * cos(alpha / 180 * pi), c * sin(alpha / 180 * pi)],
        )
        self.conv_cell = self.cell

        eta = (1 - b * cos(alpha / 180 * pi) / c) / (2 * sin(alpha / 180 * pi) ** 2)
        nu = 1 / 2 - eta * c * cos(alpha / 180 * pi) / b

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

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional lattice.
    b : float
        Length of the lattice vector of the conventional lattice.
    c : float
        Length of the lattice vector of the conventional lattice.
    alpha : float
        Angle between b and c. In degrees. Corresponds to the conventional lattice


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

    """

    _pearson_symbol = "mS"

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        if a > c:
            a, c = c, a
        if b > c:
            b, c = c, b
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.conv_a = a
        self.conv_b = b
        self.conv_c = c
        self.conv_alpha = alpha
        super().__init__(
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [
                0,
                c * cos(alpha * _toradians),
                c * sin(alpha * _toradians),
            ],
        )
        self.conv_cell = np.array(
            [
                [a, 0, 0],
                [0, b, 0],
                [
                    0,
                    c * cos(alpha * _toradians),
                    c * sin(alpha * _toradians),
                ],
            ]
        )
        # Parameters
        if self.variation in ["MCLC1", "MCLC2"]:
            zeta = (2 - b * cos(alpha * _toradians) / c) / (
                4 * sin(alpha * _toradians) ** 2
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha * _toradians) / b
            psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha * _toradians) ** 2)
            phi = psi + (3 / 4 - psi) * b * cos(alpha * _toradians) / c
        elif self.variation in ["MCLC3", "MCLC4"]:
            mu = (1 + b**2 / a**2) / 4
            delta = b * c * cos(alpha * _toradians) / (2 * a**2)
            zeta = (
                mu
                - 1 / 4
                + (1 - b * cos(alpha * _toradians) / c)
                / (4 * sin(alpha * _toradians) ** 2)
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha * _toradians) / b
            phi = 1 + zeta - 2 * mu
            psi = eta - 2 * delta
        elif self.variation == "MCLC5":
            zeta = (
                b**2 / a**2
                + (1 - b * cos(alpha * _toradians) / c) / sin(alpha * _toradians) ** 2
            ) / 4
            eta = 1 / 2 + 2 * zeta * c * cos(alpha * _toradians) / b
            mu = (
                eta / 2
                + b**2 / (4 * a**2)
                - b * c * cos(alpha * _toradians) / (2 * a**2)
            )
            nu = 2 * mu - zeta
            rho = 1 - zeta * a**2 / b**2
            omega = (
                (4 * nu - 1 - b**2 * sin(alpha * _toradians) ** 2 / a**2)
                * c
                / (2 * b * cos(alpha * _toradians))
            )
            delta = zeta * c * cos(alpha * _toradians) / b + omega / 2 - 1 / 4

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
        """

        if np.abs(self.k_gamma - 90) <= 10 * np.finfo(float).eps:
            return "MCLC2"
        elif self.k_gamma > 90:
            return "MCLC1"
        elif self.k_gamma < 90:
            expression = (
                self.conv_b * cos(self.conv_alpha * _toradians) / self.conv_c
                + self.conv_b**2
                * sin(self.conv_alpha * _toradians) ** 2
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
        super().__init__(
            [a, 0, 0],
            [b * cos(gamma * _toradians), b * sin(gamma * _toradians), 0],
            [
                c * cos(beta * _toradians),
                c
                / sin(gamma * _toradians)
                * (
                    cos(alpha * _toradians)
                    - cos(beta * _toradians) * cos(gamma * _toradians)
                ),
                c
                / sin(gamma * _toradians)
                * sqrt(
                    sin(gamma * _toradians) ** 2
                    - cos(alpha * _toradians) ** 2
                    - cos(beta * _toradians) ** 2
                    + 2
                    * cos(alpha * _toradians)
                    * cos(beta * _toradians)
                    * cos(gamma * _toradians)
                ),
            ],
        )
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


def lattice_example(
    lattice=None,
):
    r"""
    Return an example of the lattice.

    Parameters
    ----------
    lattice : str, default None
        Name of the lattice to be returned.
        For available names see documentation of each Bravais lattice class.
        Lowercased before usage.

    Returns
    -------
    lattice
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
        return TET(pi, 2 * pi)
    elif lattice in ["bct1", "bct"]:
        return BCT(2 * pi, pi)
    elif lattice == "bct2":
        return BCT(pi, 2 * pi)
    elif lattice == "orc":
        return ORC(pi, 2 * pi, 3 * pi)
    elif lattice in ["orcf1", "orcf"]:
        return ORCF(0.9 * pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orcf2":
        return ORCF(1.1 * pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orcf3":
        return ORCF(pi, 5 / 4 * pi, 5 / 3 * pi)
    elif lattice == "orci":
        return ORCI(pi, 2 * pi, 3 * pi)
    elif lattice == "orcc":
        return ORCC(pi, 2 * pi, 3 * pi)
    elif lattice == "hex":
        return HEX(pi, 2 * pi)
    elif lattice in ["rhl1", "rhl"]:
        return RHL(pi, 60)
    elif lattice == "rhl2":
        return RHL(pi, 110)
    elif lattice == "mcl":
        return MCL(pi, 2 * pi, 3 * pi, alpha=80)
    elif lattice in ["mclc1", "mclc"]:
        return MCLC(pi, 1.5 * pi, 2 * pi, 80)
    elif lattice == "mclc2":
        return MCLC(1.5 * pi * sin(80 * _toradians), 1.5 * pi, 2 * pi, 80)
    elif lattice == "mclc3":
        return MCLC(pi, pi / 2, pi, 80)
    elif lattice == "mclc4":
        b = pi
        c = 1.5 * pi
        alpha = 80 * _toradians
        a = sqrt(c * b**2 * sin(alpha) ** 2 / (c - b * cos(alpha)))
        return MCLC(a, b, c, 80)
    elif lattice == "mclc5":
        return MCLC(pi, pi, pi, 60)
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
