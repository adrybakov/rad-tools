r"""14 Bravais lattice"""

from math import sqrt, sin, cos, tan, pi
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
from rad_tools.crystal.lattice import Lattice


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
    "cub",
    "fcc",
    "bcc",
    "tet",
    "bct1",
    "bct2",
    "orc",
    "orcf1",
    "orcf2",
    "orcf3",
    "orci",
    "orcc",
    "hex",
    "rhl1",
    "rhl2",
    "mcl",
    "mclc1",
    "mclc2",
    "mclc3",
    "mclc4",
    "mclc5",
    "tri1a",
    "tri1b",
    "tri2a",
    "tri2b",
]

DISTANCE_TOLERANCE = 10e-8
ANGLE_TOLERANCE = 10e-5


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
        Length of the lattice vectors.

    Attributes
    ----------
    a : float
        Length of the lattice vectors.
    """

    _pearson_symbol = "cP"

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive_cell = self.cell
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
        Length of the lattice vector.

    Attributes
    ----------
    a : float
        Length of the lattice vector
    """

    _pearson_symbol = "cF"

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive_cell = np.array([[0, a, a], [a, 0, a], [a, a, 0]]) / 2
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
        Length of the lattice vector.

    Attributes
    ----------
    a : float
        Length of the lattice vector
    """

    _pearson_symbol = "cI"

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive_cell = np.array(
            [
                [-a / 2, a / 2, a / 2],
                [a / 2, -a / 2, a / 2],
                [a / 2, a / 2, -a / 2],
            ]
        )
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
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    c : float
        Length of the lattice vector.
    """

    _pearson_symbol = "tP"

    def __init__(self, a: float, c: float) -> None:
        self.a = a
        self.c = c
        self.cell = np.diag([a, a, c])
        self.primitive_cell = self.cell
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
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    c : float
        Length of the lattice vector.
    """

    _pearson_symbol = "tI"

    def __init__(self, a: float, c: float) -> None:
        self.PLOT_NAMES["S"] = "$\\Sigma$"
        self.PLOT_NAMES["S1"] = "$\\Sigma_1$"
        if a == c:
            raise ValueError("Are you trying to create BCC Lattice (a == c)?")
        self.a = a
        self.c = c
        self.cell = np.diag([a, a, c])
        self.primitive_cell = np.array(
            [
                [-a / 2, a / 2, c / 2],
                [a / 2, -a / 2, c / 2],
                [a / 2, a / 2, -c / 2],
            ]
        )
        if self.variant == "BCT1":
            eta = (1 + self.c**2 / self.a**2) / 4
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

        elif self.variant == "BCT2":
            eta = (1 + self.a**2 / self.c**2) / 4
            zeta = self.a**2 / (2 * self.c**2)
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
    def variant(self):
        r"""
        Two variants of the Lattice.

        :math:`\text{BCT}_1: c < a` and :math:`\text{BCT}_2: c > a`
        """
        if self.a > self.c:
            return "BCT1"
        elif self.a < self.c:
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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
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

        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive_cell = self.cell
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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
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
        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive_cell = np.array(
            [
                [0, b / 2, c / 2],
                [a / 2, 0, c / 2],
                [a / 2, b / 2, 0],
            ]
        )
        if self.variant == "ORCF1":
            eta = (1 + self.a**2 / self.b**2 + self.a**2 / self.c**2) / 4
            zeta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4
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
        elif self.variant == "ORCF2":
            eta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4
            delta = (1 + self.b**2 / self.a**2 - self.b**2 / self.c**2) / 4
            phi = (1 + self.c**2 / self.b**2 - self.c**2 / self.a**2) / 4

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
        elif self.variant == "ORCF3":
            eta = (1 + self.a**2 / self.b**2 + self.a**2 / self.c**2) / 4
            zeta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4

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
    def variant(self):
        r"""
        Three variants of the Lattice.

        :math:`\text{ORCF}_1: \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_2: \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_3: \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        """
        if 1 / self.a**2 > 1 / self.b**2 + 1 / self.c**2:
            return "ORCF1"
        elif 1 / self.a**2 < 1 / self.b**2 + 1 / self.c**2:
            return "ORCF2"
        else:
            return "ORCF3"


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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
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
        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive_cell = np.array(
            [
                [-a / 2, b / 2, c / 2],
                [a / 2, -b / 2, c / 2],
                [a / 2, b / 2, -c / 2],
            ]
        )
        zeta = (1 + self.a**2 / self.c**2) / 4
        eta = (1 + self.b**2 / self.c**2) / 4
        delta = (self.b**2 - a**2) / (4 * self.c**2)
        mu = (self.a**2 + self.b**2) / (4 * self.c**2)

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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    """

    _pearson_symbol = "oS"

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct TET Lattice (a = b = c)?")
        self.a = a
        self.b = c
        self.c = b
        self.cell = np.diag([a, b, c])
        self.primitive_cell = np.array(
            [
                [a / 2, -b / 2, 0],
                [a / 2, b / 2, 0],
                [0, 0, c],
            ]
        )
        zeta = (1 + self.a**2 / self.b**2) / 4

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
        Length of the lattice vector.
    c : float
        Length of the lattice vector.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    c : float
        Length of the lattice vector.
    """

    _pearson_symbol = "hP"

    def __init__(self, a: float, c: float) -> None:
        self.a = a
        self.c = c
        self.cell = np.array(
            [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
        )
        self.primitive_cell = self.cell

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
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    alpha : float
        Angle between b and c. In degrees.
    """

    _pearson_symbol = "hR"

    def __init__(self, a: float, alpha: float) -> None:
        if alpha == 90:
            raise ValueError("Are you trying to construct CUB Lattice (alpha == 90)?")
        if alpha >= 120:
            raise ValueError("alpha has to be < 120 degrees.")
        self.a = a
        self.alpha = alpha
        self.cell = np.array(
            [
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
            ]
        )
        self.primitive_cell = self.cell
        if self.variant == "RHL1":
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
        elif self.variant == "RHL2":
            eta = 1 / (2 * tan(self.alpha / 180 * pi / 2) ** 2)
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
    def variant(self):
        r"""
        Two variants of the Lattice.

        :math:`\text{RHL}_1 \alpha < 90^{\circ}`,
        :math:`\text{RHL}_2 \alpha > 90^{\circ}`
        """
        if self.alpha < 90:
            return "RHL1"
        elif self.alpha > 90:
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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.
    """

    _pearson_symbol = "mP"

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.cell = np.array(
            [
                [a, 0, 0],
                [0, b, 0],
                [0, c * cos(alpha / 180 * pi), c * sin(alpha / 180 * pi)],
            ]
        )
        self.primitive_cell = self.cell

        eta = (1 - self.b * cos(self.alpha / 180 * pi) / self.c) / (
            2 * sin(self.alpha / 180 * pi) ** 2
        )
        nu = 1 / 2 - eta * self.c * cos(self.alpha / 180 * pi) / self.b

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
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.
    """

    _pearson_symbol = "mS"

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        if a > c:
            a, c = c, a
        if b > c:
            b, c = c, b
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.cell = np.array(
            [
                [a, 0, 0],
                [0, b, 0],
                [
                    0,
                    c * cos(alpha / 180 * pi),
                    c * sin(alpha / 180 * pi),
                ],
            ]
        )
        self.primitive_cell = np.array(
            [
                [a / 2, b / 2, 0],
                [-a / 2, b / 2, 0],
                [
                    0,
                    c * cos(alpha / 180 * pi),
                    c * sin(alpha / 180 * pi),
                ],
            ]
        )
        # Parameters
        if self.variant in ["MCLC1", "MCLC2"]:
            zeta = (2 - b * cos(alpha / 180 * pi) / c) / (
                4 * sin(alpha / 180 * pi) ** 2
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha / 180 * pi) ** 2)
            phi = psi + (3 / 4 - psi) * b * cos(alpha / 180 * pi) / c
        elif self.variant in ["MCLC3", "MCLC4"]:
            mu = (1 + b**2 / a**2) / 4
            delta = b * c * cos(alpha / 180 * pi) / (2 * a**2)
            zeta = (
                mu
                - 1 / 4
                + (1 - b * cos(alpha / 180 * pi) / c) / (4 * sin(alpha / 180 * pi) ** 2)
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            phi = 1 + zeta - 2 * mu
            psi = eta - 2 * delta
        elif self.variant == "MCLC5":
            zeta = (
                b**2 / a**2
                + (1 - b * cos(alpha / 180 * pi) / c) / sin(alpha / 180 * pi) ** 2
            ) / 4
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            mu = (
                eta / 2
                + b**2 / (4 * a**2)
                - b * c * cos(alpha / 180 * pi) / (2 * a**2)
            )
            nu = 2 * mu - zeta
            rho = 1 - zeta * a**2 / b**2
            omega = (
                (4 * nu - 1 - b**2 * sin(alpha / 180 * pi) ** 2 / a**2)
                * c
                / (2 * b * cos(alpha / 180 * pi))
            )
            delta = zeta * c * cos(alpha / 180 * pi) / b + omega / 2 - 1 / 4

        # Path
        if self.variant == "MCLC1":
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
        elif self.variant == "MCLC2":
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
        elif self.variant == "MCLC3":
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
        elif self.variant == "MCLC4":
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
        elif self.variant == "MCLC5":
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
    def variant(self):
        r"""
        Five variant of the Lattice.

        :math:`\text{MCLC}_1: k_{\gamma} > 90^{\circ}`,
        :math:`\text{MCLC}_2: k_{\gamma} = 90^{\circ}`,
        :math:`\text{MCLC}_3: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`
        :math:`\text{MCLC}_4: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`
        :math:`\text{MCLC}_5: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`
        """

        if abs(self.k_gamma - 90) < 10e-5:
            return "MCLC2"
        elif self.k_gamma > 90:
            return "MCLC1"
        # TODO think about the ccriteria of accuracy
        elif self.k_gamma < 90:
            if (
                abs(
                    self.b * cos(self.alpha / 180 * pi) / self.c
                    + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                    - 1
                )
                < 10e-8
            ):
                return "MCLC4"
            elif (
                self.b * cos(self.alpha / 180 * pi) / self.c
                + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                < 1
            ):
                return "MCLC3"
            elif (
                self.b * cos(self.alpha / 180 * pi) / self.c
                + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                > 1
            ):
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

    Parameters
    ----------
    a : float
        Length of the lattice vector.
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.
    beta : float
        Angle between a and c. In degrees.
    gamma : float
        Angle between a and b. In degrees.


    Attributes
    ----------
    a : float
        Length of the lattice vector
    b : float
        Length of the lattice vector.
    c : float
        Length of the lattice vector.
    alpha : float
        Angle between b and c. In degrees.
    beta : float
        Angle between a and c. In degrees.
    gamma : float
        Angle between a and b. In degrees.
    """

    _pearson_symbol = "aP"

    def __init__(
        self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.cell = np.array(
            [
                [a, 0, 0],
                [b * cos(gamma / 180 * pi), b * sin(gamma / 180 * pi), 0],
                [
                    c * cos(beta / 180 * pi),
                    c
                    / sin(gamma / 180 * pi)
                    * (
                        cos(alpha / 180 * pi)
                        - cos(beta / 180 * pi) * cos(gamma / 180 * pi)
                    ),
                    c
                    / sin(gamma / 180 * pi)
                    * sqrt(
                        sin(gamma / 180 * pi) ** 2
                        - cos(alpha / 180 * pi) ** 2
                        - cos(beta / 180 * pi) ** 2
                        + 2
                        * cos(alpha / 180 * pi)
                        * cos(beta / 180 * pi)
                        * cos(gamma / 180 * pi)
                    ),
                ],
            ]
        )
        self.primitive_cell = self.cell
        if self.variant in ["TRI1a", "TRI1b"]:
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
        elif self.variant in ["TRI2a", "TRI2b"]:
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

    @property
    def variant(self):
        r"""
        Four variants of the Lattice.


        :math:`\text{TRI}_{1a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`
        :math:`\text{TRI}_{1b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`
        :math:`\text{TRI}_{2a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`
        :math:`\text{TRI}_{2b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`
        """
        if abs(self.k_gamma - 90) < 10e-5:
            if self.k_alpha > 90 and self.k_beta > 90:
                return "TRI2a"
            elif self.k_alpha < 90 and self.k_beta < 90:
                return "TRI2b"
        elif (min(self.k_gamma, self.k_beta, self.k_alpha)) > 90:
            return "TRI1a"
        elif (max(self.k_gamma, self.k_beta, self.k_alpha)) < 90:
            return "TRI1a"


# Examples
cub = CUB(pi)
fcc = FCC(pi)
bcc = BCC(pi)
tet = TET(pi, 2 * pi)
bct1 = BCT(2 * pi, pi)
bct2 = BCT(pi, 2 * pi)
orc = ORC(pi, 2 * pi, 3 * pi)
orcf1 = ORCF(0.9 * pi, 5 / 4 * pi, 5 / 3 * pi)
orcf2 = ORCF(1.1 * pi, 5 / 4 * pi, 5 / 3 * pi)
orcf3 = ORCF(pi, 5 / 4 * pi, 5 / 3 * pi)
orci = ORCI(pi, 2 * pi, 3 * pi)
orcc = ORCC(pi, 2 * pi, 3 * pi)
hex = HEX(pi, 2 * pi)
rhl1 = RHL(pi, 70)
rhl2 = RHL(pi, 110)
mcl = MCL(pi, 2 * pi, 3 * pi, alpha=80)
mclc1 = MCLC(1 * pi, 1.5 * pi, 2 * pi, 80)
mclc2 = MCLC(1.47721 * pi, 1.5 * pi, 2 * pi, 80)
mclc3 = MCLC(pi, pi / 2, pi, 80)
mclc4 = MCLC(1.06486355 * pi, pi, 1.2 * pi, 80)
mclc5 = MCLC(pi, pi, pi, 60)
tri1a = TRI(2 * pi, 3 * pi, 4 * pi, 60, 70, 80)
tri1b = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
tri2a = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
tri2b = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
examples = [
    cub,
    fcc,
    bcc,
    tet,
    bct1,
    bct2,
    orc,
    orcf1,
    orcf2,
    orcf3,
    orci,
    orcc,
    hex,
    rhl1,
    rhl2,
    mcl,
    mclc1,
    mclc2,
    mclc3,
    mclc4,
    mclc5,
    tri1a,
    tri1b,
    tri2a,
    tri2b,
]

if __name__ == "__main__":
    from math import pi

    print(
        f"BCT1 {bct1.variant}",
        f"BCT2 {bct2.variant}",
        f"ORCF1 {orcf1.variant}",
        f"ORCF2 {orcf2.variant}",
        f"ORCF3 {orcf3.variant}",
        f"RHL1 {rhl1.variant}",
        f"RHL2 {rhl2.variant}",
        f"MCLC1 {mclc1.variant}",
        f"MCLC2 {mclc2.variant}",
        f"MCLC3 {mclc3.variant}",
        f"MCLC4 {mclc4.variant}",
        f"MCLC5 {mclc5.variant}",
        f"TRI1a {tri1a.variant}",
        f"TRI1b {tri1b.variant}",
        f"TRI2a {tri2a.variant}",
        f"TRI2b {tri2b.variant}",
        sep="\n",
    )

    for e in examples:
        l = tri1a
        print(l.variant)
        l.prepare_figure()
        l.plot("brillouin_kpath")
        l.show()
        break


# TODO FIX TRI Lattice
