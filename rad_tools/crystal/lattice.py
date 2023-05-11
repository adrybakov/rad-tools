r"""
Lattice
"""
from math import pi
import numpy as np


class Lattice:
    def __init__(self, a, b, c) -> None:
        self.cell = np.array([a, b, c])

    @property
    def a(self):
        return self._cell[0]

    @property
    def b(self):
        return self._cell[1]

    @property
    def c(self):
        return self._cell[2]

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        .. math::

            V = \delta_{\vec{A}}\cdot(\delta_{\vec{B}}\times\delta_{\vec{C}})
        """

        return np.dot(self.a, np.cross(self.b, self.c))

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector.

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{b}\times\vec{c}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.b, self.c)

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector.

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{c}\times\vec{a}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.c, self.a)

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector.

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}\times\vec{b}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.a, self.b)

    def get_bravais_lattice(self):
        r"""
        Return Bravais lattice.

        """
        return
        pass
