r"""
Lattice
"""
from math import pi
import numpy as np


class Lattice:
    def __init__(self, a, b, c) -> None:
        self.cell = np.array([a, b, c])

    @property
    def reciprocal_cell(self):
        return np.array(
            [
                2 * pi / self.unit_cell_volume * np.cross(self.a2, self.a3),
                2 * pi / self.unit_cell_volume * np.cross(self.a3, self.a1),
                2 * pi / self.unit_cell_volume * np.cross(self.a1, self.a2),
            ]
        )

    @property
    def a1(self):
        return self.cell[0]

    @property
    def a2(self):
        return self.cell[1]

    @property
    def a3(self):
        return self.cell[2]

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        .. math::

            V = \delta_{\vec{A}}\cdot(\delta_{\vec{B}}\times\delta_{\vec{C}})
        """

        return np.dot(self.a1, np.cross(self.a2, self.a3))

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector.

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{b}\times\vec{c}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[0]

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector.

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{c}\times\vec{a}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[1]

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector.

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}\times\vec{b}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[2]

    def get_bravais_lattice(self):
        r"""
        Return Bravais lattice.

        """
        return
        pass
