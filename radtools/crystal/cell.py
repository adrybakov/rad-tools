from math import cos, pi, sin, sqrt

import numpy as np

from radtools.constants import TORADIANS
from radtools.geometry import angle, parallelepiped_check, volume

__all__ = [
    "reciprocal",
    "from_params",
    "params",
    "primitive",
]


def reciprocal(cell):
    r"""
    Computes reciprocal cell.

    .. versionadded:: 0.7

    Parameters
    ----------
    cell : (3, 3) |array_like|_
        Cell matrix, rows are interpreted as vectors.

    Returns
    -------
    reciprocal_cell : (3, 3) :numpy:`ndarray`
        Reciprocal cell matrix, rows are interpreted as vectors.
        :math:`cell = (\vec{v}_1, \vec{v}_2, \vec{v}_3)`, where

        .. math::

            \begin{matrix}
                \vec{b}_1 = \dfrac{2\pi}{V}\vec{a}_2\times\vec{a}_3 \\
                \vec{b}_2 = \dfrac{2\pi}{V}\vec{a}_3\times\vec{a}_1 \\
                \vec{b}_3 = \dfrac{2\pi}{V}\vec{a}_1\times\vec{a}_2 \\
            \end{matrix}

    """

    vol = volume(cell)
    reciprocal_cell = np.array(
        [
            2 * pi / vol * np.cross(cell[1], cell[2]),
            2 * pi / vol * np.cross(cell[2], cell[0]),
            2 * pi / vol * np.cross(cell[0], cell[1]),
        ]
    )
    return reciprocal_cell


def from_params(a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0):
    r"""
    Return cell from lattice parameters.

    .. versionadded:: 0.7

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

    Returns
    -------
    cell : (3, 3) :numpy:`ndarray`
        Cell matrix.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]

    Raises
    ------
    ValueError
        If parameters could not form a parallelepiped.

    See Also
    --------
    parallelepiped_check : Check if parameters could form a parallelepiped.
    """
    parallelepiped_check(a, b, c, alpha, beta, gamma, raise_error=True)
    alpha = alpha * TORADIANS
    beta = beta * TORADIANS
    gamma = gamma * TORADIANS
    return np.array(
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
        ],
        dtype=float,
    )


def params(cell):
    r"""
    Return lattice parameters from cell.

    .. versionadded:: 0.7

    Parameters
    ----------
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]

    Returns
    -------
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
    """

    return (
        np.linalg.norm(cell[0]),
        np.linalg.norm(cell[1]),
        np.linalg.norm(cell[2]),
        angle(cell[1], cell[2]),
        angle(cell[0], cell[2]),
        angle(cell[0], cell[1]),
    )


# TODO
def primitive(cell, atoms):
    r"""
    Compute primitive cell.

    .. versionadded:: 0.8

    Parameters
    ----------
    cell : (3, 3) |array_like|_
        Cell matrix, rows are interpreted as vectors.
    atoms : list of :py:class:`.Atom`
        Atoms in the cell.
        ``position`` attribute of the atom is interpreted as relative
        position in the cell.

    Returns
    -------
    primitive_cell : (3, 3) :numpy:`ndarray`
        Primitive cell matrix, rows are interpreted as vectors.
    primitive_atoms : list of :py:class:`.Atom`
        Atoms in the primitive cell.
        ``position`` attribute of the atom is interpreted as relative
        position in the primitive cell.
    """
    raise NotImplementedError
