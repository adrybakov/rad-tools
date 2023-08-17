from math import cos, sin, sqrt

import numpy as np

import radtools.crystal.cell as Cell
from radtools.constants import TORADIANS
from radtools.crystal.lattice import Lattice

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
]


# Primitive cell`s construction
def CUB(a: float, return_cell=False):
    r"""
    Construct cubic primitive lattice.

    See :ref:`guide_cub` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float or int
        Length of the all three lattice vectors of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Cubic lattice or cell.
    """

    cell = np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
    if return_cell:
        return cell
    return Lattice(cell)


def FCC(a: float, return_cell=False):
    r"""
    Construct face-centred cubic primitive lattice.

    See :ref:`guide_fcc` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the all three lattice vectors of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Face-centred cubic lattice or cell.
    """

    cell = np.array([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])
    if return_cell:
        return cell
    return Lattice(cell)


def BCC(a: float, return_cell=False):
    r"""
    Construct body-centred cubic primitive lattice.

    See :ref:`guide_bcc` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the all three lattice vectors of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Body-centred cubic lattice or cell.
    """

    cell = np.array(
        [[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def TET(a: float, c: float, return_cell=False):
    r"""
    Construct tetragonal primitive lattice.

    See :ref:`guide_tet` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the two equal lattice vectors of the conventional cell.
    c : float
        Length of the third lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Tetragonal lattice or cell.
    """

    cell = np.array([[a, 0, 0], [0, a, 0], [0, 0, c]])
    if return_cell:
        return cell
    return Lattice(cell)


def BCT(a: float, c: float, return_cell=False):
    r"""
    Construct body-centred tetragonal primitive lattice.

    See :ref:`guide_bct` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the two equal lattice vectors of the conventional cell.
    c : float
        Length of the third lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Body-centred tetragonal lattice cell.
    """

    cell = np.array(
        [[-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def ORC(a: float, b: float, c: float, return_cell=False):
    r"""
    Construct orthorhombic primitive lattice.

    See :ref:`guide_orc` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Orthorhombic lattice or cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    cell = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    if return_cell:
        return cell
    return Lattice(cell)


def ORCF(a: float, b: float, c: float, return_cell=False):
    r"""
    Construct face-centred orthorhombic primitive lattice.

    See :ref:`guide_orcf` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Face-centred orthorhombic lattice or cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    cell = np.array([[0, b / 2, c / 2], [a / 2, 0, c / 2], [a / 2, b / 2, 0]])
    if return_cell:
        return cell
    return Lattice(cell)


def ORCI(a: float, b: float, c: float, return_cell=False):
    r"""
    Construct body-centred orthorhombic primitive lattice.

    See :ref:`guide_orci` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Body-centred orthorhombic lattice or cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    cell = np.array(
        [[-a / 2, b / 2, c / 2], [a / 2, -b / 2, c / 2], [a / 2, b / 2, -c / 2]]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def ORCC(a: float, b: float, c: float, return_cell=False):
    r"""
    Construct base-centred orthorhombic primitive lattice.

    See :ref:`guide_orcc` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Base-centred orthorhombic lattice or cell.
    """

    a, b = tuple(sorted([a, b]))

    cell = np.array([[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]])
    if return_cell:
        return cell
    return Lattice(cell)


def HEX(a: float, c: float, return_cell=False):
    r"""
    Construct hexagonal primitive lattice.

    See :ref:`guide_hex` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cell.
    c : float
        Length of the lattice vector of the conventional cell.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Hexagonal lattice or cell.
    """

    cell = np.array(
        [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def RHL(a: float, alpha: float, return_cell=False):
    r"""
    Construct rhombohedral primitive lattice.

    See :ref:`guide_rhl` for the definition of primitive and conventional cells.

    Condition :math:`\alpha < 120^{\circ}` is assumed.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cell.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Rhombohedral lattice or cell.
    """

    if alpha >= 120:
        raise ValueError("alpha has to be < 120 degrees.")

    alpha *= TORADIANS
    cell = np.array(
        [
            [a * cos(alpha / 2), -a * sin(alpha / 2), 0],
            [a * cos(alpha / 2), a * sin(alpha / 2), 0],
            [
                a * cos(alpha) / cos(alpha / 2),
                0,
                a * sqrt(1 - cos(alpha) ** 2 / cos(alpha / 2) ** 2),
            ],
        ]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def MCL(a: float, b: float, c: float, alpha: float, return_cell=False):
    r"""
    Construct monoclinic primitive lattice.

    See :ref:`guide_mcl` for the definition of primitive and conventional cells.

    Order: :math:`b \le c`, :math:`\alpha < 90^{\circ}`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the first lattice vector of the conventional cell. (The one oriented along x axis)
    b : float
        Length of the shorter of the two remaining lattice vectors of the conventional cell.
    c : float
        Length of the longer of the two remaining lattice vectors of the conventional cell.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Monoclinic lattice or cell.
    """

    b, c = tuple(sorted([b, c]))
    if alpha > 90:
        alpha = 180 - alpha

    alpha *= TORADIANS
    cell = np.array([[a, 0, 0], [0, b, 0], [0, c * cos(alpha), c * sin(alpha)]])
    if return_cell:
        return cell
    return Lattice(cell)


def MCLC(a: float, b: float, c: float, alpha: float, return_cell=False):
    r"""
    Construct base-centred monoclinic primitive lattice.

    See :ref:`guide_mclc` for the definition of primitive and conventional cells.

    Order: :math:`b \le c`, :math:`\alpha < 90^{\circ}`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the first lattice vector of the conventional cell. (The one oriented along x axis)
    b : float
        Length of the shorter of the two remaining lattice vectors of the conventional cell.
    c : float
        Length of the longer of the two remaining lattice vectors of the conventional cell.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Base-centred monoclinic lattice or cell.
    """

    b, c = tuple(sorted([b, c]))
    if alpha > 90:
        alpha = 180.0 - alpha

    alpha *= TORADIANS
    cell = np.array(
        [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )
    if return_cell:
        return cell
    return Lattice(cell)


def TRI(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    reciprocal=False,
    return_cell=False,
):
    r"""
    Construct triclinic primitive lattice.

    See :ref:`guide_tri` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cell.
    b : float
        Length of the lattice vector of the conventional cell.
    c : float
        Length of the lattice vector of the conventional cell.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.
    beta : float
        Angle between vectors :math:`a_1` and :math:`a_3` of the conventional cell. In degrees.
    gamma : float
        Angle between vectors :math:`a_1` and :math:`a_2` of the conventional cell. In degrees.
    reciprocal : bool, default False
        Whether to interpret input as reciprocal parameters.
    return_cell : bool, default False
        Whether to return the cell instead of the :py:class:`.Lattice` object.

    Returns
    -------
    lattice : :py:class:`.Lattice` or (3, 3) :numpy:`ndarray`
        Triclinic lattice or cell.
    """

    cell = Cell.from_params(a, b, c, alpha, beta, gamma)
    if reciprocal:
        cell = Cell.reciprocal(cell)
    if return_cell:
        return cell
    return Lattice(cell)
