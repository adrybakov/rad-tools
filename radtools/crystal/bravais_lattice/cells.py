from math import cos, pi, sin, sqrt, tan

import numpy as np

from radtools.routines import (
    param_from_cell,
    toradians,
    compare_numerically,
    cell_from_param,
    reciprocal_cell,
)
from radtools.crystal.constants import TRANSFORM_TO_CONVENTIONAL, REL_TOL, ABS_TOL_ANGLE
from radtools.routines import volume

__all__ = ["fix_cell"]


# Primitive cell`s construction
def CUB_cell(a: float):
    r"""
    Construct cubic primitive cell.

    See :ref:`lattice-cub` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float or int
        Length of the all three lattice vectors of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])


def FCC_cell(a: float):
    r"""
    Construct face-centred cubic primitive cell.

    See :ref:`lattice-fcc` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the all three lattice vectors of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])


def BCC_cell(a: float):
    r"""
    Construct body-centred cubic primitive cell.

    See :ref:`lattice-bcc` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the all three lattice vectors of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(
        [[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]]
    )


def TET_cell(a: float, c: float):
    r"""
    Construct tetragonal primitive cell.

    See :ref:`lattice-tet` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the two equal lattice vectors of the conventional cell.
    c : float
        Length of the third lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array([[a, 0, 0], [0, a, 0], [0, 0, c]])


def BCT_cell(a: float, c: float):
    r"""
    Construct body-centred tetragonal primitive cell.

    See :ref:`lattice-bct` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the two equal lattice vectors of the conventional cell.
    c : float
        Length of the third lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(
        [[-a / 2, a / 2, c / 2], [a / 2, -a / 2, c / 2], [a / 2, a / 2, -c / 2]]
    )


def ORC_cell(a: float, b: float, c: float):
    r"""
    Construct orthorhombic primitive cell.

    See :ref:`lattice-orc` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    return np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])


def ORCF_cell(a: float, b: float, c: float):
    r"""
    Construct face-centred orthorhombic primitive cell.

    See :ref:`lattice-orcf` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    return np.array([[0, b / 2, c / 2], [a / 2, 0, c / 2], [a / 2, b / 2, 0]])


def ORCI_cell(a: float, b: float, c: float):
    r"""
    Construct body-centred orthorhombic primitive cell.

    See :ref:`lattice-orci` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    a, b, c = tuple(sorted([a, b, c]))

    return np.array(
        [[-a / 2, b / 2, c / 2], [a / 2, -b / 2, c / 2], [a / 2, b / 2, -c / 2]]
    )


def ORCC_cell(a: float, b: float, c: float):
    r"""
    Construct base-centred orthorhombic primitive cell.

    See :ref:`lattice-orcc` for the definition of primitive and conventional cells.

    Order: :math:`a < b < c`. Input is reordered if necessary.

    Parameters
    ----------
    a : float
        Length of the smallest lattice vector of the conventional cell.
    b : float
        Length of the medium lattice vector of the conventional cell.
    c : float
        Length of the largest lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    a, b = tuple(sorted([a, b]))

    return np.array([[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]])


def HEX_cell(a: float, c: float):
    r"""
    Construct hexagonal primitive cell.

    See :ref:`lattice-hex` for the definition of primitive and conventional cells.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cell.
    c : float
        Length of the lattice vector of the conventional cell.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(
        [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
    )


def RHL_cell(a: float, alpha: float):
    r"""
    Construct rhombohedral primitive cell.

    See :ref:`lattice-rhl` for the definition of primitive and conventional cells.

    Condition :math:`\alpha < 120^{\circ}` is assumed.

    Parameters
    ----------
    a : float
        Length of the lattice vector of the conventional cell.
    alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    if alpha >= 120:
        raise ValueError("alpha has to be < 120 degrees.")

    alpha *= toradians
    return np.array(
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


def MCL_cell(a: float, b: float, c: float, alpha: float):
    r"""
    Construct monoclinic primitive cell.

    See :ref:`lattice-mcl` for the definition of primitive and conventional cells.

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

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    b, c = tuple(sorted([b, c]))
    if alpha > 90:
        alpha = 180 - alpha

    alpha *= toradians
    return np.array([[a, 0, 0], [0, b, 0], [0, c * cos(alpha), c * sin(alpha)]])


def MCLC_cell(a: float, b: float, c: float, alpha: float):
    r"""
    Construct base-centred monoclinic primitive cell.

    See :ref:`lattice-mclc` for the definition of primitive and conventional cells.

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

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    b, c = tuple(sorted([b, c]))
    if alpha > 90:
        alpha = 180.0 - alpha

    alpha *= toradians
    return np.array(
        [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )


# TODO work out the order of the vectors
def TRI_cell(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
    reciprocal=False,
):
    r"""
    Construct triclinic primitive cell.

    See :ref:`lattice-tri` for the definition of primitive and conventional cells.

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

    Returns
    -------
    prim_cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    cell = cell_from_param(a, b, c, alpha, beta, gamma)
    if reciprocal:
        return reciprocal_cell(cell)
    return cell


# Cell fixers
def fix_cell(cell, correct_lattice_type, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine it
    if required to ensure the unique choice of lattice vectors.

    See :ref:`lattice` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.
    correct_lattice_type : str
        Correct lattice type.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    cell = np.array(cell)
    functions = {
        "CUB": CUB_fix_cell,
        "FCC": FCC_fix_cell,
        "BCC": BCC_fix_cell,
        "TET": TET_fix_cell,
        "BCT": BCT_fix_cell,
        "ORC": ORC_fix_cell,
        "ORCF": ORCF_fix_cell,
        "ORCI": ORCI_fix_cell,
        "ORCC": ORCC_fix_cell,
        "HEX": HEX_fix_cell,
        "RHL": RHL_fix_cell,
        "MCL": MCL_fix_cell,
        "MCLC": MCLC_fix_cell,
        "TRI": TRI_fix_cell,
    }

    return functions[correct_lattice_type](cell, eps_rel)


def CUB_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the CUB lattice conditions.

    See :ref:`lattice-cub` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.
        Ignored here, but preserved for the unification of input.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(cell)


def FCC_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the FCC lattice conditions.

    See :ref:`lattice-fcc` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.
        Ignored here, but preserved for the unification of input.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(cell)


def BCC_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the BCC lattice conditions.

    See :ref:`lattice-bcc` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.
        Ignored here, but preserved for the unification of input.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(cell)


def TET_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the TET lattice conditions.

    See :ref:`lattice-tet` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    cell = np.array(cell)

    eps = eps_rel * volume(cell) ** (1.0 / 3.0)

    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    if compare_numerically(a, "==", c, eps):
        cell = [cell[2], cell[0], cell[1]]
    elif compare_numerically(b, "==", c, eps):
        cell = [cell[1], cell[2], cell[0]]

    return cell


def BCT_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the BCT lattice conditions.

    See :ref:`lattice-bct` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """
    cell = np.array(cell)

    eps = eps_rel * volume(cell) ** (1.0 / 3.0)

    a, b, c, alpha, beta, gamma = param_from_cell(
        TRANSFORM_TO_CONVENTIONAL["BCT"] @ cell
    )

    if compare_numerically(a, "==", c, eps):
        cell = [cell[2], cell[0], cell[1]]
    elif compare_numerically(b, "==", c, eps):
        cell = [cell[1], cell[2], cell[0]]

    return cell


def ORC_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the ORC lattice conditions.

    See :ref:`lattice-orc` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """
    cell = np.array(cell)
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)

    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    if compare_numerically(a, ">", b, eps):
        # minus preserves right-hand order
        cell = [cell[1], cell[0], -cell[2]]
        a, b = b, a
    if compare_numerically(a, ">", c, eps):
        cell = [cell[2], cell[0], cell[1]]
    elif compare_numerically(b, ">", c, eps):
        # minus preserves right-hand order
        cell = [cell[0], cell[2], -cell[1]]

    return cell


def ORCF_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the ORCF lattice conditions.

    See :ref:`lattice-orcf` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """
    cell = np.array(cell)
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)
    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    if compare_numerically(a, "<", b, eps):
        # minus preserves right-hand order
        # abc - > bca
        cell = [cell[1], cell[2], cell[0]]
        if compare_numerically(b, "<", c, eps):
            # minus preserves right-hand order
            # bca -> -cba
            cell = [-cell[1], cell[0], cell[2]]
        elif compare_numerically(c, "<", a, eps):
            # minus preserves right-hand order
            # bca -> b-ac
            cell = [cell[0], -cell[2], cell[1]]
    elif compare_numerically(a, "<", c, eps):
        # abc -> cab
        cell = [cell[2], cell[0], cell[1]]
    elif compare_numerically(b, "<", c, eps):
        # minus preserves right-hand order
        # abc -> a-cb
        cell = [cell[0], -cell[2], cell[1]]

    return cell


def ORCI_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the ORCI lattice conditions.

    See :ref:`lattice-orci` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitiv4e unit cell.
    """

    cell = np.array(cell)
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)
    a, b, c, alpha, beta, gamma = param_from_cell(
        TRANSFORM_TO_CONVENTIONAL["ORCI"] @ cell
    )

    if compare_numerically(a, ">", b, eps):
        # minus preserves right-hand order
        # abc - > bca
        cell = [cell[1], cell[2], cell[0]]
        if compare_numerically(b, ">", c, eps):
            # minus preserves right-hand order
            # bca -> -cba
            cell = [-cell[1], cell[0], cell[2]]
        elif compare_numerically(c, ">", a, eps):
            # minus preserves right-hand order
            # bca -> b-ac
            cell = [cell[0], -cell[2], cell[1]]
    elif compare_numerically(a, ">", c, eps):
        # abc -> cab
        cell = [cell[2], cell[0], cell[1]]
    elif compare_numerically(b, ">", c, eps):
        # minus preserves right-hand order
        # abc -> a-cb
        cell = [cell[0], -cell[2], cell[1]]

    return cell


def ORCC_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the ORCC lattice conditions.

    See :ref:`lattice-orcc` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    a, b, c, alpha, beta, gamma = param_from_cell(cell)
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)

    # a == c
    if compare_numerically(a, "==", c, eps):
        cell = [cell[2], cell[0], cell[1]]
    # b = c
    elif compare_numerically(b, "==", c, eps):
        cell = [cell[1], cell[2], cell[0]]

    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    cell = np.array(cell)

    a, b, c, alpha, beta, gamma = param_from_cell(
        TRANSFORM_TO_CONVENTIONAL["ORCC"] @ cell
    )

    if compare_numerically(a, ">", b, eps):
        # minus preserves right-hand order
        cell = [-cell[1], cell[0], cell[2]]

    return cell


def HEX_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the HEX lattice conditions.

    See :ref:`lattice-hex` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.
        Ignored here, but preserved for the unification of input.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """
    cell = np.array(cell)
    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    # a == c
    if compare_numerically(beta, "==", 120.0, ABS_TOL_ANGLE):
        cell = [cell[2], cell[0], cell[1]]
    # b = c
    elif compare_numerically(alpha, "==", 120.0, ABS_TOL_ANGLE):
        cell = [cell[1], cell[2], cell[0]]

    return cell


def RHL_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the RHL lattice conditions.

    See :ref:`lattice-rhl` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Relative tolerance for numerical comparison.
        Ignored here, but preserved for the unification of input.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    return np.array(cell)


def MCL_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the MCL lattice conditions.

    See :ref:`lattice-mcl` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    cell = np.array(cell)
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)
    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    # beta != 90
    if compare_numerically(beta, "!=", 90.0, ABS_TOL_ANGLE):
        cell = [cell[1], cell[2], cell[0]]
    # gamma != 90
    elif compare_numerically(gamma, "!=", 90.0, ABS_TOL_ANGLE):
        cell = [cell[2], cell[0], cell[1]]

    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    # alpha > 90 (cos(alpha) < 0)
    if compare_numerically(alpha, ">", 90.0, ABS_TOL_ANGLE):
        cell = [cell[0], cell[2], -cell[1]]

    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    # b > c
    if compare_numerically(b, ">", c, eps):
        cell = [-cell[0], cell[2], cell[1]]

    return cell


def MCLC_fix_cell(cell, eps_rel=REL_TOL):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the MCLC lattice conditions.

    See :ref:`lattice-mclc` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """

    eps = eps_rel * volume(cell) ** (1.0 / 3.0)
    a, b, c, alpha, beta, gamma = param_from_cell(cell)

    # a == c
    if compare_numerically(a, "==", c, eps):
        cell = [cell[2], cell[0], cell[1]]
    # b == c
    elif compare_numerically(b, "==", c, eps):
        cell = [cell[1], cell[2], cell[0]]

    cell = np.array(cell)
    a, b, c, alpha, beta, gamma = param_from_cell(
        TRANSFORM_TO_CONVENTIONAL["MCLC"] @ cell
    )
    alpha *= toradians

    # alpha > 90 (cos(alpha) < 0)
    if compare_numerically(cos(alpha), "<", 0, eps):
        cell = [cell[0], cell[2], -cell[1]]

    cell = np.array(cell)
    a, b, c, alpha, beta, gamma = param_from_cell(
        TRANSFORM_TO_CONVENTIONAL["MCLC"] @ cell
    )

    # b > c
    if compare_numerically(b, ">", c, eps):
        cell = [-cell[0], cell[2], cell[1]]

    return cell


# TODO
def TRI_fix_cell(cell, eps_rel=REL_TOL, resiprocal=False):
    r"""
    Analyse arbitrary cell and redefine vectors if required to satisfy the TRI lattice conditions.

    See :ref:`lattice-tri` for the details.

    Parameters
    ----------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    eps_rel : float, default ``REL_TOL``
        Tolerance for numerical comparison.
    resiprocal : bool, default False
        Whether to interpret input as reciprocal cell.

    Returns
    -------
    cell : (3,3) :numpy:`ndarray`
        Primitive unit cell.
    """
    eps = eps_rel * volume(cell) ** (1.0 / 3.0)

    return np.array(cell)
