from math import cos, pi, sin

from radtools.crystal.identify import lepage
from radtools.crystal.lattice import Lattice
from radtools.routines import param_from_cell, toradians
from radtools.crystal.constants import EPS_REL, BRAVAIS_LATTICE_VARIATIONS

from radtools.crystal.bravais_lattice.cells import (
    CUB_cell,
    FCC_cell,
    BCC_cell,
    TET_cell,
    BCT_cell,
    ORC_cell,
    ORCF_cell,
    ORCI_cell,
    ORCC_cell,
    HEX_cell,
    RHL_cell,
    MCL_cell,
    MCLC_cell,
    TRI_cell,
    CUB_fix_cell,
    FCC_fix_cell,
    BCC_fix_cell,
    TET_fix_cell,
    BCT_fix_cell,
    ORC_fix_cell,
    ORCF_fix_cell,
    ORCI_fix_cell,
    ORCC_fix_cell,
    HEX_fix_cell,
    RHL_fix_cell,
    MCL_fix_cell,
    MCLC_fix_cell,
    TRI_fix_cell,
)

__all__ = [
    "bravais_lattice_from_param",
    "bravais_lattice_from_cell",
    "lattice_example",
]


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
        return Lattice(CUB_cell(a))
    if lattice_type == "FCC":
        a = (a + b + c) / 3
        return Lattice(FCC_cell(a))
    if lattice_type == "BCC":
        a = (a + b + c) / 3
        return Lattice(BCC_cell(a))
    if lattice_type == "TET":
        if a == b:
            return Lattice(TET_cell((a + b) / 2, c))
        elif a == c:
            return Lattice(TET_cell((a + c) / 2, b))
        elif b == c:
            return Lattice(TET_cell((b + c) / 2, a))
    if lattice_type == "BCT":
        if a == b:
            return Lattice(BCT_cell((a + b) / 2, c))
        elif a == c:
            return Lattice(BCT_cell((a + c) / 2, b))
        elif b == c:
            return Lattice(BCT_cell((b + c) / 2, a))
    if lattice_type == "ORC":
        return Lattice(ORC_cell(a, b, c))
    if lattice_type == "ORCF":
        return Lattice(ORCF_cell(a, b, c))
    if lattice_type == "ORCC":
        return Lattice(ORCC_cell(a, b, c))
    if lattice_type == "ORCI":
        return Lattice(ORCI_cell(a, b, c))
    if lattice_type == "HEX":
        if abs(a - b) < abs(a - c) and abs(a - b) < abs(b - c):
            return Lattice(HEX_cell((a + b) / 2, c))
        elif abs(a - c) < abs(a - b) and abs(a - c) < abs(b - c):
            return Lattice(HEX_cell((a + c) / 2, b))
        elif abs(b - c) < abs(a - c) and abs(b - c) < abs(a - b):
            return Lattice(HEX_cell((b + c) / 2, a))
    if lattice_type == "RHL":
        return Lattice(RHL_cell((a + b + c) / 3, (alpha + beta + gamma) / 2))
    if lattice_type == "MCL":
        alpha, beta, gamma = [alpha, beta, gamma].sort()
        return Lattice(MCL_cell(a, b, c, alpha))
    if lattice_type == "MCLC":
        alpha, beta, gamma = [alpha, beta, gamma].sort()
        return Lattice(MCLC_cell(a, b, c, alpha))
    return Lattice(TRI_cell(a, b, c, alpha, beta, gamma))


def bravais_lattice_from_cell(cell, eps_rel=EPS_REL) -> Lattice:
    r"""
    Create Bravais lattice from cell matrix.

    Orientation of the crystal is respected, however the lattice vectors are renamed
    with respect to [1]_.

    This function may change the angles and the lengths of the cell.
    It preserve the volume, right- or left- handedness, lattice type and variation
    of the cell.

    The lattice vector`s lengths are preserved as a set.

    The angles between the lattice vectors are preserved as a set with possible
    changes of the form: :math:`angle \rightarrow 180 - angle`.

    The returned cell may not be the same as the input one, but it is translationally
    equivalent.

    Parameters
    ----------
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]

    eps_rel : float, default 1e-5
        Relative tolerance for determining the lattice type.

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

    lattice_type = lepage(*param_from_cell(cell), eps_rel=eps_rel)

    if lattice_type == "CUB":
        return Lattice(CUB_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "FCC":
        return Lattice(FCC_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "BCC":
        return Lattice(BCC_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "TET":
        return Lattice(TET_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "BCT":
        return Lattice(BCT_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "ORC":
        return Lattice(ORC_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "ORCF":
        return Lattice(ORCF_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "ORCC":
        return Lattice(ORCC_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "ORCI":
        return Lattice(ORCI_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "HEX":
        return Lattice(HEX_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "RHL":
        return Lattice(RHL_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "MCL":
        return Lattice(MCL_fix_cell(cell=cell, eps_rel=eps_rel))
    if lattice_type == "MCLC":
        return Lattice(MCLC_fix_cell(cell=cell, eps_rel=eps_rel))
    return Lattice(TRI_fix_cell(cell=cell, eps_rel=eps_rel))


def lattice_example(
    lattice_name: str = None,
):
    r"""
    Return an example of the lattice.

    Parameters
    ----------
    lattice_name : str, optional
        Name of the lattice to be returned.
        For available names see documentation of each Bravais lattice class.
        Lowercased before usage.

    Returns
    -------
    lattice : Lattice or list
        :py:class:`.Lattice` class is returned.
        If no math found a list with available examples is returned.
    """

    if (
        not isinstance(lattice_name, str)
        or lattice_name not in BRAVAIS_LATTICE_VARIATIONS
    ):
        message = (
            f"There is no {lattice_name} Bravais lattice. Available examples are:\n"
        )
        for name in BRAVAIS_LATTICE_VARIATIONS:
            message += f"  * {name}\n"
        raise ValueError(message)

    lattice_name = lattice_name.lower()

    if lattice_name == "cub":
        return Lattice(CUB_cell(pi))
    elif lattice_name == "fcc":
        return Lattice(FCC_cell(pi))
    elif lattice_name == "bcc":
        return Lattice(BCC_cell(pi))
    elif lattice_name == "tet":
        return Lattice(TET_cell(pi, 1.5 * pi))
    elif lattice_name in ["bct1", "bct"]:
        return Lattice(BCT_cell(1.5 * pi, pi))
    elif lattice_name == "bct2":
        return Lattice(BCT_cell(pi, 1.5 * pi))
    elif lattice_name == "orc":
        return Lattice(ORC_cell(pi, 1.5 * pi, 2 * pi))
    elif lattice_name in ["orcf1", "orcf"]:
        return Lattice(ORCF_cell(0.7 * pi, 5 / 4 * pi, 5 / 3 * pi))
    elif lattice_name == "orcf2":
        return Lattice(ORCF_cell(1.2 * pi, 5 / 4 * pi, 5 / 3 * pi))
    elif lattice_name == "orcf3":
        return Lattice(ORCF_cell(pi, 5 / 4 * pi, 5 / 3 * pi))
    elif lattice_name == "orci":
        return Lattice(ORCI_cell(pi, 1.3 * pi, 1.7 * pi))
    elif lattice_name == "orcc":
        return Lattice(ORCC_cell(pi, 1.3 * pi, 1.7 * pi))
    elif lattice_name == "hex":
        return Lattice(HEX_cell(pi, 2 * pi))
    elif lattice_name in ["rhl1", "rhl"]:
        # If alpha = 60 it is effectively FCC!
        return Lattice(RHL_cell(pi, 70))
    elif lattice_name == "rhl2":
        return Lattice(RHL_cell(pi, 110))
    elif lattice_name == "mcl":
        return Lattice(MCL_cell(pi, 1.3 * pi, 1.6 * pi, alpha=75))
    elif lattice_name in ["mclc1", "mclc"]:
        return Lattice(MCLC_cell(pi, 1.4 * pi, 1.7 * pi, 80))
    elif lattice_name == "mclc2":
        return Lattice(
            MCLC_cell(1.4 * pi * sin(75 * toradians), 1.4 * pi, 1.7 * pi, 75)
        )
    elif lattice_name == "mclc3":
        b = pi
        x = 1.1
        alpha = 78
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha) * 1.8
        a = x * b * sin(ralpha)
        return Lattice(MCLC_cell(a, b, c, alpha))
    elif lattice_name == "mclc4":
        b = pi
        x = 1.2
        alpha = 65
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha)
        a = x * b * sin(ralpha)
        return Lattice(MCLC_cell(a, b, c, alpha))
    elif lattice_name == "mclc5":
        b = pi
        x = 1.4
        alpha = 53
        ralpha = alpha * toradians
        c = b * (x**2) / (x**2 - 1) * cos(ralpha) * 0.9
        a = x * b * sin(ralpha)
        return Lattice(MCLC_cell(a, b, c, alpha))
    elif lattice_name in ["tri1a", "tri1", "tri", "tria"]:
        return Lattice(TRI_cell(1, 1.5, 2, 120, 110, 100, reciprocal=True))
    elif lattice_name in ["tri2a", "tri2"]:
        return Lattice(TRI_cell(1, 1.5, 2, 120, 110, 90, reciprocal=True))
    elif lattice_name in ["tri1b", "trib"]:
        return Lattice(TRI_cell(1, 1.5, 2, 60, 70, 80, reciprocal=True))
    elif lattice_name == "tri2b":
        return Lattice(TRI_cell(1, 1.5, 2, 60, 70, 90, reciprocal=True))
