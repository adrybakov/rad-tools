r"""
Collection of small routines and constants, 
which may be used across the whole package.
"""

import sys
from math import asin, cos, pi, sin, sqrt

import numpy as np
from termcolor import cprint

__all__ = [
    "print_2D_array",
    "volume",
    "angle",
    "absolute_to_relative",
    "get_permutation",
    "cell_from_param",
    "reciprocal_cell",
]

RED = "#FF4D67"
GREEN = "#58EC2E"
ORANGE = "#F7CB3D"
BLUE = "#274DD1"
PURPLE = "#DC5CFF"

TOLERANCE = 1e-4
TOL_BASE = 4

_todegrees = 180 / pi
_toradians = pi / 180


def atom_mark_to_latex(mark):
    r"""
    Latexifier for atom marks.

    Cr12 -> Cr\ :sub:`12`\.

    Parameters
    ----------
    mark : str
        Mark of atom.

    Returns
    -------
    new_mark : str
        Latex version of the mark.
    """
    numbers = "0123456789"
    new_mark = "$"
    insert_underline = False
    for symbol in mark:
        if symbol in numbers and not insert_underline:
            insert_underline = True
            new_mark += "_{"
        new_mark += symbol
    new_mark += "}$"
    return new_mark


def rot_angle(x, y, dummy=False):
    r"""
    Rotational angle from 2D vector.

    Mathematically positive => counterclockwise.
    From [0 to 360)

    Parameters
    ----------
    x : float or int
        x coordinate of a vector.
    y : float or int
        y coordinate of a vector.
    """
    # rot_cos = x / (x ** 2 + y ** 2) ** 0.5
    # rot_angle = m.acos(rot_cos) / m.pi * 180
    try:
        sin = abs(y) / sqrt(x**2 + y**2)
    except ZeroDivisionError:
        raise ValueError("Angle is ill defined (x = y = 0).")
    if x > 0:
        if y > 0:
            return asin(sin) / pi * 180
        elif y == 0:
            return 0
        elif y < 0:
            if not dummy:
                return -asin(sin) / pi * 180
            return 360 - asin(sin) / pi * 180
    elif x == 0:
        if y > 0:
            return 90
        elif y == 0:
            raise ValueError("Angle is ill defined (x = y = 0).")
        elif y < 0:
            if not dummy:
                return 90
            return 270
    elif x < 0:
        if y > 0:
            if not dummy:
                return -asin(sin) / pi * 180
            return 180 - asin(sin) / pi * 180
        elif y == 0:
            if not dummy:
                return 0
            return 180
        elif y < 0:
            if not dummy:
                return asin(sin) / pi * 180
            return 180 + asin(sin) / pi * 180


def absolute_to_relative(cell, absolute):
    r"""
    Compute relative coordinates with respect to the unit cell.

    Parameters
    ----------
    cell : 3 x 3 array
        Lattice vectors.
    absolute : (3,) |array_like|_
        Absolute coordinates.

    Returns
    -------
    relative : 1 x 3 array
        Relative coordinates.
    """

    a = np.array(cell[0], dtype=float)
    b = np.array(cell[1], dtype=float)
    c = np.array(cell[2], dtype=float)
    v = np.array(absolute, dtype=float)
    if (v == np.zeros(3)).all():
        return np.zeros(3)
    B = np.array([np.dot(a, v), np.dot(b, v), np.dot(c, v)])
    A = np.array(
        [
            [np.dot(a, a), np.dot(a, b), np.dot(a, c)],
            [np.dot(b, a), np.dot(b, b), np.dot(b, c)],
            [np.dot(c, a), np.dot(c, b), np.dot(c, c)],
        ]
    )
    relative = np.linalg.solve(A, B)
    return relative


def volume(*args):
    r"""
    Computes volume.

    .. versionadded:: 0.7

    Three type of arguments are expected:

    * One argument.
        Matrix ``cell``.
        Volume is computed as:

        .. math::
            V = v_1 \cdot (v_2 \times v_3)
    * Three arguments.
        Vectors ``v1``, ``v2``, ``v3``.
        Volume is computed as:

        .. math::
            V = v_1 \cdot (v_2 \times v_3)
    * Six arguments.
        Parameters ``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma``.
        Volume is computed as:

        .. math::
            V = abc\sqrt(1+2\cos(\alpha)\cos(\beta)\cos(\gamma)-\cos^2(\alpha)-\cos^2(\beta)-\cos^2(\gamma))

    Parameters
    ----------
    v1 : (3,) |array_like|_
        First vector.
    v2 : (3,) |array_like|_
        Second vector.
    v3 : (3,) |array_like|_
        Third vector.
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.
    a : float, default 1
        Length of the :math:`v_1` vector.
    b : float, default 1
        Length of the :math:`v_2` vector.
    c : float, default 1
        Length of the :math:`v_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`v_2` and :math:`v_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`v_1` and :math:`v_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`v_1` and :math:`v_2`. In degrees.

    Returns
    -------
    volume: float
        Volume of corresponding region in space.
    """
    if len(args) == 1:
        v1, v2, v3 = np.array(args[0])
    elif len(args) == 3:
        v1, v2, v3 = np.array(args)
    elif len(args) == 6:
        a, b, c, alpha, beta, gamma = args
        alpha = alpha * _toradians
        beta = beta * _toradians
        gamma = gamma * _toradians
        return (
            a
            * b
            * c
            * sqrt(
                1
                + 2 * cos(alpha) * cos(beta) * cos(gamma)
                - cos(alpha) ** 2
                - cos(beta) ** 2
                - cos(gamma) ** 2
            )
        )
    else:
        raise ValueError(
            "Unable to identify input. "
            + "Supported: one (3,3) array_like, or three (3,) array_like, or 6 floats."
        )

    return v1 @ np.cross(v2, v3)


def winwait():
    r"""
    Add "Press Enter to continue" behaviour to Windows.
    """
    if sys.platform == "win32":
        cprint("Press Enter to continue", "green")
        input()


def print_2D_array(array, fmt="5.2f", posneg=False):
    r"""
    Print 1D and 2D arrays in the terminal.

    .. versionadded:: 0.7

    Parameters
    ----------
    array : (N,) or (N, M) |array_like|_
        Array to be printed. Passed to :numpy:`array`\ () before any action on it.
    fmt : str
        Format string.
    posneg : bool, default False
        Whether to highlight positive and negative values.
        Only works for real-valued arrays.

    Returns
    -------
    string : String with the array, ready to be printed.

    Examples
    --------
    Real-valued array:

    .. doctest::

        >>> import radtools as rad
        >>> array = [[1, 2], [3, 4], [5, 6]]
        >>> rad.print_2D_array(array)
        ┌───────┬───────┐
        │  1.00 │  2.00 │
        ├───────┼───────┤
        │  3.00 │  4.00 │
        ├───────┼───────┤
        │  5.00 │  6.00 │
        └───────┴───────┘

    Custom formatting:

    .. doctest::

        >>> import radtools as rad
        >>> rad.print_2D_array(array, fmt="10.2f")
        ┌────────────┬────────────┐
        │       1.00 │       2.00 │
        ├────────────┼────────────┤
        │       3.00 │       4.00 │
        ├────────────┼────────────┤
        │       5.00 │       6.00 │
        └────────────┴────────────┘

    Scientific notation:

    .. doctest::

        >>> import radtools as rad
        >>> array = [[1, 2], [3, 4], [52414345345, 6]]
        >>> rad.print_2D_array(array, fmt="10.2E")
        ┌────────────┬────────────┐
        │   1.00E+00 │   2.00E+00 │
        ├────────────┼────────────┤
        │   3.00E+00 │   4.00E+00 │
        ├────────────┼────────────┤
        │   5.24E+10 │   6.00E+00 │
        └────────────┴────────────┘

    Complex values:

    .. doctest::

        >>> import radtools as rad
        >>> array = [[1, 2 + 1j], [3, 4], [52, 6]]
        >>> rad.print_2D_array(array)
        ┌────────────────┬────────────────┐
        │  1.00          │  2.00 + i1.00  │
        ├────────────────┼────────────────┤
        │  3.00          │  4.00          │
        ├────────────────┼────────────────┤
        │ 52.00          │  6.00          │
        └────────────────┴────────────────┘
        >>> rad.print_2D_array(array, fmt="4.2E")
        ┌──────────────────────┬──────────────────────┐
        │ 1.00E+00             │ 2.00E+00 + i1.00E+00 │
        ├──────────────────────┼──────────────────────┤
        │ 3.00E+00             │ 4.00E+00             │
        ├──────────────────────┼──────────────────────┤
        │ 5.20E+01             │ 6.00E+00             │
        └──────────────────────┴──────────────────────┘
        >>> array = [[1, 2 - 1j], [3, 4], [52, 6]]
        >>> rad.print_2D_array(array)
        ┌────────────────┬────────────────┐
        │  1.00          │  2.00 - i1.00  │
        ├────────────────┼────────────────┤
        │  3.00          │  4.00          │
        ├────────────────┼────────────────┤
        │ 52.00          │  6.00          │
        └────────────────┴────────────────┘

    Empty arrays:

    .. doctest::

        >>> import radtools as rad
        >>> rad.print_2D_array([])
        None
        >>> rad.print_2D_array([[]])
        None
    """

    array = np.array(array)
    if array.shape[0] != 0 and array.shape[1] != 0:
        if len(array.shape) == 1:
            array = np.array([array])
        N = len(array)
        M = len(array[0])
        n = max(
            len(f"{np.amax(array.real):{fmt}}"), len(f"{np.amin(array.real):{fmt}}")
        )
        n = max(
            n, len(f"{np.amax(array.imag):{fmt}}"), len(f"{np.amin(array.imag):{fmt}}")
        )
        if "E" in fmt or "e" in fmt:
            n = int(fmt.split(".")[0])
        else:
            n = max(n, int(fmt.split(".")[0]))
        fmt = f"{n}.{fmt.split('.')[1]}"
        complex_values = not np.isreal(array).all()
        if complex_values:
            nn = 2 * n + 6
            if "E" in fmt or "e" in fmt:
                nn += 8
        else:
            nn = n + 2
        print("┌" + (M - 1) * f"{nn*'─'}┬" + f"{nn*'─'}┐")
        for i in range(0, N):
            print("│", end="")
            for j in range(0, M):
                if complex_values:
                    if np.iscomplex(array[i][j]):
                        if array.imag[i][j] >= 0:
                            sign = "+"
                        else:
                            sign = "-"
                        print(
                            f" {array.real[i][j]:{fmt}} {sign} i{abs(array.imag[i][j]):<{fmt}} │",
                            end="",
                        )
                    elif "E" in fmt or "e" in fmt:
                        print(
                            f" {array.real[i][j]:{fmt}}{(n + 8)*' '} │",
                            end="",
                        )
                    else:
                        print(
                            f" {array.real[i][j]:{fmt}}{(n + 4)*' '} │",
                            end="",
                        )
                else:
                    if posneg:
                        if array.real[i][j] > 0:
                            cprint(f" {array.real[i][j]:{fmt}}", "red", end="")
                            print(" │", end="")
                        elif array.real[i][j] < 0:
                            cprint(f" {array.real[i][j]:{fmt}}", "blue", end="")
                            print(" │", end="")
                        else:
                            cprint(f" {array.real[i][j]:{fmt}}", "green", end="")
                            print(" │", end="")
                    else:
                        print(f" {array.real[i][j]:{fmt}} │", end="")
            print()
            if i != N - 1:
                print("├" + (M - 1) * f"{nn*'─'}┼" + f"{nn*'─'}┤")

        print("└" + (M - 1) * f"{nn*'─'}┴" + f"{nn*'─'}┘")
    else:
        print("None")


def angle(v1, v2, radians=False):
    r"""
    Angle between two vectors.

    .. versionadded:: 0.7

    Parameters
    ----------
    v1 : (3,) |array_like|_
        First vector.
    v2 : (3,) |array_like|_
        Second vector.
    radians : bool, default False
        Whether to return value in radians. Return value in degrees by default.

    Returns
    -------
    angle: float
        Angle in degrees or radians (see ``radians``).
    """

    v1 = np.array(v1) / np.linalg.norm(v1)
    v2 = np.array(v2) / np.linalg.norm(v2)
    alpha = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    if radians:
        return alpha
    return alpha * _todegrees


def reciprocal_cell(cell):
    r"""
    Computes reciprocal cell.

    .. versionadded:: 0.7

    Parameters
    ----------
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.

    Returns
    -------
    reciprocal_cell : (3,3) :numpy:`ndarray`
        Reciprocal cell matrix, rows are interpreted as vectors.
        :math:`cell = (\vec{v}_1, \vec{v}_2, \vec{v}_3)`, where

        .. math::

            \begin{matrix}
                \vec{b}_1 = \frac{2\pi}{V}\vec{a}_2\times\vec{a}_3 \\
                \vec{b}_2 = \frac{2\pi}{V}\vec{a}_3\times\vec{a}_1 \\
                \vec{b}_3 = \frac{2\pi}{V}\vec{a}_1\times\vec{a}_2 \\
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


def cell_from_param(a=1, b=1, c=1, alpha=90, beta=90, gamma=90):
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
    cell : (3,3) :numpy:`ndarray`
        Cell matrix.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    """

    alpha = alpha * _toradians
    beta = beta * _toradians
    gamma = gamma * _toradians
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
        ]
    )


def param_from_cell(cell):
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


def get_permutation(n, k):
    r"""
    Return array of index permutations

    .. versionadded:: 0.7

    Parameters
    ----------
    n : int
        Amount of index to be used:

        .. code-block::

            range(0, n)
    k : int
        Length of the permuted arrays.

    Returns
    -------
    permutations : list
        List of permutations. If N permutations found,
        the it is a list of N lists of length k.

    """
    if k == n:
        return [[i for i in range(n)]]
    elif n < k:
        raise ValueError("Permutations: n < k")
    else:
        result = [[i] for i in range(n)]
        true_k = k
        k = 0
        while k < true_k - 1:
            new_result = []
            for i in range(len(result)):
                for j in range(result[i][-1] + 1, n):
                    new_result.append(result[i] + [j])
            result = new_result
            k += 1
        return result
