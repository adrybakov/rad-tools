r"""
Collection of small routines and constants, 
which are used across the whole package.

It's purpose is to serve as an "other" folder.
"""

import sys
from math import asin, cos, pi, sin, sqrt

import numpy as np
from termcolor import cprint, colored

__all__ = [
    "print_2d_array",
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

# Constants are usually defined in uppercase
# but these two are intentionally defined in lowercase
todegrees = 180 / pi
toradians = pi / 180


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
    cell : (3, 3) |array_like|_
        Lattice vectors.
    absolute : (3,) |array_like|_
        Absolute coordinates.

    Returns
    -------
    relative : (3,) :numpy:`ndarray`
        Relative coordinates.
    """

    # Three vectors of the cell
    v1 = np.array(cell[0], dtype=float)
    v2 = np.array(cell[1], dtype=float)
    v3 = np.array(cell[2], dtype=float)

    v = np.array(absolute, dtype=float)
    if (v == np.zeros(3)).all():
        return np.zeros(3)

    # Compose system of linear equations
    B = np.array([np.dot(v1, v), np.dot(v2, v), np.dot(v3, v)])
    A = np.array(
        [
            [np.dot(v1, v1), np.dot(v1, v2), np.dot(v1, v3)],
            [np.dot(v2, v1), np.dot(v2, v2), np.dot(v2, v3)],
            [np.dot(v3, v1), np.dot(v3, v2), np.dot(v3, v3)],
        ]
    )

    # Solve and return
    return np.linalg.solve(A, B)


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
            V = abc\sqrt{1+2\cos(\alpha)\cos(\beta)\cos(\gamma)-\cos^2(\alpha)-\cos^2(\beta)-\cos^2(\gamma)}

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
    volume : float
        Volume of corresponding region in space.
    """

    if len(args) == 1:
        v1, v2, v3 = np.array(args[0])
    elif len(args) == 3:
        v1, v2, v3 = np.array(args)
    elif len(args) == 6:
        a, b, c, alpha, beta, gamma = args
        alpha = alpha * toradians
        beta = beta * toradians
        gamma = gamma * toradians
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

    Its a hotfix for Window`s terminal behaviour.
    """
    if sys.platform == "win32":
        cprint("Press Enter to continue", "green")
        input()


def print_2d_array(
    array, fmt=".2f", highlight=False, print_result=True, borders=True, shift=0
):
    r"""
    Print 1D and 2D arrays in the terminal.

    .. versionadded:: 0.7

    .. versionchanged:: 0.7.11 Renamed from ``print_2D_array``

    Parameters
    ----------
    array : (N,) or (N, M) |array_like|_
        Array to be printed.
    fmt : str
        Format string.
    highlight : bool, default False
        Whether to highlight positive and negative values.
        Only works for real-valued arrays.

        .. versionchanged:: 0.7.11 Renamed from ``posneg``

    print_result : bool, default True
        Whether to print the result or return it as a string.
    borders : bool, default True
        Whether to print borders around the array.

        .. versionadded:: 0.7.11

    shift : int, default 0
        Shifts the array to the right by ``shift`` columns.

        .. versionadded:: 0.7.11

    Returns
    -------
    string : str
        String representation of the array.
        Returned only if ``print_result`` is False.
    """

    array = np.array(array)
    if (len(array.shape) == 1 and array.shape[0] != 0) or (
        len(array.shape) == 2 and array.shape[1] != 0
    ):
        # Convert 1D array to 2D array
        if len(array.shape) == 1:
            array = np.array([array])

        # Array dimensions
        N = len(array)
        M = len(array[0])

        # Find the longest number, used for string formatting
        n = max(
            len(f"{np.amax(array.real):{fmt}}"), len(f"{np.amin(array.real):{fmt}}")
        )
        n = max(
            n, len(f"{np.amax(array.imag):{fmt}}"), len(f"{np.amin(array.imag):{fmt}}")
        )

        # Check if input fmt exceeds the longest number
        try:
            n = max(n, int(fmt.split(".")[0]))
        except ValueError:
            pass

        fmt = f"{n}.{fmt.split('.')[1]}"

        def print_number(number, fmt, highlight=False, condition=None):
            if condition is None:
                condition = number
            string = ""
            # Highlight positive and negative values with colours
            if highlight:
                if condition > 0:
                    string += colored(f"{number:{fmt}}", "red", attrs=["bold"])
                elif condition < 0:
                    string += colored(f"{number:{fmt}}", "blue", attrs=["bold"])
                else:
                    string += colored(f"{number:{fmt}}", "green", attrs=["bold"])
            # Print without colours
            else:
                string += f"{number:{fmt}}"
            return string

        def print_complex(number, fmt, highlight):
            string = ""
            if number.real != 0:
                string += print_number(number.real, fmt, highlight)
            else:
                string += " " * len(print_number(number.real, fmt))

            if number.imag > 0:
                sign = "+"
            else:
                sign = "-"

            if number.imag != 0:
                string += f" {sign} i"
                string += print_number(
                    abs(number.imag), f"<{fmt}", highlight, condition=number.imag
                )
            else:
                string += " " * (
                    len(print_number(abs(number.imag), fmt, condition=number.imag)) + 4
                )
            return string

        def print_border(symbol_start, symbol_middle, symbol_end, n):
            string = symbol_start
            for j in range(0, M):
                # If at least one complex value is present in the column
                if np.iscomplex(array[:, j]).any():
                    string += f"{(2*n + 6)*'─'}"
                else:
                    string += f"{(n + 2)*'─'}"
                if j == M - 1:
                    string += symbol_end + "\n"
                else:
                    string += symbol_middle
            return string

        string = ""
        for i in range(0, N):
            substring = " " * shift
            if borders:
                substring += "│"

            for j in range(0, M):
                # Print complex values
                if np.iscomplex(array[:, j]).any():
                    # Print complex part if it is non-zero
                    substring += " " + print_complex(array[i][j], fmt, highlight)
                # Print real values
                else:
                    # Highlight positive and negative values with colours
                    substring += " " + print_number(array[i][j].real, fmt, highlight)
                if borders:
                    substring += " │"
            if i != N - 1 or borders:
                substring += "\n"

            if borders:
                # Header of the table
                if i == 0:
                    symbol_start = "┌"
                    symbol_middle = "┬"
                    symbol_end = "┐"
                    substring = (
                        " " * shift
                        + print_border(symbol_start, symbol_middle, symbol_end, n)
                        + substring
                    )

                # Footer of the table
                if i == N - 1:
                    symbol_start = "└"
                    symbol_middle = "┴"
                    symbol_end = "┘"
                    substring += (
                        " " * shift
                        + print_border(symbol_start, symbol_middle, symbol_end, n)[:-1]
                    )
                # Middle of the table
                else:
                    symbol_start = "├"
                    symbol_middle = "┼"
                    symbol_end = "┤"
                    substring += " " * shift + print_border(
                        symbol_start, symbol_middle, symbol_end, n
                    )

            string += substring
        if print_result:
            print(string)
        else:
            return string
    else:
        if print_result:
            print(None)
        else:
            return None


def angle(v1, v2, radians=False):
    r"""
    Angle between two vectors.

    .. versionadded:: 0.7

    .. math::

        \alpha = \dfrac{(\vec{v_1} \cdot \vec{v_2})}{\vert\vec{v_1}\vert\cdot\vert\vec{v_2}\vert}

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

    # Normalize vectors
    v1 = np.array(v1) / np.linalg.norm(v1)
    v2 = np.array(v2) / np.linalg.norm(v2)

    alpha = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    if radians:
        return alpha
    return alpha * todegrees


def reciprocal_cell(cell):
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
    cell : (3, 3) :numpy:`ndarray`
        Cell matrix.

        .. code-block:: python

            cell = [[a1_x, a1_y, a1_z],
                    [a2_x, a2_y, a2_z],
                    [a3_x, a3_y, a3_z]]
    """
    alpha = alpha * toradians
    beta = beta * toradians
    gamma = gamma * toradians
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
        Length of the array to be permuted.

        .. code-block::

            range(0, n)
    k : int
        Length of the permutation.

    Returns
    -------
    permutations : list
        List of permutations. If N permutations are found,
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


if __name__ == "__main__":
    print_2d_array([[1, 2], [3, 4], [5, 6]])
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt="10.2f")
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt=".2f")
    print_2d_array([[1, 2], [3, 4], [52414345345, 6]], fmt="10.2E")
    print_2d_array([[1, 2 + 1j], [3, 4], [52, 6]])
    print_2d_array([[1, 2 - 1j], [3, 4], [52, 6]])
    print_2d_array([[1, -1j], [3, 4], [52, 6]])
    print_2d_array([])

    print_2d_array([[1, 2], [3, 4], [5, 6]], highlight=True)
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt="10.2f", highlight=True)
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt=".2f", highlight=True)
    print_2d_array([[1, 2], [3, 4], [52414345345, 6]], fmt="10.2E", highlight=True)
    print_2d_array([[1, 2 + 1j], [3, 4], [52, 6]], highlight=True)
    print_2d_array([[1, 2 - 1j], [3, 4], [52, 6]], highlight=True)
    print_2d_array([[1, -1j], [3, 4], [52, 6]], highlight=True)

    print_2d_array([[1, 2], [3, 4], [5, 6]], highlight=True, borders=False)
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt="10.2f", highlight=True, borders=False)
    print_2d_array([[1, 2], [3, 4], [5, 6]], fmt=".2f", highlight=True, borders=False)
    print_2d_array(
        [[1, 2], [3, 4], [52414345345, 6]], fmt="10.2E", highlight=True, borders=False
    )
    print_2d_array([[1, 2 + 1j], [3, 4], [52, 6]], highlight=True, borders=False)
    print_2d_array([[1, 2 - 1j], [3, 4], [52, 6]], highlight=True, borders=False)
    print_2d_array([[1, -1j], [3, 4], [52, 6]], highlight=True, borders=False)

    print_2d_array([[1, 2], [3, 4], [5, 6]], highlight=True, borders=False, shift=8)
    print_2d_array([[1, 2], [3, 4], [5, 6]], highlight=True, shift=8)
