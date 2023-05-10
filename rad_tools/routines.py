r"""
Collection of small routines and constants, 
which may be used across the whole package.
"""

import sys

from math import asin, pi, sqrt

import numpy as np

# Terminal colours
BLACK = "\u001b[30m"
"""
ANSI escape code for black color of the text.
"""
RED = "\u001b[31m"
"""
ANSI escape code for red color of the text.
"""
GREEN = "\u001b[32m"
"""
ANSI escape code for green color of the text.
"""
YELLOW = "\u001b[33m"
"""
ANSI escape code for yellow color of the text.
"""
BLUE = "\u001b[34m"
"""
ANSI escape code for blue color of the text.
"""
MAGENTA = "\u001b[35m"
"""
ANSI escape code for magenta color of the text.
"""
CYAN = "\u001b[36m"
"""
ANSI escape code for cyan color of the text.
"""
WHITE = "\u001b[37m"
"""
ANSI escape code for white color of the text.
"""
RESET = "\u001b[0m"
"""
ANSI escape code for resetting color to default.
"""
WARNING = YELLOW
"""
ANSI escape code for warnings.
"""
OK = GREEN
"""
ANSI escape code for ok messages.
"""
ERROR = RED
"""
ANSI escape code for errors.
"""

__all__ = [
    "get_named_colours",
    "get_256_colours",
    "cprint",
    "atom_mark_to_latex",
    "rot_angle",
    "absolute_to_relative",
    "winwait",
    "print_2D_array",
]


def get_named_colours(colour: str):
    r"""
    Get the colours for the terminal based on name.

    Parameters
    ----------
    colour : str
        Name of the colour.
        Choices: black, red, green, yellow, blue, magenta, cyan, white.
    """
    colour_dict = {
        "black": BLACK,
        "red": RED,
        "green": GREEN,
        "yellow": YELLOW,
        "blue": BLUE,
        "magenta": MAGENTA,
        "cyan": CYAN,
        "white": WHITE,
    }
    colour = colour.lower()
    if colour not in colour_dict:
        return ""
    return colour_dict[colour]


def get_256_colours(n):
    r"""
    ANSI escape codes for terminal color with 256-colours support
    (see: |ANSI|_).

    Parameters
    ----------
    n : int
        Integer from 0 to 255 to be mapped to colours.

    Returns
    -------
    str
        string with the ANSI escape code.
    """

    if type(n) != int or not 0 <= n <= 255:
        raise ValueError(
            f"Integer n have to be in range 0 <= n <= 255. "
            f"You provided n = {n}, type<n> = {type(n)}"
        )
    return f"\033[38:5:{n}m"


def cprint(*args, colour=None, **kwargs):
    r"""
    Add colour argument to the standard ``print()``.

    Parameters
    ----------
    colour : str or int
        Name or number for a colour. Number is used with the base of 256.
        Name should comply with :py:func:`.get_named_colours`.

    Example
    -------
    >>> import rad_tools as rad
    >>> rad.cprint("Hellow world!", colour="green")
    Hellow world!
    >>> rad.cprint("Hellow world!", colour = 30)
    Hellow world!
    >>> rad.cprint("Hellow world!", colour=1342)
    Hellow world!
    """

    if isinstance(colour, int):
        colour = get_256_colours(colour % 256)
    elif isinstance(colour, str):
        colour = get_named_colours(colour)
    else:
        raise ValueError("Colour must be a string or integer")
    if colour is not None:
        print(colour, end="")
    print(*args, **kwargs)
    if colour is not None:
        print(RESET, end="")


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


def absolute_to_relative(cell, x, y, z):
    r"""
    Compute relative coordinates with respect to the unit cell.

    Parameters
    ----------
    cell : 3 x 3 array
        Lattice vectors.
    x : float
        x coordinate.
    y : float
        y coordinate.
    z : float
        z coordinate.

    Returns
    -------
    relative : 1 x 3 array
        Relative coordinate.
    """

    a = np.array(cell[0], dtype=float)
    b = np.array(cell[1], dtype=float)
    c = np.array(cell[2], dtype=float)
    v = np.array([x, y, z], dtype=float)
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


def winwait():
    r"""
    Add "Press Enter to continue" behaviour to Windows.
    """
    if sys.platform == "win32":
        cprint("Press Enter to continue", colour="green")
        input()


def print_2D_array(array, fmt="5.2f"):
    r"""
    Print 2D array in the terminal.

    Parameters
    ----------
    array : array
        Array to be printed. Passed to ``np.array()`` before any action on it.
    fmt : str
        Format string.

    Returns
    -------
    string : String with the array, ready to be printed.

    Examples
    --------
    Real-valued array:

    .. doctest::

        >>> import rad_tools as rad
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

        >>> import rad_tools as rad
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

        >>> import rad_tools as rad
        >>> array = [[1, 2], [3, 4], [52435345345, 6]]
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

        >>> import rad_tools as rad
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
    """

    array = np.array(array)
    N = len(array)
    M = len(array[0])
    n = max(len(f"{np.amax(array.real):{fmt}}"), len(f"{np.amin(array.real):{fmt}}"))
    n = max(n, len(f"{np.amax(array.imag):{fmt}}"), len(f"{np.amin(array.imag):{fmt}}"))
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
                print(f" {array.real[i][j]:{fmt}} │", end="")
        print()
        if i != N - 1:
            print("├" + (M - 1) * f"{nn*'─'}┼" + f"{nn*'─'}┤")

    print("└" + (M - 1) * f"{nn*'─'}┴" + f"{nn*'─'}┘")


if __name__ == "__main__":
    a = [[1, 2], [3, 4], [5, 6]]
    print_2D_array(a)
    print_2D_array(a, fmt="10.2f")
    a = [[1, 2], [3, 4], [52435345345, 6]]
    print_2D_array(a, fmt="10.2E")
    a = [[1, 2 + 1j], [3, 4], [52, 6]]
    print_2D_array(a)
    print_2D_array(a, fmt="4.2E")
