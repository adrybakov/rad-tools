r"""
Collection of small routines and constants, 
which may be used across the whole package.
"""

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


def get_named_colours(colour: str):
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
    (see: :ANSI:`wiki <>`).

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
