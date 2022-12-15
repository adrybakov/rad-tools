r"""
Collection of small routines and constants, 
which may be used across the whole package.
"""

import numpy as np
from math import asin, sqrt, pi

# Teminal colors
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
        raise ValueError(f'Integer n have to be in range 0 <= n <= 255. '
                         f'You provided n = {n}, type<n> = {type(n)}')
    return f'\033[38:5:{n}m'


def exchange_to_matrix(iso=None, aniso=None, dmi=None):
    r"""
    Combine isotropic, anisotropic and dmi exchange into exchange matrix.

    Parameters
    ----------
    iso : float
        Value of isotropic exchange parameter in meV.
    aniso : 3 x 3 array_like
        Matrix of symmetric anisotropic exchange in meV.
    dmi : 1 x 3 array_like
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV.

    Returns
    -------
    matrix : 3 x 3 np.ndarray of floats
        Exchange matrix in meV.
    """
    matrix = np.zeros((3, 3), dtype=float)
    if aniso is not None:
        matrix += np.array(aniso, dtype=float)
    if iso is not None:
        matrix += iso * np.identity(3, dtype=float)
    if dmi is not None:
        matrix += np.array([[0, dmi[2], -dmi[1]],
                            [-dmi[2], 0, dmi[0]],
                            [dmi[1], -dmi[0], 0]],
                           dtype=float)
    return matrix


def exchange_from_matrix(matrix):
    r"""
    Decompose matrix into isotropic, anisotropic and dmi exchange.

    Parameters
    ----------
    matrix : 3 x 3 array_like
        Exchange matrix in meV.

    Returns
    -------
    iso : float
        Value of isotropic exchange parameter in meV.
    aniso : 3 x 3 np.ndarray of floats
        Matrix of symmetric anisotropic exchange in meV.
    dmi : 1 x 3 np.ndarray of floats
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV.
    """
    matrix = np.array(matrix, dtype=float)
    symm = (matrix + matrix.T) / 2
    assym = (matrix - matrix.T) / 2
    dmi = np.array([assym[1][2], assym[2][0], assym[0][1]],
                   dtype=float)
    iso = np.trace(symm) / 3
    aniso = symm - iso * np.identity(3, dtype=float)
    return iso, aniso, dmi


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
    numbers = '0123456789'
    new_mark = '$'
    insert_underline = False
    for symbol in mark:
        if symbol in numbers and not insert_underline:
            insert_underline = True
            new_mark += '_{'
        new_mark += symbol
    new_mark += '}$'
    return new_mark


def rot_angle(x, y, dummy=False):
    r"""
    Rotational ange from 2D vector.

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
        raise ValueError('Angle is ill defined (x = y = 0).')
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
            raise ValueError('Angle is ill defined (x = y = 0).')
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


def strip_digits(line: str):
    r"""
    Remove all digits from the string

    Parameters
    ----------
    line : str
        Input string.

    Returns
    -------
    new_line : str
        ``line`` without digits
    """

    new_line = ""
    numbers = "1234567890"
    for char in line:
        if char not in numbers:
            new_line += char
    return new_line


def two_points_distance(point1, point2):
    r"""
    Compute distance between two points.

    Parameters
    ----------
    point1 : array
        Coordinates of the first point. 

        .. code-block:: python

            [x1, y1, z1]

    point2 : array
        Coordinates of the second point. 

        .. code-block:: python

            [x2, y2, z2]

    Returns
    -------
    distance : float
        Distance between two points.
    """

    point1 = np.array(point1)
    point2 = np.array(point2)
    return sqrt(np.sum((point1 - point2)**2))


def search_on_atoms(centre, atoms):
    r"""
    Search the closest atom to the given centre position.

    Parameters
    ----------
    centre : array
        xyz coordinates of the centre.
    atoms : list
        List of atom names and coordinates of the form: ::

            [(name, [x, y, z]), ...]

    Returns
    -------
    min_span : float
        Distance from the centre to the closest atom.
    name : str
        Name of the closest atom.
    """

    min_span = 10000
    name = "None"
    for atom, a_coord in atoms:
        if sqrt(np.sum((centre[1] - a_coord)**2)) < min_span:
            min_span = sqrt(np.sum((centre[1] - a_coord)**2))
            name = atom
    return min_span, name


def search_between_atoms(centre, atoms):
    r"""
    Search the closest bond centre to the given centre position.

    Parameters
    ----------
    centre : array
        xyz coordinates of the centre.
    atoms : list
        List of atom names and coordinates of the form: ::

            [(name, [x, y, z]), ...]

    Returns
    -------
    min_span : float
        Distance from the centre to the bond`s centre.
    name : str
        Name of the closest bond`s centre (atom1-atom2).
    """

    pairs = []
    for i, atom in enumerate(atoms):
        for j in range(i+1, len(atoms)):
            pair = f"{atom[0]}-{atoms[j][0]}"
            p_coord = (atom[1]+atoms[j][1])/2
            pairs.append([pair, p_coord])
    return search_on_atoms(centre, pairs)
