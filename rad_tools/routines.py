from os import mkdir
from os.path import split, isdir, join
from typing import Union
import numpy as np

# Teminal colors
BLACK = '\u001b[30m'
RED = '\u001b[31m'
GREEN = '\u001b[32m'
YELLOW = '\u001b[33m'
BLUE = '\u001b[34m'
MAGENTA = '\u001b[35m'
CYAN = '\u001b[36m'
WHITE = '\u001b[37m'
RESET = ' \u001b[0m'

WARNING = YELLOW
OK = GREEN
ERROR = RED


def get_256_colours(n: int):
    """
    ANSI escape codes for terminal color with 256-colours support
    (see: https://en.wikipedia.org/wiki/ANSI_escape_code).

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


def check_make_dir(path: str):
    """
    Check if directory exist, create one if not.

    If parent dir does not exist this script will create one, recursively. 

    Parameters
    ----------
    path : str
        path to the desired directory
    """
    head = path
    tail = ''
    dirs_to_make = []
    while not isdir(head) and head:
        try:
            mkdir(head)
        except FileNotFoundError:
            head, tail = split(head)
            dirs_to_make.append(tail)

    dirs_to_make.reverse()
    for dir in dirs_to_make:
        head = join(head, dir)
        mkdir(head)


def exchange_to_matrix(iso=None, aniso=None, dmi=None):
    """
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
    """
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
