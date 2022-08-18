from os import mkdir
from typing import Union


class TerminalCoulours:
    """
    Class with constants for ANSI escape codes
    """

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


def check_make(path: str):
    """
    Check if directory exist, create one if not.

    Parameters
    ----------
    path : str
        path to the desired directory
    """

    try:
        mkdir(path)
    except FileExistsError:
        pass


def matrix_exchange(J_iso: Union[float, int] = None,
                    J_aniso: list = None,
                    dmi: Union[tuple, list] = None):
    """
    Recompute different types of exchange into a common matrix form.

    Matrix from of isotropic exchange:
    |J_iso   0     0  |
    |  0   J_iso   0  |
    |  0     0   J_iso|

    Symmetric anisotropic exchange already in the matrix form:
    |J_xx J_xy J_xz|
    |J_yx J_yy J_yz|
    |J_zx J_zy J_zz|
    Note: J_ij = J_ji

    Matrix from of DMI:
    |  0   D_z -D_y|
    |-D_z   0   D_x|
    | D_y -D_x   0 |

    Parameters
    ----------
    J_iso : float | int
        Isotropic exchange parameter.
    J_aniso : list of list of float
        Symmetric anisotropic exchange 3x3 matrix.
    dmi : tuple | list of float
        Dzyaloshinskii-Moriya interaction vector.

    Returns
    -------
    J_matrix : list of list of float
        Matrix form of exchange parameter.
    """

    J_matrix = [[0., 0., 0.],
                [0., 0., 0.],
                [0., 0., 0.]]
    if J_iso is not None:
        for i in range(0, 3):
            J_matrix[i][i] += J_iso

    if J_aniso is not None:
        for i in range(0, 3):
            for j in range(0, 3):
                J_matrix[i][j] += J_aniso[i][j]
    if dmi is not None:
        J_matrix[0][1] += dmi[2]  # D_z
        J_matrix[0][2] -= dmi[1]  # -D_y
        J_matrix[1][0] -= dmi[2]  # -D_z
        J_matrix[1][2] += dmi[0]  # D_x
        J_matrix[2][0] += dmi[1]  # D_y
        J_matrix[2][1] -= dmi[0]  # -D_x
    return J_matrix
