import numpy as np


class ColpaFailed(Exception):
    r"""
    Raise when Diagonalization via Colpa fails.
    """

    def __init__(self):
        self.message = "Diagonalization via Colpa failed."

    def __str__(self):
        return self.message


def solve_via_colpa(h):
    r"""
    Diagonalize grand-dynamical matrix following the method of Colpa [1]_.

    Parameters
    ----------
    h : (2N, 2N) |array_like|_

    Returns
    -------
    omegas : (2N,) :numpy:`ndarray`
        Eigenvalues first N eigenvalues are the same as second N.
    U : (2N, 2N) : :numpy:`ndarray`
        Eigenvectors.

    References
    ----------
    .. [1] Colpa, J.H.P., 1978.
        Diagonalization of the quadratic boson hamiltonian.
        Physica A: Statistical Mechanics and its Applications,
        93(3-4), pp.327-353.
    """

    try:
        K = np.linalg.cholesky(h)
    except:
        raise ColpaFailed

    N = len(h) // 2
    g = np.eye(2 * N)
    for i in range(N, 2 * N):
        g[i][i] = -1

    omegas, U = np.linalg.eig(K @ g @ np.conjugate(K).T)
    return g @ omegas, U
