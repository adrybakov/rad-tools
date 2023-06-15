import numpy as np
from numpy.linalg import LinAlgError


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
    omega : (2N,) :numpy:`ndarray`
        The eigenvalues, each repeated according to its multiplicity.
        The eigenvalues are sorted in descending order.
        See: :numpy:`linalg.eig` for more details.

    U : (2N, 2N) :numpy:`ndarray`
        The normalized (unit “length”) eigenvectors, such that the column
        U[:,i] is the eigenvector corresponding to the eigenvalue omega[i].
        See: :numpy:`linalg.eig` for more details.

    References
    ----------
    .. [1] Colpa, J.H.P., 1978.
        Diagonalization of the quadratic boson hamiltonian.
        Physica A: Statistical Mechanics and its Applications,
        93(3-4), pp.327-353.
    """

    try:
        K = np.linalg.cholesky(h)
    except LinAlgError:
        raise ColpaFailed

    N = len(h) // 2
    g = np.diag(np.concatenate((np.ones(N), -np.ones(N))))

    omegas, U = np.linalg.eig(K @ g @ np.conjugate(K).T)

    # sort with respect to omegas, in descending order
    U = np.concatenate((omegas[:, None], U.T), axis=1).T
    U = U[:, np.argsort(U[0])]
    omegas = U[0, ::-1]
    U = U[1:, ::-1]

    return omegas, U
