import numpy as np
from numpy.linalg import LinAlgError

__all__ = ["solve_via_colpa"]


class ColpaFailed(Exception):
    r"""
    Raised when Diagonalization via Colpa fails.
    """

    def __init__(self):
        self.message = "Diagonalization via Colpa failed."

    def __str__(self):
        return self.message


def solve_via_colpa(D):
    r"""
    Diagonalize grand-dynamical matrix following the method of Colpa [1]_.

    Algorithm itself is described in section 3, Remark 1 of [1]_.

    Parameters
    ----------
    D : (2N, 2N) |array_like|_
        Grand dynamical matrix. Must be Hermitian and positive-defined.

    Returns
    -------
    E : (2N,) :numpy:`ndarray`
        The eigenvalues, each repeated according to its multiplicity.
        First N eigenvalues are sorted in descending order.
        Last N eigenvalues are sorted in ascending order.
        In the case of diagonalization of the magnon Hamiltonian
        first N eigenvalues are the same as last N eigenvalues, but
        in reversed order.

    G : (2N, 2N) :numpy:`ndarray`
        Transformation matrix, which change the basis from the original set of bosonic 
        operators :math:`\boldsymbol{a}_{\boldsymbol{k}}`to the set of 
        new bosonic operators :math:`\boldsymbol{c}_{\boldsymbol{k}}` which diagonalize
        the Hamiltonian which corresponds to the grand-dynamical matrix ``D``:
        
        .. math::

            \boldsymbol{c}_{\boldsymbol{k}} = \boldsymbol{G} \boldsymbol{a}_{\boldsymbol{k}}

        where (same for :math:`\boldsymbol{c}_{\boldsymbol{k}}`)

        .. math:: 

            \boldsymbol{a}_{\boldsymbol{k}} =
            \begin{pmatrix}
                \boldsymbol{\alpha}_{\boldsymbol{k}} \\
                \boldsymbol{\alpha}_{-\boldsymbol{k}}^{\dagger}
            \end{pmatrix}

    Notes
    -----

    Let :math:`\boldsymbol{E}` be the diagonal matrix of eigenvalues ``E``.

    .. math::

        \boldsymbol{E} = (\boldsymbol{G}^{\dagger})^{-1} \boldsymbol{D} \boldsymbol{G}^{-1}



    References
    ----------
    .. [1] Colpa, J.H.P., 1978.
        Diagonalization of the quadratic boson hamiltonian.
        Physica A: Statistical Mechanics and its Applications,
        93(3-4), pp.327-353.
    """

    D = np.array(D)

    N = len(D) // 2
    g = np.diag(np.concatenate((np.ones(N), -np.ones(N))))

    try:
        # In Colpa article decomposition is K^{\dag}K, while numpy gives KK^{\dag}
        K = np.conjugate(np.linalg.cholesky(D)).T
    except LinAlgError:
        raise ColpaFailed

    L, U = np.linalg.eig(K @ g @ np.conjugate(K).T)

    # Sort with respect to L, in descending order
    U = np.concatenate((L[:, None], U.T), axis=1).T
    U = U[:, np.argsort(U[0])]
    L = U[0, ::-1]
    U = U[1:, ::-1]

    E = g @ L

    G_minus_one = np.linalg.inv(K) @ U @ np.sqrt(np.diag(E))

    # Compute G from G^-1 following Colpa, see equation (3.7) for details
    G = np.conjugate(G_minus_one).T
    G[:N, N:] *= -1
    G[N:, :N] *= -1

    return E, G
