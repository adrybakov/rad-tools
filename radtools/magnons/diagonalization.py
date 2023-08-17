import numpy as np
from numpy.linalg import LinAlgError

from radtools.exceptions import ColpaFailed

__all__ = ["solve_via_colpa"]


def solve_via_colpa(D):
    r"""
    Diagonalize grand-dynamical matrix following the method of Colpa [1]_.

    Algorithm itself is described in section 3, Remark 1 of [1]_.

    Solves the Bogoliubov Hamiltonian of the form:

    .. math::

        \hat{H} = \sum_{r^{\prime}, r = 1}^m 
        \hat{\alpha}_{r^{\prime}}^{\dagger}\boldsymbol{\Delta}_1^{r^{\prime}r}\hat{\alpha}_r +
        \hat{\alpha}_{r^{\prime}}^{\dagger}\boldsymbol{\Delta}_2^{r^{\prime}r}\hat{\alpha}_{m+r}^{\dagger} + 
        \hat{\alpha}_{m+r^{\prime}}^{\dagger}\boldsymbol{\Delta}_3^{r^{\prime}r}\hat{\alpha}_r + 
        \hat{\alpha}_{m+r^{\prime}}^{\dagger}\boldsymbol{\Delta}_4^{r^{\prime}r}\hat{\alpha}_{m+r}^{\dagger}

    In a matrix form the Hamiltonian is:

    .. math::

        \hat{H} = \boldsymbol{\hat{a}}^{\dagger} \boldsymbol{D} \boldsymbol{\hat{a}}
    
    where

    .. math:: 

        \boldsymbol{\hat{a}} =
        \begin{pmatrix}
            \hat{\alpha}_1 \\
            \cdots \\
            \hat{\alpha}_m \\
            \hat{\alpha}_{m+1} \\
            \cdots \\
            \hat{\alpha}_{2m} \\
        \end{pmatrix}

    After diagonalization the Hamiltonian is:

    .. math::

        \hat{H} = \boldsymbol{\hat{c}}^{\dagger} \boldsymbol{E} \boldsymbol{\hat{c}}

    Parameters
    ----------
    D : (2N, 2N) |array_like|_
        Grand dynamical matrix. Must be Hermitian and positive-defined.

        .. math::

            \boldsymbol{D} = \begin{pmatrix}
                \boldsymbol{\Delta_1} & \boldsymbol{\Delta_2} \\
                \boldsymbol{\Delta_3} & \boldsymbol{\Delta_4}
            \end{pmatrix}

    Returns
    -------
    E : (2N,) :numpy:`ndarray`
        The eigenvalues, each repeated according to its multiplicity.
        First N eigenvalues are sorted in descending order.
        Last N eigenvalues are sorted in ascending order.
        In the case of diagonalization of the magnon Hamiltonian
        first N eigenvalues are the same as last N eigenvalues, but
        in reversed order. It is an array of the diagonal elements of the 
        diagonal matrix :math:`\boldsymbol{E}` from the diagonalized Hamiltonian.

    G : (2N, 2N) :numpy:`ndarray`
        Transformation matrix, which change the basis from the original set of bosonic 
        operators :math:`\boldsymbol{\hat{a}}` to the set of 
        new bosonic operators :math:`\boldsymbol{\hat{c}}` which diagonalize
        the Hamiltonian:
        
        .. math::

            \boldsymbol{\hat{c}} = \boldsymbol{G} \boldsymbol{\hat{a}}

    Notes
    -----

    Let :math:`\boldsymbol{E}` be the diagonal matrix of eigenvalues ``E``, then:

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
