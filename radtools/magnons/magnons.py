from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.routines import span_orthonormal_set
from radtools.crystal.atom import Atom
from radtools.magnons.diagonalization import solve_via_colpa, ColpaFailed
from radtools.routines import print_2d_array

from copy import deepcopy


from scipy.spatial.transform import Rotation
import numpy as np
from math import sqrt, pi

import matplotlib.pyplot as plt

__all__ = ["MagnonDispersion"]


class MagnonDispersion:
    r"""
    Magnon dispersion wrapper.

    Parameters
    ----------
    model : :py:class:`.ExchangeHamiltonian`
        Exchange Hamiltonian.
    Q : (3,) |array_like|_
        Ordering wave vector of the spin-spiral.
        In relative coordinates with respect to the model`s reciprocal cell.
        It rotates spins from unit cell to unit cell,
        but not from atom to atom in (0, 0, 0) unit cell.
    n : (3,) |array_like|_, optional
        Global rotational axis. If None provided, then it is set to the direction of ``Q``.

    Attributes
    ----------
    Q : (3,) :numpy:`ndarray`
        Ordering wave vector of the spin-spiral. in absolute coordinates in reciprocal space.
    n : (3,) :numpy:`ndarray`
        Global rotational axis.
    N : int
        Number of magnetic atoms.
    J : (M,) :numpy:`ndarray`
        Exchange parameters.
    i : (M,) :numpy:`ndarray`
        Indices of the first atom in the exchange pair.
    j : (M,) :numpy:`ndarray`
        Indices of the second atom in the exchange pair.
    d : (M, 3) :numpy:`ndarray`
        Vectors from the first atom to the second atom in the exchange pair.
    S : (N, 3) :numpy:`ndarray`
        Spin vectors.
    u : (N, 3) :numpy:`ndarray`
        Defined from local spin directions.
    v : (N, 3) :numpy:`ndarray`
        Defined from local spin directions.
    """

    def __init__(self, model: ExchangeHamiltonian, Q=None, n=None):
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        if Q is None:
            Q = [0, 0, 0]
        self.Q = np.array(Q, dtype=float) @ self._model.crystal.lattice.reciprocal_cell

        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
            if np.allclose([0, 0, 0], Q):
                self.n = np.array([0, 0, 1])
            else:
                self.n = Q / np.linalg.norm(Q)
        else:
            self.n = np.array(n, dtype=float) / np.linalg.norm(n)

        # Get the number of magnetic atoms
        self.N = len(self._model.magnetic_atoms)

        # Get the exchange parameters, indices and vectors form the ExchangeHamiltonian
        self.J, self.i, self.j, self.d = self._model.input_for_magnons()

        # Initialize spin vector, u vector and v vector arrays
        self.S = np.zeros((self.N, 3), dtype=float)
        self.u = np.zeros((self.N, 3), dtype=complex)
        self.v = np.zeros((self.N, 3), dtype=complex)

        # Get spin vectors and compute u and v vectors from local spin directions
        for a_i, atom in enumerate(self._model.magnetic_atoms):
            try:
                self.S[a_i] = atom.spin_vector
                e1, e2, e3 = span_orthonormal_set(self.S[a_i])
                self.v[a_i] = e3
                self.u[a_i] = e1 + 1j * e2
            except ValueError:
                raise ValueError(
                    f"Spin vector is not defined for {atom.fullname} atom."
                )

        # Rotate exchange matrices
        for i in range(len(self.J)):
            rotvec = self.n * (Q @ self.d[i])
            R_nm = Rotation.from_rotvec(rotvec).as_matrix()
            self.J[i] = self.J[i] @ R_nm

    def omega(self, k, return_none=False):
        r"""
        Computes magnon energies.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        return_none : bool, default=False
            If True, then return None if Colpa fails.

        Returns
        -------
        omegas : (N,) :numpy:`ndarray`
            Magnon energies for the vector ``k``.
        """

        # Initialize matrices
        A_1 = np.zeros((self.N, self.N), dtype=complex)
        A_2 = np.zeros((self.N, self.N), dtype=complex)
        B = np.zeros((self.N, self.N), dtype=complex)
        C = np.zeros((self.N, self.N), dtype=complex)

        # Compute A(k) (A_1), A(-k)^* (A_2), B and C matrices
        for index in range(len(self.J)):
            i = self.i[index]
            j = self.j[index]
            A_1[i][j] += (
                sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                / 2
                * (self.u[i] @ self.J[index] @ np.conjugate(self.u[j]))
                * np.exp(1j * (k @ self.d[index]))
            )
            A_2[i][j] += (
                sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                / 2
                * (self.u[i] @ self.J[index] @ np.conjugate(self.u[j]))
                * np.exp(-1j * (k @ self.d[index]))
            )
            B[i][j] += (
                sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                / 2
                * (self.u[i] @ self.J[index] @ self.u[j])
                * np.exp(1j * (k @ self.d[index]))
            )
            C[i, i] += np.linalg.norm(self.S[j]) * (
                self.v[i] @ self.J[index] @ self.v[j]
            )

        # Compute h matrix
        left = np.concatenate((2 * A_1 - 2 * C, 2 * np.conjugate(B).T), axis=0)
        right = np.concatenate((2 * B, 2 * np.conjugate(A_2) - 2 * C), axis=0)
        h = np.concatenate((left, right), axis=1)

        # Diagonalize h matrix via Colpa
        try:
            omegas, U = solve_via_colpa(h)
            omegas = omegas.real[: self.N]
        except ColpaFailed:
            if return_none:
                omegas = np.array([None] * self.N)
            else:
                omegas = np.zeros(self.N)

        return omegas
