from copy import deepcopy
from math import sqrt

import numpy as np
from scipy.spatial.transform import Rotation

from radtools.crystal.kpoints import Kpoints
from radtools.geometry import span_orthonormal_set
from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa
from radtools.spinham.hamiltonian import SpinHamiltonian

__all__ = ["MagnonDispersion"]


class MagnonDispersion:
    r"""
    Magnon dispersion wrapper.

    Parameters
    ----------
    model : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.
    Q : (3,) |array_like|_
        Ordering wave vector of the spin-spiral.
        In relative coordinates with respect to the model`s reciprocal cell.
        It rotates spins from unit cell to unit cell,
        but not from atom to atom in (0, 0, 0) unit cell.
    n : (3,) |array_like|_, optional
        Global rotational axis. If None provided, then it is set to the direction of ``Q``.
    nodmi : bool, default=False
        If True, then DMI is not included in the dispersion.
    noaniso : bool, default=False
        If True, then anisotropy is not included in the dispersion.
    custom_mask : func
        Custom mask for the exchange parameter. Function which take (3,3) :numpy:`ndarray`
        as an input and returns (3,3) :numpy:`ndarray` as an output.

    Attributes
    ----------
    Q : (3,) :numpy:`ndarray`
        Ordering wave vector of the spin-spiral. in absolute coordinates in reciprocal space.
    n : (3,) :numpy:`ndarray`
        Global rotational axis.
    N : int
        Number of magnetic atoms.
    J_matrices : (M,) :numpy:`ndarray`
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

    def __init__(
        self,
        model: SpinHamiltonian,
        Q=None,
        n=None,
        nodmi=False,
        noaniso=False,
        custom_mask=None,
    ):
        self._C = None
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        if Q is None:
            Q = [0, 0, 0]
        self.Q = np.array(Q, dtype=float) @ self._model.reciprocal_cell

        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
            if np.allclose([0, 0, 0], Q):
                self.n = np.array([0, 0, 1])
            else:
                self.n = self.Q / np.linalg.norm(self.Q)
        else:
            self.n = np.array(n, dtype=float) / np.linalg.norm(n)

        # Get the number of magnetic atoms
        self.N = len(self._model.magnetic_atoms)

        # Get the exchange parameters, indices and vectors form the SpinHamiltonian
        (
            self.J_matrices,
            self.indices_i,
            self.indices_j,
            self.dis_vectors,
        ) = self._model.input_for_magnons(
            nodmi=nodmi, noaniso=noaniso, custom_mask=custom_mask
        )

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
        for i in range(len(self.J_matrices)):
            rotvec = self.n * (self.Q @ self.dis_vectors[i])
            R_nm = Rotation.from_rotvec(rotvec).as_matrix()
            self.J_matrices[i] = self.J_matrices[i] @ R_nm

    def J(self, k):
        r"""
        Computes J(k) matrix.

        .. math::

            \boldsymbol{J}_{i,j}(\boldsymbol{k}) = \sum_{\boldsymbol{d}}\boldsymbol{J}_{i,j}(\boldsymbol{d})e^{-i\boldsymbol{k}\boldsymbol{d}}

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N, 3, 3), dtype=complex)
        k = np.array(k)
        # Compute J(k)
        for index in range(len(self.J_matrices)):
            i = self.indices_i[index]
            j = self.indices_j[index]

            result[i][j] += self.J_matrices[index] * np.exp(
                -1j * (k @ self.dis_vectors[index])
            )
        return result

    def A(self, k):
        r"""
        Computes A(k) matrix.

        .. math::

            A(\boldsymbol{k})^{ij} = \dfrac{\sqrt{S_i\cdot S_j}}{2}\boldsymbol{u}^T_i\boldsymbol{J}_{i,j}(-\boldsymbol{k})\overline{\boldsymbol{u}}_j

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N), dtype=complex)
        k = np.array(k)

        # Compute A(k)
        J = self.J(-k)
        for i in range(len(J)):
            for j in range(len(J[i])):
                result[i][j] += (
                    sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                    / 2
                    * (self.u[i] @ J[i][j] @ np.conjugate(self.u[j]))
                )
        return result

    def B(self, k):
        r"""
        Computes B(k) matrix.

        .. math::

            B(\boldsymbol{k})^{ij} = \dfrac{\sqrt{S_i\cdot S_j}}{2}\boldsymbol{u}^T_i\boldsymbol{J}_{i,j}(-\boldsymbol{k})\boldsymbol{u}_j

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N), dtype=complex)
        k = np.array(k)

        # Compute B(k)
        J = self.J(-k)
        for i in range(len(J)):
            for j in range(len(J[i])):
                result[i][j] += (
                    sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                    / 2
                    * (self.u[i] @ J[i][j] @ self.u[j])
                )
        return result

    def C(self):
        r"""
        Computes C matrix.

        .. math::

            C^{i,j} = \delta_{i,j}\sum_{l}S_l \boldsymbol{v}^T_i\boldsymbol{J}_{i, l}(\boldsymbol{0})\boldsymbol{v}_l

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.
        """

        if self._C is None:
            self._C = np.zeros((self.N, self.N), dtype=complex)

            # Compute C matrix, note: sum over l is hidden here
            J = self.J(np.zeros(3))
            for i in range(len(J)):
                for j in range(len(J[i])):
                    self._C[i][i] += (
                        np.linalg.norm(self.S[j]) * self.v[i] @ J[i][j] @ self.v[j]
                    )
        return self._C

    def h(self, k):
        k = np.array(k)
        # Compute h matrix
        left = np.concatenate(
            (2 * self.A(k) - 2 * self.C(), 2 * np.conjugate(self.B(k)).T), axis=0
        )
        right = np.concatenate(
            (2 * self.B(k), 2 * np.conjugate(self.A(-k)) - 2 * self.C()), axis=0
        )
        h = np.concatenate((left, right), axis=1)
        return h

    def omega(self, k, zeros_to_none=False):
        r"""
        Computes magnon energies.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        zeros_to_none : bool, default=False
            If True, then return ``None`` instead of 0 if Colpa fails.

        Returns
        -------
        omegas : (N,) :numpy:`ndarray`
            Magnon energies for the vector ``k``.
        """
        # Diagonalize h matrix via Colpa
        k = np.array(k)
        h = self.h(k)
        try:
            omegas, U = solve_via_colpa(h)
            omegas = omegas.real[: self.N]
        except ColpaFailed:
            # Try to solve for positive semidefinite matrix
            try:
                omegas, U = solve_via_colpa(h + np.diag(1e-8 * np.ones(2 * self.N)))
                omegas = omegas.real[: self.N]
            except ColpaFailed:
                # Try to solve for negative defined matrix
                try:
                    omegas, U = solve_via_colpa(-h)
                    omegas = omegas.real[: self.N] * -1
                except ColpaFailed:
                    # Try to solve for negative semidefinite matrix
                    try:
                        omegas, U = solve_via_colpa(
                            -h - np.diag(1e-8 * np.ones(h.shape[0]))
                        )
                        omegas = omegas.real[: self.N] * -1
                    except ColpaFailed:
                        # If all fails, return None or 0
                        if zeros_to_none:
                            omegas = np.array([None] * self.N)
                        else:
                            omegas = np.zeros(self.N)
        omegas[np.abs(omegas) <= 1e-8] = 0
        return omegas

    def omegas(self, kpoints, zeros_to_none=False):
        r"""
        Dispersion spectra.

        Parameters
        ----------
        kpoints : (M, 3) |array_like|_
            K points in absolute coordinates.
        zeros_to_none : bool, default=False
            If True, then return ``None`` instead of 0 if Colpa fails.
        """

        data = []

        if isinstance(kpoints, Kpoints):
            kpoints = kpoints.points()

        for point in kpoints:
            data.append(self.omega(point, zeros_to_none=zeros_to_none))

        return np.array(data).T

    def __call__(self, *args, **kwargs):
        return self.omegas(*args, **kwargs)
