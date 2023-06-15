from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.routines import span_orthonormal_set
from radtools.crystal.atom import Atom
from radtools.magnons.diagonalization import solve_via_colpa, ColpaFailed
from radtools.routines import print_2D_array

from copy import deepcopy


from scipy.spatial.transform import Rotation
import numpy as np
from math import sqrt, pi

import matplotlib.pyplot as plt


class MagnonDispersion:
    r"""
    Magnon dispersion wrapper.

    Parameters
    ----------
    model : :py:class:`.ExchangeHamiltonian`
        Exchange Hamiltonian.
    magnetic atoms : list
        List of magnetic atoms. Indices from data are used on this list.
    Q : (3,) |array_like|_
        Ordering wave vector of the spin-spiral.
        In relative coordinates with respect to the model`s reciprocal cell.
        It rotates spins from unit cell to unit cell,
        but not from atom to atom in (0, 0, 0) unit cell.
    n : (3,) |array_like|_, optional
        Global rotational axis. If None provided, then it is set to the direction of ``Q``.
    """

    def __init__(self, model: ExchangeHamiltonian, Q, n=None):
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        self.Q = np.array(Q, dtype=float) @ self._model.crystal.lattice.reciprocal_cell

        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
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
            In relative coordinates with respect to the model`s reciprocal cell.
        return_none : bool, default=False
            If True, then return None if Colpa fails.

        Returns
        -------
        omegas : (N,) :numpy:`ndarray`
            Magnon energies for the vector ``k``.
        """

        # Convert k to absolute coordinates
        k = k @ self._model.crystal.lattice.reciprocal_cell

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


if __name__ == "__main__":
    model = ExchangeHamiltonian()
    model.add_atom(Atom("Fe", (0, 0, 0), spin=[0, 0, 1]))
    # model.add_atom("Fe1", 0.75, 0.25, 0)
    # model.add_atom("Fe2", 0.25, 0.75, 0)
    bond1 = ExchangeParameter(iso=1)
    model.add_bond(bond1, "Fe", "Fe", (1, 0, 0))
    model.add_bond(bond1, "Fe", "Fe", (-1, 0, 0))
    model.add_bond(bond1, "Fe", "Fe", (0, 1, 0))
    model.add_bond(bond1, "Fe", "Fe", (0, -1, 0))
    model.add_bond(bond1, "Fe", "Fe", (0, 0, 1))
    model.add_bond(bond1, "Fe", "Fe", (0, 0, -1))
    bond2 = ExchangeParameter(iso=1)
    bond3 = ExchangeParameter(iso=1, aniso=np.diag([0, 0.1, 0.5]))
    # model.add_bond(bond1, "Fe1", "Fe1", (1, 0, 0))
    # model.add_bond(bond1, "Fe1", "Fe1", (-1, 0, 0))
    # model.add_bond(bond1, "Fe2", "Fe2", (1, 0, 0))
    # model.add_bond(bond1, "Fe2", "Fe2", (-1, 0, 0))

    # model.add_bond(bond2, "Fe1", "Fe2", (0, 0, 0))
    # model.add_bond(bond2, "Fe1", "Fe2", (1, 0, 0))
    # model.add_bond(bond2, "Fe1", "Fe2", (1, -1, 0))
    # model.add_bond(bond2, "Fe1", "Fe2", (0, -1, 0))

    # model.add_bond(bond2, "Fe2", "Fe1", (0, 0, 0))
    # model.add_bond(bond2, "Fe2", "Fe1", (0, 1, 0))
    # model.add_bond(bond2, "Fe2", "Fe1", (-1, 0, 0))
    # model.add_bond(bond2, "Fe2", "Fe1", (-1, 1, 0))

    # model.add_bond(bond3, "Fe1", "Fe1", (0, 1, 0))
    # model.add_bond(bond3, "Fe1", "Fe1", (0, -1, 0))
    # model.add_bond(bond3, "Fe2", "Fe2", (0, 1, 0))
    # model.add_bond(bond3, "Fe2", "Fe2", (0, -1, 0))
    from radtools.io import read_tb2j_model

    model = read_tb2j_model(
        "/Volumes/work-backup/projects (closed)/2022 Nanoletters CrSBr/Calculations/UniaxialA/100.3/TB2J_results/exchange.out"
    )
    # model.filter(max_distance=5)
    model.crystal.get_atom("Cr1").spin_direction = [0, 0.7, 1]
    model.crystal.get_atom("Cr2").spin_direction = [0, 0.7, 1]
    model.crystal.get_atom("Cr1").spin = 3.2857 / 2
    model.crystal.get_atom("Cr2").spin = 3.2857 / 2

    # model = read_tb2j_model(
    #     "/Users/rybakov.ad/Projects/rad-tools/debug/niI2_v2/exchange.out"
    # )
    # model.crystal.get_atom("Ni1").spin_vector = [0, 0, 1]

    dispersion = MagnonDispersion(model, [0, 0, 0], [0, 0, 1])

    kpoints = np.linspace(0, 0.5, 100)
    omegas = []
    print(dispersion.omega([0, 0, 0], return_none=True))

    def plot_graph(kpoints, omegas, model, Q, n):
        for i in kpoints:
            omegas.append(dispersion.omega([i, 0, 0], return_none=True))
        for i in kpoints:
            omegas.append(dispersion.omega([0.5, i, 0], return_none=True))
        for i in kpoints:
            omegas.append(dispersion.omega([(0.5 - i), 0.5, 0], return_none=True))
        for i in kpoints:
            omegas.append(dispersion.omega([0, (0.5 - i), 0], return_none=True))
        omegas = np.array(omegas).T

        fig, ax = plt.subplots()
        kpoints = np.concatenate((kpoints, 0.5 + kpoints, 1 + kpoints, 1.5 + kpoints))

        for i in range(len(omegas)):
            ax.plot(kpoints, omegas[i].real, label=f"{i+1}")
        ax.set_xlabel(R"k, $0 \rightarrow \pi$")
        ax.set_ylabel("E, meV")
        ax.legend()
        plt.savefig("test.png", dpi=400, bbox_inches="tight")

    plot_graph(kpoints, omegas, model, [0, 0, 0], [0, 0, 1])
