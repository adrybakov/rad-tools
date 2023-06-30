from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.routines import span_orthonormal_set
from radtools.crystal.atom import Atom
from radtools.magnons.diagonalization import solve_via_colpa, ColpaFailed
from radtools.routines import print_2d_array, plot_horizontal_lines

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
        self._C = None
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        if Q is None:
            Q = [0, 0, 0]
        self.Q = np.array(Q, dtype=float) @ self._model.crystal.lattice.reciprocal_cell

        print("Q: ", self.Q)

        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
            if np.allclose([0, 0, 0], Q):
                self.n = np.array([0, 0, 1])
            else:
                self.n = Q / np.linalg.norm(Q)
        else:
            self.n = np.array(n, dtype=float) / np.linalg.norm(n)

        print("n: ", self.n)

        # Get the number of magnetic atoms
        self.N = len(self._model.magnetic_atoms)
        print("N: ", self.N)

        # Get the exchange parameters, indices and vectors form the ExchangeHamiltonian
        (
            self.J_matrices,
            self.indices_i,
            self.indices_j,
            self.dis_vectors,
        ) = self._model.input_for_magnons()

        for i in range(0, len(self.J_matrices)):
            print(self.indices_i[i], self.indices_j[i], self.dis_vectors[i])
            print_2d_array(self.J_matrices[i], fmt=".4f")

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
            rotvec = self.n * (Q @ self.dis_vectors[i])
            R_nm = Rotation.from_rotvec(rotvec).as_matrix()
            print_2d_array(R_nm)
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
        result = np.zeros((3, 3), dtype=complex)

        # Compute J(k)
        for index in range(len(self.J_matrices)):
            result += self.J_matrices[index] * np.exp(
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

        # Compute A(k)
        for index in range(len(self.J)):
            i = self.indices_i[index]
            j = self.indices_j[index]
            result[i][j] += (
                sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                / 2
                * (self.u[i] @ self.J(-k) @ np.conjugate(self.u[j]))
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

        # Compute B(k)
        for index in range(len(self.J_matrices)):
            i = self.indices_i[index]
            j = self.indices_j[index]
            result[i][j] += (
                sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                / 2
                * (self.u[i] @ self.J(-k) @ self.u[j])
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
            for index in range(len(self.J_matrices)):
                i = self.indices_i[index]
                j = self.indices_j[index]
                self._C[i, i] += (
                    np.linalg.norm(self.S[j]) * self.v[i] @ self.J(0) @ self.v[j]
                )
        return self._C

    def h(self, k):
        # Compute h matrix
        left = np.concatenate(
            (2 * self.A(k) - 2 * self.C(), 2 * np.conjugate(self.B(k)).T), axis=0
        )
        right = np.concatenate(
            (2 * self.B(k), 2 * np.conjugate(self.A(-k)) - 2 * self.C()), axis=0
        )
        h = np.concatenate((left, right), axis=1)
        print(k)
        print_2d_array(h)
        return h

    # TODO Create compute() method and move all dispersion-related storage in the class. Make it iterable over omegas (not k points, but eigenvalues)
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
        try:
            omegas, U = solve_via_colpa(self.h(k))
            omegas = omegas.real[: self.N]
        except ColpaFailed:
            if zeros_to_none:
                omegas = np.array([None] * self.N)
            else:
                omegas = np.zeros(self.N)

        return omegas


if __name__ == "__main__":
    from radtools.io.tb2j import read_tb2j_model
    from termcolor import cprint

    model = read_tb2j_model(
        "/Users/rybakov.ad/Desktop/exchange.out",
        quiet=False,
        bravais_type="HEX",
    )
    model.filter(max_distance=8)
    cprint(f"{model.variation} crystal detected", "green")
    cprint(f"Notation is {model.notation}", "green")
    model.crystal.lattice.add_kpoint("Mprime", [0, 0.5, 0], "M$^{\prime}$")
    model.crystal.lattice.add_kpoint("Kprime", [-1 / 3, 2 / 3, 0], "K$^{\prime}$")
    kp = model.crystal.lattice.get_kpoints()  # Set custom k path
    # kp.path = "G-K-M-G"
    kp.path = "Mprime-G-M-K-Mprime-Kprime-G-K"
    kp.n = 40

    spin = ["Ni1", 0, 1, 0]
    if spin is not None:
        for i in range(len(spin) // 4):
            atom_name = spin[4 * i]
            atom = model.crystal.get_atom(atom_name)
            atom_spin = list(map(float, spin[4 * i + 1 : 4 * i + 4]))
            atom.spin_vector = atom_spin

    # Get the magnon dispersion
    dispersion = MagnonDispersion(model, Q=(0.138, 0, 0))
    dispersion2 = MagnonDispersion(model)

    fig, ax = plt.subplots()

    omegas = []
    omegas_fm = []
    for point in kp.points:
        omegas.append(dispersion.omega(point, zeros_to_none=0))
        omegas_fm.append(dispersion2.omega(point, zeros_to_none=0))

    omegas = np.array(omegas).T
    omegas_fm = np.array(omegas_fm).T

    ax.set_xticks(kp.coordinates, kp.labels, fontsize=15)
    ax.set_ylabel("E, meV", fontsize=15)
    ax.vlines(
        kp.coordinates,
        0,
        1,
        transform=ax.get_xaxis_transform(),
        colors="black",
        alpha=0.5,
        lw=1,
        ls="dashed",
    )
    ax.plot(kp.flatten_points, omegas[0], label="helix")
    ax.plot(kp.flatten_points, omegas_fm[0], label="fm")
    ax.legend()

    ax.set_xlim(kp.flatten_points[0], kp.flatten_points[-1])
    ax.set_ylim(-1, None)
    plot_horizontal_lines(ax, 0)

    plt.savefig(
        f"magnon_dispersion.png",
        bbox_inches="tight",
        dpi=600,
    )
