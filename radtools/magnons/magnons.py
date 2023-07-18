from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.routines import span_orthonormal_set, absolute_to_relative
from radtools.crystal.atom import Atom
from radtools.magnons.diagonalization import solve_via_colpa, ColpaFailed
from radtools.routines import print_2d_array, plot_horizontal_lines, plot_vertical_lines

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
    nodmi : bool, default=False
        If True, then DMI is not included in the dispersion.
    noaniso : bool, default=False
        If True, then anisotropy is not included in the dispersion.
    custom_mask : func
        Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
        as an input and returns (3,3) numpy:`ndarray` as an output.

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

    def __init__(
        self,
        model: ExchangeHamiltonian,
        Q=None,
        n=None,
        nodmi=False,
        noaniso=False,
        custom_mask=None,
    ):
        self._omegas = None
        self._C = None
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        if Q is None:
            Q = [0, 0, 0]
        self.Q = np.array(Q, dtype=float) @ self._model.reciprocal_cell
        print("Q:", self.Q)
        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
            if np.allclose([0, 0, 0], Q):
                self.n = np.array([0, 0, 1])
            else:
                self.n = self.Q / np.linalg.norm(self.Q)
        else:
            self.n = np.array(n, dtype=float) / np.linalg.norm(n)
        print("n:", self.n)
        # Get the number of magnetic atoms
        self.N = len(self._model.magnetic_atoms)

        # Get the exchange parameters, indices and vectors form the ExchangeHamiltonian
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

        for index in range(len(self.S)):
            print(f"S ({index}) : {self.S[index]}")
            print(f"u ({index}) : {self.u[index]}")
            print(f"v ({index}) : {self.v[index]}")

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

            #     print(
            #         i,
            #         j,
            #         k,
            #         f"{k @ self.dis_vectors[index]:.4f}",
            #         self.dis_vectors[index],
            #     )
            #     print(result[i][j], np.exp(-1j * (k @ self.dis_vectors[index])))
            result[i][j] += self.J_matrices[index] * np.exp(
                -1j * (k @ self.dis_vectors[index])
            )
        # print()
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
        try:
            omegas, U = solve_via_colpa(self.h(k))
            omegas = omegas.real[: self.N]
        except ColpaFailed:
            try:
                omegas, U = solve_via_colpa(-self.h(k))
                omegas = omegas.real[: self.N] * -1
            except ColpaFailed:
                if zeros_to_none:
                    omegas = np.array([None] * self.N)
                else:
                    omegas = np.zeros(self.N)
            # if zeros_to_none:
            #     omegas = np.array([None] * self.N)
            # else:
            #     omegas = np.zeros(self.N)

        return omegas

    def omegas(self):
        r"""
        Dispersion spectra.
        """

        if self._omegas is None:
            self.compute()

        return self._omegas

    def compute(self, kpoints, zeros_to_none=False):
        r"""
        Computes magnon energies with respect to the class settings.

        Parameters
        ----------
        zeros_to_none : bool, default=False
            If True, then return ``None`` instead of 0 if Colpa fails.
        """

        self._omegas = []

        for point in kpoints:
            self._omegas.append(self.omega(point, zeros_to_none=zeros_to_none))

        self._omegas = np.array(self._omegas).T


if __name__ == "__main__":
    # from radtools.io.tb2j import read_tb2j_model

    # model = read_tb2j_model("/Users/rybakov.ad/Desktop/exchange.out")
    # model.get_atom("Ni1").spin_vector = [0, 0, 1]
    # model.get_atom("Ni2").spin_vector = [0, 0, -1]
    # model.filter(max_distance=5)
    # for a1, a2, R, J in model:
    #     print(R, model.get_distance(a1, a2, R))

    # dispersion = MagnonDispersion(model)

    # for dis in dispersion.dis_vectors:
    #     print(dis)
    # for m in dispersion.J_matrices:
    #     print_2d_array(m, ".4f")
    # KPOINT = np.array([0.5, 0, 0]) @ model.reciprocal_cell
    # print(dispersion.J(KPOINT).shape)
    # for i in [0, 1]:
    #     for j in [0, 1]:
    #         print(f"{i} {j}")
    #         print_2d_array(dispersion.J(KPOINT)[i][j], ".4f")

    # print("h")
    # print_2d_array(dispersion.h(KPOINT), ".4f")

    # print("A")
    # print_2d_array(dispersion.A(KPOINT), ".4f")

    # print("B")
    # print_2d_array(dispersion.B(KPOINT), ".4f")

    # print("C")
    # print_2d_array(dispersion.C(), ".4f")

    # print("omega")
    # print(dispersion.omega(KPOINT))

    from radtools.io.tb2j import read_tb2j_model
    from termcolor import cprint

    model = read_tb2j_model(
        "debug/magnons/exchange.out",
        bravais_type="HEX",
    )
    model.filter(max_distance=8)
    cprint(f"{model.variation} crystal detected", "green")
    cprint(f"Notation is {model.notation}", "green")
    model.kpoints.add_hs_point("Mprime", [0, 0.5, 0], "M$^{\prime}$")
    model.kpoints.add_hs_point("Kprime", [-1 / 3, 2 / 3, 0], "K$^{\prime}$")
    kp = model.kpoints  # Set custom k path
    kp.path = "G-M-K-G"
    # kp.path = "Mprime-G-M-K-Mprime-Kprime-G-K"
    kp.n = 100

    spin = ["Ni1", 0, 0, 1]
    if spin is not None:
        for i in range(len(spin) // 4):
            atom_name = spin[4 * i]
            atom = model.get_atom(atom_name)
            atom_spin = list(map(float, spin[4 * i + 1 : 4 * i + 4]))
            print("here")
            atom.spin_vector = atom_spin

    # # Get the magnon dispersion
    dispersion = MagnonDispersion(model, Q=(0.12568743, 0.12568743, 0.0), n=[0, 1, 0])
    dispersion2 = MagnonDispersion(model)
    KPOINT = [0.28282347, 0.16325446, 0.0]
    print(dispersion2.omega(KPOINT))
    print("h")
    print_2d_array(dispersion2.h(KPOINT), ".4f")
    print("A")
    print_2d_array(dispersion2.A(KPOINT), ".4f")
    print("B")
    print_2d_array(dispersion2.B(KPOINT), ".4f")
    print("C")
    print_2d_array(dispersion2.C(), ".4f")
    print("J")
    print_2d_array(dispersion2.J(KPOINT)[0][0], ".4f")

    dispersion.compute(kp.points())
    # A = []
    # B = []
    # C = []
    # h = []

    # for point in kp.points():
    #     A.append(dispersion.A(point))
    #     B.append(dispersion.B(point))
    #     C.append(dispersion.C())
    #     h.append(dispersion.h(point))
    # h = np.array(h)

    # fig, ax = plt.subplots(15, 1, figsize=(5, 10))

    # fig.subplots_adjust(hspace=0)
    # ax[0].plot(kp.flatten_points(), np.array(A).real[:, 0, 0], label="A real")
    # ax[1].plot(kp.flatten_points(), np.array(A).imag[:, 0, 0], label="A imag")
    # ax[2].plot(kp.flatten_points(), np.array(B).real[:, 0, 0], label="B real")
    # ax[3].plot(kp.flatten_points(), np.array(B).imag[:, 0, 0], label="B imag")
    # ax[4].plot(kp.flatten_points(), np.array(C).real[:, 0, 0], label="C real")
    # ax[5].plot(kp.flatten_points(), np.array(C).imag[:, 0, 0], label="C imag")
    # ax[6].plot(kp.flatten_points(), h.real[:, 0, 0], color="red", label="h 0 0 real")
    # ax[7].plot(kp.flatten_points(), h.real[:, 0, 1], color="green", label="h 0 1 real")
    # ax[8].plot(kp.flatten_points(), h.real[:, 1, 0], color="black", label="h 1 0 real")
    # ax[9].plot(
    #     kp.flatten_points(), h.real[:, 1, 1], color="magenta", label="h 1 1 real"
    # )
    # ax[10].plot(kp.flatten_points(), h.imag[:, 0, 0], color="red", label="h 0 0 imag")
    # ax[11].plot(kp.flatten_points(), h.imag[:, 0, 1], color="green", label="h 0 1 imag")
    # ax[12].plot(kp.flatten_points(), h.imag[:, 1, 0], color="black", label="h 1 0 imag")
    # ax[13].plot(
    #     kp.flatten_points(), h.imag[:, 1, 1], color="magenta", label="h 1 1 imag"
    # )
    # ax[14].plot(
    #     kp.flatten_points(),
    #     np.array(A).real[:, 0, 0] - np.array(B).real[:, 0, 0],
    #     label="A - B real",
    # )
    # ax[0].set_ylabel("A real")
    # ax[1].set_ylabel("A imag")
    # ax[2].set_ylabel("B real")
    # ax[3].set_ylabel("B imag")
    # ax[4].set_ylabel("C real")
    # ax[5].set_ylabel("C imag")
    # ax[6].set_ylabel()
    # ax[7].set_ylabel()
    # ax[8].set_ylabel()
    # ax[9].set_ylabel()
    # ax[10].set_ylabel()
    # ax[11].set_ylabel()
    # ax[12].set_ylabel()
    # ax[13].set_ylabel()
    # for num, i in enumerate(ax):
    #     i.set_xlim(kp.coordinates()[0], kp.coordinates()[-1])
    #     i.set_xticks(kp.coordinates(), kp.labels, fontsize=15)
    #     plot_vertical_lines(i, kp.coordinates())
    #     plot_horizontal_lines(i, 0)
    #     i.legend(fontsize=8, loc="upper right")
    #     if num != 14:
    #         i.get_xaxis().set_visible(False)
    # plt.savefig("test.png", dpi=600, bbox_inches="tight")
    # plt.close()

    fig, ax = plt.subplots()

    dispersion2.compute(kp.points())

    omega_plus = []
    omega_minus = []
    for point in kp.points():
        omega_plus.append(dispersion.omega(point + dispersion.Q)[0])
        omega_minus.append(dispersion.omega(point - dispersion.Q)[0])

    # omega_plus = np.array(omega_plus)
    # omega_minus = np.array(omega_minus)

    ax.set_xticks(kp.coordinates(), kp.labels, fontsize=15)
    ax.set_ylabel("E, meV", fontsize=15)
    minval = 0
    minval_i = 0
    for i in range(len(kp.points())):
        if (
            dispersion2.omegas()[0][i] is not None
            and dispersion2.omegas()[0][i] < minval
        ):
            minval = dispersion2.omegas()[0][i]
            minval_i = i
    print(kp.points()[minval_i])
    print(absolute_to_relative(model.reciprocal_cell, kp.points()[minval_i]))

    plot_vertical_lines(ax, kp.coordinates())
    ax.plot(kp.flatten_points(), dispersion.omegas()[0], label="helix")
    # ax.plot(kp.flatten_points(), omega_plus, label="helix + k")
    # ax.plot(kp.flatten_points(), omega_minus, label="helix - k")
    ax.plot(kp.flatten_points(), dispersion2.omegas()[0], label="fm")
    ax.legend()

    ax.set_xlim(kp.flatten_points()[0], kp.flatten_points()[-1])
    # ax.set_ylim(-1, None)
    plot_horizontal_lines(ax, 0)

    plt.savefig(
        f"magnon_dispersion.png",
        bbox_inches="tight",
        dpi=600,
    )


# # problem point = [0.28524077, 0.49403396, 0.]
# # for fm [0.28282347, 0.16325446, 0.]
