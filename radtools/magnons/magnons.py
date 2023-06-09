from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.io.tb2j import read_tb2j_model
from radtools.routines import print_2D_array
from radtools.crystal.atom import Atom

from scipy.spatial.transform import Rotation
import numpy as np
from math import sqrt, pi

import matplotlib.pyplot as plt

verbose = False


class MagnonDispersion:
    def __init__():
        pass


def input_from_ham(model: ExchangeHamiltonian):
    r"""
    Input from Exchange model.
    """


def prepare_dicts(model: ExchangeHamiltonian, Q, n, S, k):
    r"""

    Parameters
    ----------
    model : :py:class:`.ExchangeHamiltonian`
    Q : array
        Q-vector of incommensurate rotation.
        In relative coordinates with respect to the reciprocal lattice vectors.
    n : array
        Global rotation axis.
    S : list
        List of the spin directions inside the unit cell.
    k : array
        k-points
        In relative coordinates with respect to the reciprocal lattice vectors.
    """

    Q = np.array(Q)
    Q = Q[0] * model.b1 + Q[1] * model.b2 + Q[2] * model.b3
    k = k[0] * model.b1 + k[1] * model.b2 + k[2] * model.b3

    n = np.array(n)

    for i in range(0, len(S)):
        S[i] = np.array(S[i])

    J = {}
    d = {}
    indexing = {}
    magnetic_atoms = [atom for atom in model.magnetic_atoms]

    for i, atom in enumerate(magnetic_atoms):
        indexing[atom] = i

    J = []
    d = []
    ij = []

    for atom1, atom2, R, param in model:
        i, j = indexing[atom1], indexing[atom2]
        J.append(np.array(param.matrix, dtype=complex))
        d.append(np.array(model.get_bond_vector(atom1, atom1, R)))
        ij.append((i, j))

    for y, (i, j) in enumerate(ij):
        rotvec = -n * np.linalg.norm(Q * d[y])
        R = Rotation.from_rotvec(rotvec).as_matrix()
        J[y] = np.matmul(J[y], R)

    if len(S) != model.number_spins_in_unit_cell:
        raise ValueError("Different amount of spins and exchanges.")
    N = len(S)

    u = np.zeros((N, 3), dtype=complex)
    v = np.zeros((N, 3), dtype=complex)
    for i in range(0, N):
        angle = np.linalg.norm(
            Q * model.get_distance(magnetic_atoms[0], magnetic_atoms[i], (0, 0, 0))
        )
        rotvec = n * angle
        R = Rotation.from_rotvec(rotvec).as_matrix()
        v[i] = R.T[2]
        u[i] = R.T[0] + 1j * R.T[1]

    # All good
    A_1 = np.zeros((N, N), dtype=complex)
    A_2 = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    C = np.zeros((N, N), dtype=complex)
    for y, (i, j) in enumerate(ij):
        if verbose:
            print(20 * "=")
            print(i, j)
            print(np.matmul(u[i], np.matmul(J[y], np.conjugate(u[j]))))
            print(np.matmul(k, d[y]))
            print_2D_array([k, d[y]])
            print_2D_array([[np.exp(-1j * np.matmul(k, d[y]))]])
            print_2D_array([[np.exp(1j * np.matmul(k, d[y]))]])
            print(20 * "=")
        A_1[i][j] += (
            sqrt(np.linalg.norm(S[i]) * np.linalg.norm(S[j]))  # good
            / 2
            * np.matmul(u[i], np.matmul(J[y], np.conjugate(u[j])))
            * np.exp(1j * np.matmul(k, d[y]))
        )
        A_2[i][j] += (
            sqrt(np.linalg.norm(S[i]) * np.linalg.norm(S[j]))  # good
            / 2
            * np.matmul(u[i], np.matmul(J[y], np.conjugate(u[j])))
            * np.exp(-1j * np.matmul(k, d[y]))
        )
        B[i][j] += (
            sqrt(np.linalg.norm(S[i]) * np.linalg.norm(S[j]))
            / 2
            * np.matmul(u[i], np.matmul(J[y], u[j]))
            * np.exp(1j * np.matmul(k, d[y]))
        )
        C[i, i] += np.linalg.norm(S[j]) * np.matmul(v[i], np.matmul(J[y], v[j]))

    left = np.concatenate((2 * A_1 - 2 * C, 2 * np.conjugate(B).T), axis=0)
    right = np.concatenate((2 * B, 2 * np.conjugate(A_2) - 2 * C), axis=0)
    h = -np.concatenate((left, right), axis=1)
    if verbose:
        print("A_1")
        print_2D_array(A_1)
        print("A_2")
        print_2D_array(A_2)
        print("B")
        print_2D_array(B)
        print("C")
        print_2D_array(C)
        print("h")
        print_2D_array(h)
        print(np.linalg.eigvals(h))
    # return A_1.flatten()
    try:
        K = np.linalg.cholesky(h)
        print(f"{k} You are good")
    except:
        print(f"{k[0]:5.2f} {k[1]:5.2f} {k[2]:5.2f} Fuck you")
        return np.zeros(2 * N)

    g = np.eye(2 * N)
    for i in range(N, 2 * N):
        g[i][i] = -1

    omegas, U = np.linalg.eig(np.matmul(np.matmul(K, g), np.conjugate(K).T))
    omegas = list(omegas.real)
    omegas.sort(key=lambda x: -x)
    return omegas


if __name__ == "__main__":
    model = ExchangeHamiltonian()
    model.add_atom(Atom("Fe", (0, 0, 0)))
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
    # model = read_exchange_model(
    #     "/Volumes/work-backup/projects (closed)/2022 Nanoletters CrSBr/Calculations/UniaxialA/100.3/TB2J_results/exchange.out"
    # )
    model = read_tb2j_model(
        "/Users/rybakov.ad/Projects/rad-tools/debug/niI2_v2/exchange.out"
    )

    # print(model.magnetic_atoms)
    kpoints = np.linspace(0, 0.5, 100)
    omegas = []
    S = [0, 0, 1]

    def plot_graph(kpoints, omegas, model, S, Q, n):
        for i in kpoints:
            omegas.append(prepare_dicts(model, Q, n, S, [i, 0, 0]))
        for i in kpoints:
            omegas.append(
                prepare_dicts(
                    model,
                    Q,
                    n,
                    [S],
                    [0.5, i, 0],
                )
            )
        for i in kpoints:
            omegas.append(
                prepare_dicts(
                    model,
                    Q,
                    n,
                    [S],
                    [(0.5 - i), 0.5, 0],
                )
            )
        for i in kpoints:
            omegas.append(
                prepare_dicts(
                    model,
                    Q,
                    n,
                    [S],
                    [0, (0.5 - i), 0],
                )
            )
        omegas = np.array(omegas).T

        fig, ax = plt.subplots()
        kpoints = np.concatenate((kpoints, 0.5 + kpoints, 1 + kpoints, 1.5 + kpoints))

        for i in range(len(omegas)):
            ax.plot(kpoints, omegas[i].real, label=f"{i+1}")
        ax.set_xlabel(R"k, $0 \rightarrow \pi$")
        ax.set_ylabel("E, meV")
        ax.legend()
        plt.savefig("test.png", dpi=400, bbox_inches="tight")

    plot_graph(kpoints, omegas, model, [S], [0, 0, 0], [0, 0, 1])
