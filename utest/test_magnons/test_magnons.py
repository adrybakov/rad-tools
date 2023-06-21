import pytest

import numpy as np

from radtools.magnons.magnons import MagnonDispersion
from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.crystal.bravais_lattice import lattice_example
from radtools.crystal.crystal import Crystal
from radtools.crystal.atom import Atom
from radtools.routines import print_2d_array
from math import cos


def test_ferromagnetic():
    model = ExchangeHamiltonian()
    model.crystal.lattice = lattice_example("cub")
    model.add_atom(Atom("Fe", (0, 0, 0), spin=[0, 0, 1]))
    J = ExchangeParameter(iso=1)
    model.add_bond(J, "Fe", "Fe", (1, 0, 0))
    model.add_bond(J, "Fe", "Fe", (0, 1, 0))
    model.add_bond(J, "Fe", "Fe", (0, 0, 1))
    model.notation = "standard"
    print()
    for a1, a2, R, J in model:
        print(a1, a2, R)
        print(J)

    dispersion = MagnonDispersion(model)
    print()
    for i in dispersion.J:
        print(i)
    print()
    for i in dispersion.S:
        print(i)
    print()
    for i in dispersion.u:
        print(i)
    print()
    for i in dispersion.v:
        print(i)
    print()
    for i in dispersion.d:
        print(i)
    kp = model.get_kpoints()
    computed_omegas = []
    analytical_omegas = []
    for point in kp.points:
        computed_omegas.append(dispersion.omega(point)[0])
        analytical_omegas.append(
            -12
            * (
                (
                    cos(point[0] * model.a)
                    + cos(point[1] * model.b)
                    + cos(point[2] * model.c)
                )
                / 3
                - 1
            )
        )
    test = np.array(
        [
            kp.points[:, 0],
            kp.points[:, 1],
            kp.points[:, 2],
            computed_omegas,
            analytical_omegas,
        ]
    ).T

    assert np.allclose(computed_omegas, analytical_omegas)

    print_2d_array(test)

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(kp.flatten_points, analytical_omegas, label="analytical", alpha=0.7)
    ax.plot(
        kp.flatten_points,
        computed_omegas,
        label="computed",
        alpha=0.5,
        linestyle="dashed",
        lw=3,
    )
    ax.plot(
        kp.flatten_points,
        0.5 * test[:, 3],
        label="SpinW",
        alpha=0.5,
        linestyle="dashed",
    )

    ax.set_xticks(kp.coordinates, kp.labels)
    ax.vlines(
        kp.coordinates,
        0,
        1,
        transform=ax.get_xaxis_transform(),
        linestyles="dashed",
        color="black",
        lw=0.5,
        alpha=0.5,
    )
    ax.hlines(
        [8, 16, 24],
        0,
        1,
        transform=ax.get_yaxis_transform(),
        linestyles="dashed",
        color="black",
        lw=0.5,
        alpha=0.5,
    )
    ax.set_ylim(0, None)
    ax.set_xlim(kp.coordinates[0], kp.coordinates[-1])

    ax.legend()
    plt.savefig("test.png", dpi=400, bbox_inches="tight")
