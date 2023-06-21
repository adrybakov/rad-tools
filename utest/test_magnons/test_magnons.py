import pytest

import numpy as np

from radtools.magnons.magnons import MagnonDispersion
from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter
from radtools.crystal.bravais_lattice import lattice_example
from radtools.crystal.atom import Atom
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

    dispersion = MagnonDispersion(model)

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

    assert np.allclose(computed_omegas, analytical_omegas)
