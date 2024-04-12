# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from math import cos

import numpy as np
import pytest

from radtools.crystal.atom import Atom
from radtools.crystal.bravais_lattice import lattice_example
from radtools.magnons.dispersion import MagnonDispersion
from radtools.spinham.hamiltonian import SpinHamiltonian
from radtools.spinham.parameter import ExchangeParameter


def test_ferromagnetic():
    model = SpinHamiltonian(lattice=lattice_example("CUB"))
    model.add_atom(Atom("Fe", (0, 0, 0), spin=[0, 0, 1]))
    J = ExchangeParameter(iso=1)
    model.add_bond("Fe", "Fe", (1, 0, 0), J=J)
    model.add_bond("Fe", "Fe", (0, 1, 0), J=J)
    model.add_bond("Fe", "Fe", (0, 0, 1), J=J)
    model.notation = "standard"

    dispersion = MagnonDispersion(model)

    kp = model.kpoints
    computed_omegas = []
    analytical_omegas = []
    for point in kp.points():
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
