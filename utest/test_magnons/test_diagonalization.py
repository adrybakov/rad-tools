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

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa


@pytest.mark.parametrize(
    "D",
    [
        ([[1, 0], [0, 1]]),
        (np.eye(20)),
        (np.diag(np.concatenate((np.linspace(1, 10, 10), np.linspace(10, 1, 10))))),
        ([[2, 1], [1, 2]]),
        ([[2, 0.1, 1, 0], [0.1, 2, 0, 1], [1, 0, 2, 0.1], [0, 1, 0.1, 2]]),
    ],
)
def test_solve_via_colpa(D):
    N = len(D) // 2
    E, G = solve_via_colpa(D)
    assert np.allclose(
        np.diag(E), np.linalg.inv(np.conjugate(G).T) @ D @ np.linalg.inv(G), rtol=1e-5
    )
    assert np.allclose(E[:N], E[: N - 1 : -1], rtol=1e-5)
    assert (E > 0).all()
    assert np.allclose(
        np.sum(G[:N, :N] ** 2, axis=1) - np.sum(G[:N, N:] ** 2, axis=1),
        np.ones(N),
        rtol=1e-5,
    )
    assert np.allclose(
        np.sum(G[N:, :N] ** 2, axis=1) - np.sum(G[N:, N:] ** 2, axis=1),
        -np.ones(N),
        rtol=1e-5,
    )


@pytest.mark.parametrize("D", [[[1, 0], [0, 0]], [[1, -1], [-1, -1]]])
def test_fail_via_colpa(D):
    with pytest.raises(ColpaFailed):
        solve_via_colpa(D)
