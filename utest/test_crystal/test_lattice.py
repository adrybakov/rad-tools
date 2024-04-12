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

from math import acos, pi, sqrt

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis import target
from scipy.spatial.transform import Rotation

import radtools.crystal.cell as Cell
from radtools.crystal.bravais_lattice.examples import lattice_example
from radtools.crystal.constants import (
    ABS_TOL,
    ABS_TOL_ANGLE,
    BRAVAIS_LATTICE_NAMES,
    BRAVAIS_LATTICE_VARIATIONS,
    MAX_LENGTH,
    MIN_ANGLE,
    MIN_LENGTH,
    PEARSON_SYMBOLS,
    REL_TOL,
    REL_TOL_ANGLE,
)
from radtools.crystal.lattice import Lattice
from radtools.geometry import parallelepiped_check

n_order = 5


def shuffle(cell, order):
    if order == 0:
        return [cell[2], cell[0], cell[1]]
    if order == 1:
        return [cell[1], cell[2], cell[0]]
    if order == 2:
        return cell
    if order == 3:
        return [cell[1], cell[0], cell[2]]
    if order == 4:
        return [cell[2], cell[1], cell[0]]
    if order == 5:
        return [cell[0], cell[2], cell[1]]


def rotate(cell, r1, r2, r3):
    R = Rotation.from_rotvec([r1, r2, r3]).as_matrix()
    return cell @ R.T


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_Lattice_cell_attributes(r1, r2, r3, a, b, c, alpha, beta, gamma):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(a, b, c) / max(a, b, c) > REL_TOL
        and parallelepiped_check(a, b, c, alpha, beta, gamma)
    ):
        # Prepare cell
        cell = rotate(
            Cell.from_params(a, b, c, alpha, beta, gamma),
            r1,
            r2,
            r3,
        )
        l = Lattice(cell, standardize=False)
        # TODO: fix accuracy in Cell.from_params (sin, cos)
        # assert np.allclose(
        #     sorted([l.a, l.b, l.c]), sorted([a, b, c]), rtol=REL_TOL, atol=ABS_TOL
        # )
        # assert np.allclose(
        #     sorted([l.alpha, l.beta, l.gamma]),
        #     sorted([alpha, beta, gamma]),
        #     rtol=REL_TOL_ANGLE,
        #     atol=ABS_TOL_ANGLE,
        # )
        assert np.allclose(
            [l.a, l.b, l.c], list(l.parameters[:3]), rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            [l.alpha, l.beta, l.gamma],
            list(l.parameters[3:]),
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )
        assert np.allclose(
            l.unit_cell_volume,
            abs(np.linalg.det(l.cell)),
            rtol=REL_TOL,
            atol=ABS_TOL,
        )


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_Lattice_reciprocal_cell_attributes(r1, r2, r3, a, b, c, alpha, beta, gamma):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(a, b, c) / max(a, b, c) > REL_TOL
        and parallelepiped_check(a, b, c, alpha, beta, gamma)
    ):
        # Prepare cell
        cell = rotate(
            Cell.from_params(a, b, c, alpha, beta, gamma),
            r1,
            r2,
            r3,
        )
        l = Lattice(cell, standardize=False)

        assert np.allclose(
            [l.k_a, l.k_b, l.k_c],
            list(l.reciprocal_parameters[:3]),
            rtol=REL_TOL,
            atol=ABS_TOL,
        )
        assert np.allclose(
            [l.k_alpha, l.k_beta, l.k_gamma],
            list(l.reciprocal_parameters[3:]),
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )
        r_vol = (2 * pi) ** 3 / abs(np.linalg.det(l.cell))
        assert np.allclose(
            l.reciprocal_cell_volume,
            r_vol,
            rtol=REL_TOL,
            atol=ABS_TOL,
        )


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=0.0, max_value=1 / ABS_TOL),
)
def test_Lattice_eps(r1, r2, r3, a, b, c, alpha, beta, gamma, eps_rel):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(a, b, c) / max(a, b, c) > REL_TOL
        and parallelepiped_check(a, b, c, alpha, beta, gamma)
    ):
        # Prepare cell
        cell = rotate(
            Cell.from_params(a, b, c, alpha, beta, gamma),
            r1,
            r2,
            r3,
        )
        l = Lattice(cell, standardize=False)
        l.eps_rel = eps_rel
        assert l.eps_rel == eps_rel
        assert np.allclose(
            l.eps,
            eps_rel * abs(l.unit_cell_volume) ** (1.0 / 3.0),
            rtol=REL_TOL,
            atol=ABS_TOL,
        )


@pytest.mark.parametrize(
    "variation", BRAVAIS_LATTICE_VARIATIONS, ids=BRAVAIS_LATTICE_VARIATIONS
)
def test_Lattice_classifications(variation: str):
    lattice = lattice_example(variation)
    type_name = variation.translate(str.maketrans("", "", "12345ab"))
    assert lattice.type() == type_name
    assert lattice.variation == variation
    assert lattice.name == BRAVAIS_LATTICE_NAMES[type_name]
    assert lattice.pearson_symbol == PEARSON_SYMBOLS[type_name]
    assert lattice.crystal_family == lattice.pearson_symbol[0]
    assert lattice.centring_type == lattice.pearson_symbol[1]


# legacy tests

l = Lattice([1, 0, 0], [0, 2, 0], [0, 0, 3])


def test_extended_init():
    with pytest.raises(ValueError):
        l = Lattice(1, 2)
    l = Lattice(1, 1, 1, 90, 90, 60)
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10
    l = Lattice([1, 0, 0], [1 / 2, sqrt(3) / 2, 0], [0, 0, 1])
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10
    l = Lattice([[1, 0, 0], [1 / 2, sqrt(3) / 2, 0], [0, 0, 1]])
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10


def test_cell():
    assert (l.cell == np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])).all()
    assert (l.a1 == np.array([1, 0, 0])).all()
    assert (l.a2 == np.array([0, 2, 0])).all()
    assert (l.a3 == np.array([0, 0, 3])).all()
    assert l.unit_cell_volume == 6
    assert l.a == 1
    assert l.b == 2
    assert l.c == 3
    assert l.alpha == l.beta == l.gamma == 90


def test_reciprocal():
    assert (
        l.reciprocal_cell == np.array([[2 * pi, 0, 0], [0, pi, 0], [0, 0, 2 / 3 * pi]])
    ).all()
    assert np.allclose(l.b1, [2 * pi, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(l.b2, [0, pi, 0], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(l.b3, [0, 0, 2 / 3 * pi], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(
        l.reciprocal_cell_volume, 4 / 3 * pi**3, rtol=REL_TOL, atol=ABS_TOL
    )
    assert np.allclose(
        [l.k_a, l.k_b, l.k_c], [2 * pi, pi, 2 / 3 * pi], rtol=REL_TOL, atol=ABS_TOL
    )

    assert np.allclose(
        [l.k_alpha, l.k_beta, l.k_gamma], [90.0, 90.0, 90.0], rtol=REL_TOL, atol=ABS_TOL
    )


def test_variation():
    assert l.variation == "ORC"
