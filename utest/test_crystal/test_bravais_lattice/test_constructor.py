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

from math import cos, sin, sqrt

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from radtools.constants import TORADIANS
from radtools.crystal.bravais_lattice.constructor import (
    BCC,
    BCT,
    CUB,
    FCC,
    HEX,
    MCL,
    MCLC,
    ORC,
    ORCC,
    ORCF,
    ORCI,
    RHL,
    TET,
    TRI,
)
from radtools.crystal.constants import (
    ABS_TOL,
    MAX_LENGTH,
    MIN_ANGLE,
    MIN_LENGTH,
    REL_TOL,
)


@given(st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH))
def test_CUB(a):
    cell = CUB(a, return_cell=True)
    assert np.allclose(cell, np.eye(3) * a, rtol=REL_TOL, atol=ABS_TOL)


@given(st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH))
def test_FCC(a):
    cell = FCC(a, return_cell=True)
    assert np.allclose(
        cell, (np.ones((3, 3)) - np.eye(3)) * a / 2, rtol=REL_TOL, atol=ABS_TOL
    )


@given(st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH))
def test_BCC(a):
    cell = BCC(a, return_cell=True)
    assert np.allclose(
        cell, (np.ones((3, 3)) - 2 * np.eye(3)) * a / 2, rtol=REL_TOL, atol=ABS_TOL
    )


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_TET(a, c):
    cell = TET(a, c, return_cell=True)
    assert np.allclose(cell, np.diag([a, a, c]), rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_BCT(a, c):
    cell = BCT(a, c, return_cell=True)
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, :2] *= a
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_ORC(a, b, c):
    cell = ORC(a, b, c, return_cell=True)
    a, b, c = sorted([a, b, c])
    assert np.allclose(cell, np.diag([a, b, c]), rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_ORCF(a, b, c):
    cell = ORCF(a, b, c, return_cell=True)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_ORCI(a, b, c):
    cell = ORCI(a, b, c, return_cell=True)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_ORCC(a, b, c):
    cell = ORCC(a, b, c, return_cell=True)
    a, b = sorted([a, b])
    correct_cell = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]]) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
)
def test_HEX(a, c):
    cell = HEX(a, c, return_cell=True)
    correct_cell = np.array(
        [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
    )
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=120.0 - MIN_ANGLE),
)
def test_RHL(a, alpha):
    cell = RHL(a, alpha, return_cell=True)
    alpha *= TORADIANS
    correct_cell = np.array(
        [
            [a * cos(alpha / 2), -a * sin(alpha / 2), 0],
            [a * cos(alpha / 2), a * sin(alpha / 2), 0],
            [
                a * cos(alpha) / cos(alpha / 2),
                0,
                a * sqrt(1 - cos(alpha) ** 2 / cos(alpha / 2) ** 2),
            ],
        ]
    )
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)

    with pytest.raises(ValueError):
        RHL(a, 120)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_MCL(a, b, c, alpha):
    cell = MCL(a, b, c, alpha, return_cell=True)
    if alpha > 90:
        alpha = 180 - alpha
    alpha *= TORADIANS
    b, c = sorted([b, c])
    correct_cell = np.array(
        [
            [a, 0, 0],
            [0, b, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_MCLC(a, b, c, alpha):
    cell = MCLC(a, b, c, alpha, return_cell=True)
    if alpha > 90:
        alpha = 180.0 - alpha
    alpha *= TORADIANS
    b, c = sorted([b, c])
    correct_cell = np.array(
        [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


# # TODO Test trigonal
# full_set = []


# @pytest.mark.parametrize("a, b, c, alpha, beta, gamma, reciprocal", full_set)
# def test_TRI(a, b, c, alpha, beta, gamma, reciprocal):
#     cell = TRI(a, b, c, alpha, beta, gamma, reciprocal, return_cell = True)
#     pass
