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

from math import cos, pi, sin

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays
from scipy.spatial.transform import Rotation

from radtools.constants import TORADIANS
from radtools.crystal.constants import (
    ABS_TOL,
    ABS_TOL_ANGLE,
    MAX_LENGTH,
    MIN_ANGLE,
    MIN_LENGTH,
)
from radtools.geometry import (
    absolute_to_relative,
    angle,
    parallelepiped_check,
    span_orthonormal_set,
    volume,
)
from radtools.numerical import compare_numerically

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
    return R.T @ cell


@pytest.mark.parametrize(
    "cell, absolute, relative",
    [
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0], [0, 0, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 1], [0, 0, 1]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 1, 0], [0, 1, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 0, 0], [1, 0, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.5, 0], [0.5, 0.5, 0]),
        ([[1, 1, 0], [0, 1, 0], [0, 0, 1]], [0.5, 1, 0], [0.5, 0.5, 0]),
        ([[2, 1, 0], [1, 1, 0], [0, 0, 1]], [0.9, 0.7, 0.4], [0.2, 0.5, 0.4]),
    ],
)
def test_absolute_to_relative(cell, absolute, relative):
    new_relative = absolute_to_relative(cell, absolute)
    assert (new_relative == relative).all()


@given(
    harrays(float, 3, elements=st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH)),
    harrays(
        float,
        3,
        elements=st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    ),
)
def test_angle(v1, v2):
    if (
        abs(np.linalg.norm(v1)) > np.finfo(float).eps
        and abs(np.linalg.norm(v2)) > np.finfo(float).eps
    ):
        result_degrees = angle(v1, v2)
        result_radians = angle(v1, v2, radians=True)
        assert 0 <= result_degrees <= 180
        assert 0 <= result_radians <= pi


@given(st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE))
def test_angle_values(alpha):
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([cos(alpha * TORADIANS), sin(alpha * TORADIANS), 0.0])
    assert abs(angle(v1, v2) - alpha) < ABS_TOL_ANGLE


def test_angle_raises():
    with pytest.raises(ValueError):
        angle([0, 0, 0], [0, 0, 0])
    with pytest.raises(ValueError):
        angle([0, 0, 0], [1.0, 0, 0])


@pytest.mark.parametrize(
    "args, result, eps", [((4, 4.472, 4.583, 79.03, 64.13, 64.15), 66.3840797, 1e-8)]
)
def test_volume_example(args, result, eps):
    assert volume(*args) - result < eps


# no need to test the vectors - they take the same route as the cell
# it the vectors will move from rows to columns it is not a problem as well
# since the volume of the cell is its determinant
@example(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
@example(np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]))
@given(
    harrays(
        float,
        (3, 3),
        elements=st.floats(
            min_value=-MAX_LENGTH,
            max_value=MAX_LENGTH,
        ),
    ),
)
def test_volume_with_cell(cell):
    # Its an "or" condition
    if ((np.abs(cell) > MIN_LENGTH) + (cell == np.zeros(cell.shape))).all():
        if np.linalg.det(cell) > 0:
            assert volume(cell) >= 0
        else:
            assert volume(cell) <= 0


@given(
    st.floats(
        min_value=MIN_LENGTH,
        max_value=MAX_LENGTH,
    ),
    st.floats(
        min_value=MIN_LENGTH,
        max_value=MAX_LENGTH,
    ),
    st.floats(
        min_value=MIN_LENGTH,
        max_value=MAX_LENGTH,
    ),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_volume_parameters(a, b, c, alpha, beta, gamma):
    if parallelepiped_check(a, b, c, alpha, beta, gamma):
        assert volume(a, b, c, alpha, beta, gamma) > 0


@pytest.mark.parametrize(
    "e3",
    [
        ([0, 0, 1]),
        ([0, 1, 0]),
        ([1, 0, 0]),
        ([1, 1, 0]),
        ([1, 0, 1]),
        ([0, 1, 1]),
    ],
)
def test_span_orthonormal_set(e3):
    e1p, e2p, e3p = span_orthonormal_set(e3)

    # Check if the third vector is the same as the normalized input vector
    assert np.allclose(e3p, np.array(e3) / np.linalg.norm(e3))

    # Check if all vectors are orthogonal to each other
    assert np.allclose(np.dot(e1p, e2p), 0)
    assert np.allclose(np.dot(e1p, e3p), 0)
    assert np.allclose(np.dot(e2p, e3p), 0)

    # Check if all vectors are normalized
    assert np.allclose(np.linalg.norm(e1p), 1)
    assert np.allclose(np.linalg.norm(e2p), 1)
    assert np.allclose(np.linalg.norm(e3p), 1)

    # Check if the system is right-handed
    assert np.allclose(np.dot(e1p, np.cross(e2p, e3p)), 1)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_parallelepiped_check(a, b, c, alpha, beta, gamma):
    assert parallelepiped_check(a, b, c, alpha, beta, gamma) == (
        compare_numerically(a, ">", 0.0, ABS_TOL)
        and compare_numerically(b, ">", 0.0, ABS_TOL)
        and compare_numerically(c, ">", 0.0, ABS_TOL)
        and compare_numerically(alpha, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(beta, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(alpha, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(beta, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, "<", alpha + beta, ABS_TOL_ANGLE)
        and compare_numerically(alpha + beta, "<", 360.0 - gamma, ABS_TOL_ANGLE)
        and compare_numerically(beta, "<", alpha + gamma, ABS_TOL_ANGLE)
        and compare_numerically(alpha + gamma, "<", 360.0 - beta, ABS_TOL_ANGLE)
        and compare_numerically(alpha, "<", beta + gamma, ABS_TOL_ANGLE)
        and compare_numerically(beta + gamma, "<", 360.0 - alpha, ABS_TOL_ANGLE)
    )
