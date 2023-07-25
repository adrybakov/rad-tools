from math import cos, pi, sin

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays
from scipy.spatial.transform import Rotation

from radtools.crystal.constants import (
    ABS_TOL,
    ABS_TOL_ANGLE,
    MAX_LENGTH,
    MIN_ANGLE,
    MIN_LENGTH,
    REL_TOL,
    REL_TOL_ANGLE,
)
from radtools.utils import (
    absolute_to_relative,
    angle,
    atom_mark_to_latex,
    compare_numerically,
    parallelepiped_check,
    rot_angle,
    span_orthonormal_set,
    toradians,
    volume,
)


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
    "mark, new_mark",
    [
        ("Cr1", "$Cr_{1}$"),
        ("Cr11", "$Cr_{11}$"),
    ],
)
def test_atom_mark_to_latex(mark, new_mark):
    assert atom_mark_to_latex(mark) == new_mark


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
    v2 = np.array([cos(alpha * toradians), sin(alpha * toradians), 0.0])
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
# it ht e vectors will move from rows to columns it is not a problem as well
# since the volume of the cell is its determinant
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
    if ((np.abs(cell) > MIN_LENGTH) + (cell == 0)).all():
        assert volume(cell) >= 0


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


class TestRotAngle:
    @pytest.mark.parametrize(
        "x, y, angle",
        [
            (1, 0, 0),
            (1, 1, 45),
            (0, 1, 90),
            (-1, 1, 135),
            (-1, 0, 180),
            (-1, -1, 225),
            (0, -1, 270),
            (1, -1, 315),
        ],
    )
    def test_dummy(self, x, y, angle):
        assert round(rot_angle(x, y, dummy=True), 4) == round(angle, 4)

    @pytest.mark.parametrize(
        "x, y, angle",
        [
            (1, 0, 0),
            (1, 1, 45),
            (0, 1, 90),
            (-1, 1, -45),
            (-1, 0, 0),
            (-1, -1, 45),
            (0, -1, 90),
            (1, -1, -45),
        ],
    )
    def test_not_dummy(self, x, y, angle):
        assert round(rot_angle(x, y), 4) == round(angle, 4)

    def test_ill_case(self):
        with pytest.raises(ValueError):
            rot_angle(0, 0)


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
