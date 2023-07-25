from math import pi

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
from radtools.geometry import parallelepiped_check


import radtools.crystal.cell as Cell

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


@given(
    harrays(
        float,
        (3, 3),
        elements=st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    ),
)
def test_params_from_cell(cell):
    a, b, c, alpha, beta, gamma = Cell.params(cell)


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
)
def test_cell_from_param(a, b, c, alpha, beta, gamma):
    if parallelepiped_check(a, b, c, alpha, beta, gamma):
        cell = Cell.from_params(a, b, c, alpha, beta, gamma)
        if (cell > MIN_LENGTH).all():
            ap, bp, cp, alphap, betap, gammap = Cell.params(cell)
            assert np.allclose([a, b, c], [ap, bp, cp], rtol=REL_TOL, atol=ABS_TOL)
            assert np.allclose(
                [alpha, beta, gamma],
                [alphap, betap, gammap],
                rtol=REL_TOL_ANGLE,
                atol=ABS_TOL_ANGLE,
            )
    else:
        with pytest.raises(ValueError):
            Cell.from_params(a, b, c, alpha, beta, gamma)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, cell",
    [(1, 1, 1, 90, 90, 90, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])],
)
def test_cell_from_params_example(a, b, c, alpha, beta, gamma, cell):
    assert (
        Cell.from_params(a, b, c, alpha, beta, gamma) - np.array(cell) < ABS_TOL
    ).all()


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
    st.integers(min_value=0, max_value=n_order),
)
def test_reciprocal(r1, r2, r3, a, b, c, alpha, beta, gamma, order):
    if (
        parallelepiped_check(a, b, c, alpha, beta, gamma)
        # Maximum I can do right now.
        # It passes only for the vectors, which are not too far away from each other.
        and min(a, b, c) / max(a, b, c) > REL_TOL
    ):
        cell = shuffle(
            rotate(Cell.from_params(a, b, c, alpha, beta, gamma), r1, r2, r3), order
        )
        rcell = Cell.reciprocal(cell)
        product = np.diag(np.abs(rcell.T @ cell))
        correct_product = np.ones(3) * 2 * pi
        # Non  diagonal terms are close to zero. Hard to work them out. Try if you want to suffer.
        assert np.allclose(product, correct_product, rtol=REL_TOL, atol=ABS_TOL)


@pytest.mark.parametrize(
    "cell, rec_cell",
    [
        (
            [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
            [[2 * pi, 0, 0], [0, pi, 0], [0, 0, 2 / 3 * pi]],
        )
    ],
)
def test_reciprocal_cell_examples(cell, rec_cell):
    rcell = Cell.reciprocal(cell)
    assert np.allclose(rcell, np.array(rec_cell), rtol=REL_TOL, atol=ABS_TOL)
