# NOTE Relative size of lattice parameters have to be consistent with the REL_TOL, in order for the tests to be correct.

from hypothesis import given, example, strategies as st

from radtools.crystal.bravais_lattice.cells import (
    CUB_cell,
    FCC_cell,
    BCC_cell,
    TET_cell,
    BCT_cell,
    ORC_cell,
    ORCF_cell,
    ORCI_cell,
    ORCC_cell,
    HEX_cell,
    RHL_cell,
    MCL_cell,
    MCLC_cell,
    TRI_cell,
    CUB_fix_cell,
    FCC_fix_cell,
    BCC_fix_cell,
    TET_fix_cell,
    BCT_fix_cell,
    ORC_fix_cell,
    ORCF_fix_cell,
    ORCI_fix_cell,
    ORCC_fix_cell,
    HEX_fix_cell,
    RHL_fix_cell,
    MCL_fix_cell,
    MCLC_fix_cell,
    TRI_fix_cell,
)

from scipy.spatial.transform import Rotation
from radtools.routines import print_2d_array, param_from_cell, todegrees

import numpy as np
from math import sqrt, cos, sin, acos
import pytest

from radtools.crystal.constants import (
    MAX_LENGTH,
    MIN_LENGTH,
    MIN_ANGLE,
    REL_TOL,
    ABS_TOL,
    ABS_TOL_ANGLE,
    REL_TOL_ANGLE,
)

from radtools.routines import toradians


@given(st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH))
def test_CUB_cell(a):
    cell = CUB_cell(a)
    assert np.allclose(cell, np.eye(3) * a, rtol=REL_TOL, atol=ABS_TOL)


@given(st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH))
def test_FCC_cell(a):
    cell = FCC_cell(a)
    assert np.allclose(
        cell, (np.ones((3, 3)) - np.eye(3)) * a / 2, rtol=REL_TOL, atol=ABS_TOL
    )


@given(st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH))
def test_BCC_cell(a):
    cell = BCC_cell(a)
    assert np.allclose(
        cell, (np.ones((3, 3)) - 2 * np.eye(3)) * a / 2, rtol=REL_TOL, atol=ABS_TOL
    )


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_TET_cell(a, c):
    cell = TET_cell(a, c)
    assert np.allclose(cell, np.diag([a, a, c]), rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_BCT_cell(a, c):
    cell = BCT_cell(a, c)
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, :2] *= a
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_ORC_cell(a, b, c):
    cell = ORC_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    assert np.allclose(cell, np.diag([a, b, c]), rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_ORCF_cell(a, b, c):
    cell = ORCF_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_ORCI_cell(a, b, c):
    cell = ORCI_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_ORCC_cell(a, b, c):
    cell = ORCC_cell(a, b, c)
    a, b = sorted([a, b])
    correct_cell = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]]) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
)
def test_HEX_cell(a, c):
    cell = HEX_cell(a, c)
    correct_cell = np.array(
        [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
    )
    assert np.allclose(cell, correct_cell, rtol=REL_TOL, atol=ABS_TOL)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=120.0 - MIN_ANGLE,
        exclude_max=True,
    ),
)
def test_RHL_cell(a, alpha):
    cell = RHL_cell(a, alpha)
    alpha *= toradians
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
        RHL_cell(a, 120)


@given(
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=180.0 - MIN_ANGLE,
        exclude_max=True,
    ),
)
def test_MCL_cell(a, b, c, alpha):
    cell = MCL_cell(a, b, c, alpha)
    if alpha > 90:
        alpha = 180 - alpha
    alpha *= toradians
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
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=180.0 - MIN_ANGLE,
        exclude_max=True,
    ),
)
def test_MCLC_cell(a, b, c, alpha):
    cell = MCLC_cell(a, b, c, alpha)
    if alpha > 90:
        alpha = 180.0 - alpha
    alpha *= toradians
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
# def test_TRI_cell(a, b, c, alpha, beta, gamma, reciprocal):
#     cell = TRI_cell(a, b, c, alpha, beta, gamma, reciprocal)
#     pass

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
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_CUB_fix_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(CUB_cell(conv_a), r1, r2, r3), order)
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = CUB_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose([a, b, c], [a, a, a], rtol=REL_TOL, atol=ABS_TOL)
        assert np.allclose(
            [alpha, beta, gamma],
            [90.0, 90.0, 90.0],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_FCC_fix_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(FCC_cell(conv_a), r1, r2, r3), order)
        prim_a = conv_a * sqrt(2) / 2
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = FCC_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_a, prim_a], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            [alpha, beta, gamma],
            [60.0, 60.0, 60.0],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_BCC_fix_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(BCC_cell(conv_a), r1, r2, r3), order)
        angle = acos(-1 / 3) * todegrees
        prim_a = conv_a * sqrt(3) / 2
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = BCC_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_a, prim_a], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            [alpha, beta, gamma],
            [angle, angle, angle],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_TET_fix_cell(r1, r2, r3, conv_a, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(TET_cell(conv_a, conv_c), r1, r2, r3), order)
        angle = 90
        prim_a = conv_a
        prim_c = conv_c
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = TET_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c],
            [prim_a, prim_a, prim_c],
            rtol=REL_TOL,
            atol=ABS_TOL,
        )
        assert np.allclose(
            [alpha, beta, gamma],
            [angle, angle, angle],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_BCT_fix_cell(r1, r2, r3, conv_a, conv_c, order):
    if (
        not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
        and min(conv_a, conv_c) / max(conv_a, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(BCT_cell(conv_a, conv_c), r1, r2, r3), order)
        prim = sqrt(2 * conv_a**2 + conv_c**2) / 2
        angle12 = acos((conv_c**2 - 2 * conv_a**2) / 4 / prim**2)
        angle = acos(-(conv_c**2) / 4 / prim**2)
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = BCT_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        alpha *= toradians
        beta *= toradians
        gamma *= toradians
        assert np.allclose([a, b, c], [prim, prim, prim], rtol=REL_TOL, atol=ABS_TOL)
        assert np.allclose(
            [alpha, beta, gamma],
            [angle, angle, angle12],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORC_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if (
        not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
        and min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(ORC_cell(conv_a, conv_b, conv_c), r1, r2, r3), order)
        prim_a, prim_b, prim_c = sorted([conv_a, conv_b, conv_c])
        angle = 90
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORC_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            [alpha, beta, gamma],
            [angle, angle, angle],
            rtol=REL_TOL_ANGLE,
            atol=ABS_TOL_ANGLE,
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCF_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(ORCF_cell(conv_a, conv_b, conv_c), r1, r2, r3), order)
        conv_a, conv_b, conv_c = sorted([conv_a, conv_b, conv_c])
        prim_a = sqrt(conv_b**2 + conv_c**2) / 2.0
        prim_b = sqrt(conv_a**2 + conv_c**2) / 2.0
        prim_c = sqrt(conv_a**2 + conv_b**2) / 2.0
        prim_alpha = acos(conv_a**2 / 4.0 / prim_b / prim_c) * todegrees
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = acos(conv_b**2 / 4.0 / prim_a / prim_c) * todegrees
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = acos(conv_c**2 / 4.0 / prim_a / prim_b) * todegrees
        prim_gamma_twin = 180.0 - prim_gamma
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORCF_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(alpha, prim_alpha_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(
            beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(beta, prim_beta_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(
            gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(gamma, prim_gamma_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCI_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(ORCI_cell(conv_a, conv_b, conv_c), r1, r2, r3), order)
        conv_a, conv_b, conv_c = sorted([conv_a, conv_b, conv_c])
        prim = sqrt(conv_a**2 + conv_b**2 + conv_c**2) / 2
        prim_alpha = (
            acos((conv_a**2 - conv_b**2 - conv_c**2) / 4.0 / prim**2)
            * todegrees
        )
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = (
            acos((-(conv_a**2) + conv_b**2 - conv_c**2) / 4.0 / prim**2)
            * todegrees
        )
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = (
            acos((-(conv_a**2) - conv_b**2 + conv_c**2) / 4.0 / prim**2)
            * todegrees
        )
        prim_gamma_twin = 180.0 - prim_gamma
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORCI_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose([a, b, c], [prim, prim, prim], rtol=REL_TOL, atol=ABS_TOL)
        assert np.allclose(
            alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(alpha, prim_alpha_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(
            beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(beta, prim_beta_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(
            gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(gamma, prim_gamma_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCC_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        cell = shuffle(rotate(ORCC_cell(conv_a, conv_b, conv_c), r1, r2, r3), order)

        conv_a, conv_b = sorted([conv_a, conv_b])
        prim_a = sqrt(conv_a**2 + conv_b**2) / 2
        prim_b = prim_a
        prim_c = conv_c
        prim_alpha = 90.0
        prim_beta = prim_alpha
        prim_gamma = (
            acos((conv_a**2 - conv_b**2) / 4.0 / prim_a / prim_b) * todegrees
        )

        old_det = np.linalg.det(cell)

        cell = ORCC_fix_cell(cell, REL_TOL)
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            alpha, prim_alpha, rtol=REL_TOL, atol=ABS_TOL
        ) or np.allclose(alpha, 180.0 - prim_alpha, rtol=REL_TOL, atol=ABS_TOL)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL, atol=ABS_TOL) or np.allclose(
            beta, 180.0 - prim_beta, rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            gamma, prim_gamma, rtol=REL_TOL, atol=ABS_TOL
        ) or np.allclose(gamma, 180.0 - prim_gamma, rtol=REL_TOL, atol=ABS_TOL)

        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_HEX_fix_cell(r1, r2, r3, conv_a, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_c) / max(conv_a, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(HEX_cell(conv_a, conv_c), r1, r2, r3), order)
        prim_a = conv_a
        prim_b = conv_a
        prim_c = conv_c
        prim_alpha = 90.0
        prim_beta = 90.0
        prim_gamma = 120.0
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = HEX_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=120.0 - MIN_ANGLE,
        exclude_max=True,
    ),
    st.integers(min_value=0, max_value=n_order),
)
def test_RHL_fix_cell(r1, r2, r3, conv_a, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(RHL_cell(conv_a, conv_alpha), r1, r2, r3), order)
        prim_a = conv_a
        prim_b = conv_a
        prim_c = conv_a
        prim_alpha = conv_alpha
        prim_beta = conv_alpha
        prim_gamma = conv_alpha
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = RHL_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=180.0 - MIN_ANGLE,
        exclude_max=True,
    ),
    st.integers(min_value=0, max_value=n_order),
)
def test_MCL_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(MCL_cell(conv_a, conv_b, conv_c, conv_alpha), r1, r2, r3), order
        )
        prim_a = conv_a
        prim_b, prim_c = sorted([conv_b, conv_c])
        prim_beta = 90.0
        if conv_alpha > 90.0:
            prim_alpha = 180.0 - conv_alpha
        else:
            prim_alpha = conv_alpha
        prim_gamma = 90.0
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = MCL_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, exclude_min=True, max_value=MAX_LENGTH),
    st.floats(
        min_value=MIN_ANGLE,
        exclude_min=True,
        max_value=180.0 - MIN_ANGLE,
        exclude_max=True,
    ),
    st.integers(min_value=0, max_value=n_order),
)
def test_MCLC_fix_cell(r1, r2, r3, conv_a, conv_b, conv_c, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(MCLC_cell(conv_a, conv_b, conv_c, conv_alpha), r1, r2, r3), order
        )
        conv_b, conv_c = sorted([conv_b, conv_c])
        if conv_alpha > 90.0:
            conv_alpha = 180.0 - conv_alpha
        else:
            conv_alpha = conv_alpha
        prim_a = sqrt(conv_a**2 + conv_b**2) / 2
        prim_b = sqrt(conv_a**2 + conv_b**2) / 2
        prim_c = conv_c
        prim_alpha = (
            acos(conv_b * conv_c * cos(conv_alpha * toradians) / 2.0 / prim_b / prim_c)
            * todegrees
        )
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = (
            acos(conv_b * conv_c * cos(conv_alpha * toradians) / 2.0 / prim_a / prim_c)
            * todegrees
        )
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = (
            acos((conv_b**2 - conv_a**2) / 4.0 / prim_a / prim_b) * todegrees
        )
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = MCLC_fix_cell(cell, REL_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = param_from_cell(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(
            alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(alpha, prim_alpha_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(
            beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
        ) or np.allclose(beta, prim_beta_twin, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


# @given(
#     st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
#     st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
#     st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
#     st.floats(min_value=0.1, exclude_min=True, max_value=MAX_LENGTH),
#     st.floats(min_value=0.1, exclude_min=True, max_value=MAX_LENGTH),
#     st.floats(min_value=0.1, exclude_min=True, max_value=MAX_LENGTH),
#     st.floats(min_value=0.1, exclude_min=True, max_value=179.9, exclude_max=True),
#     st.floats(min_value=0.1, exclude_min=True, max_value=179.9, exclude_max=True),
#     st.floats(min_value=0.1, exclude_min=True, max_value=179.9, exclude_max=True),
#     st.integers(min_value=0, max_value=n_order),
# )
# def test_TRI_fix_cell(
#     r1, r2, r3, conv_a, conv_b, conv_c, conv_alpha, conv_beta, conv_gamma, order
# ):
#     if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL,atol=ABS_TOL):
#         cell = shuffle(
#             rotate(
#                 TRI_cell(conv_a, conv_b, conv_c, conv_alpha, conv_beta, conv_gamma),
#                 r1,
#                 r2,
#                 r3,
#             ),
#             order,
#         )

#         old_det = np.linalg.det(cell)
#         # TODO: Decide the logic of parameter`s order
#         cell = TRI_fix_cell(cell, REL_TOL)
#         a, b, c, alpha, beta, gamma = param_from_cell(cell)

#         assert np.allclose(
#             [a, b, c, alpha, beta, gamma],
#             [conv_a, conv_b, conv_c, conv_alpha, conv_beta, conv_gamma],
#         , rtol=REL_TOL,atol=ABS_TOL)

#         assert np.linalg.det(cell) * old_det > 0
