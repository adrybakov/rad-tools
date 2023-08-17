# NOTE Relative size of lattice parameters have to be consistent with the REL_TOL, in order for the tests to be correct.

from math import acos, cos, sqrt, pi

import numpy as np
from hypothesis import given
from hypothesis import strategies as st
from scipy.spatial.transform import Rotation

from radtools.crystal.bravais_lattice.standardize import (
    BCC_standardize_cell,
    BCT_standardize_cell,
    CUB_standardize_cell,
    FCC_standardize_cell,
    HEX_standardize_cell,
    MCL_standardize_cell,
    MCLC_standardize_cell,
    ORC_standardize_cell,
    ORCC_standardize_cell,
    ORCF_standardize_cell,
    ORCI_standardize_cell,
    RHL_standardize_cell,
    TET_standardize_cell,
    TRI_standardize_cell,
)
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
    ABS_TOL_ANGLE,
    MAX_LENGTH,
    MIN_ANGLE,
    MIN_LENGTH,
    REL_TOL,
    REL_TOL_ANGLE,
)
from radtools.constants import TODEGREES, TORADIANS

import radtools.crystal.cell as Cell

from radtools.numerical import compare_numerically
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
    st.integers(min_value=0, max_value=n_order),
)
def test_CUB_standardize_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(CUB(conv_a, return_cell=True), r1, r2, r3), order)
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = CUB_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_FCC_standardize_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(FCC(conv_a, return_cell=True), r1, r2, r3), order)
        prim_a = conv_a * sqrt(2) / 2
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = FCC_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_BCC_standardize_cell(r1, r2, r3, conv_a, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(BCC(conv_a, return_cell=True), r1, r2, r3), order)
        angle = acos(-1 / 3) * TODEGREES
        prim_a = conv_a * sqrt(3) / 2
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = BCC_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_TET_standardize_cell(r1, r2, r3, conv_a, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(rotate(TET(conv_a, conv_c, return_cell=True), r1, r2, r3), order)
        angle = 90
        prim_a = conv_a
        prim_c = conv_c
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = TET_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_BCT_standardize_cell(r1, r2, r3, conv_a, conv_c, order):
    if (
        not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
        and min(conv_a, conv_c) / max(conv_a, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(BCT(conv_a, conv_c, return_cell=True), r1, r2, r3), order)
        prim = sqrt(2 * conv_a**2 + conv_c**2) / 2
        angle12 = acos((conv_c**2 - 2 * conv_a**2) / 4 / prim**2)
        angle = acos(-(conv_c**2) / 4 / prim**2)
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = BCT_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
        alpha *= TORADIANS
        beta *= TORADIANS
        gamma *= TORADIANS
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORC_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if (
        not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
        and min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(ORC(conv_a, conv_b, conv_c, return_cell=True), r1, r2, r3), order
        )
        prim_a, prim_b, prim_c = sorted([conv_a, conv_b, conv_c])
        angle = 90
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORC_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)
        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCF_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(ORCF(conv_a, conv_b, conv_c, return_cell=True), r1, r2, r3), order
        )
        conv_a, conv_b, conv_c = sorted([conv_a, conv_b, conv_c])
        prim_a = sqrt(conv_b**2 + conv_c**2) / 2.0
        prim_b = sqrt(conv_a**2 + conv_c**2) / 2.0
        prim_c = sqrt(conv_a**2 + conv_b**2) / 2.0
        prim_alpha = acos(conv_a**2 / 4.0 / prim_b / prim_c) * TODEGREES
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = acos(conv_b**2 / 4.0 / prim_a / prim_c) * TODEGREES
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = acos(conv_c**2 / 4.0 / prim_a / prim_b) * TODEGREES
        prim_gamma_twin = 180.0 - prim_gamma
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORCF_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCI_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(ORCI(conv_a, conv_b, conv_c, return_cell=True), r1, r2, r3), order
        )
        conv_a, conv_b, conv_c = sorted([conv_a, conv_b, conv_c])
        prim = sqrt(conv_a**2 + conv_b**2 + conv_c**2) / 2
        prim_alpha = (
            acos((conv_a**2 - conv_b**2 - conv_c**2) / 4.0 / prim**2)
            * TODEGREES
        )
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = (
            acos((-(conv_a**2) + conv_b**2 - conv_c**2) / 4.0 / prim**2)
            * TODEGREES
        )
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = (
            acos((-(conv_a**2) - conv_b**2 + conv_c**2) / 4.0 / prim**2)
            * TODEGREES
        )
        prim_gamma_twin = 180.0 - prim_gamma
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = ORCI_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_ORCC_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        cell = shuffle(
            rotate(ORCC(conv_a, conv_b, conv_c, return_cell=True), r1, r2, r3), order
        )

        conv_a, conv_b = sorted([conv_a, conv_b])
        prim_a = sqrt(conv_a**2 + conv_b**2) / 2
        prim_b = prim_a
        prim_c = conv_c
        prim_alpha = 90.0
        prim_beta = prim_alpha
        prim_gamma = (
            acos((conv_a**2 - conv_b**2) / 4.0 / prim_a / prim_b) * TODEGREES
        )

        old_det = np.linalg.det(cell)

        cell = ORCC_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.integers(min_value=0, max_value=n_order),
)
def test_HEX_standardize_cell(r1, r2, r3, conv_a, conv_c, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_c) / max(conv_a, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(rotate(HEX(conv_a, conv_c, return_cell=True), r1, r2, r3), order)
        prim_a = conv_a
        prim_b = conv_a
        prim_c = conv_c
        prim_alpha = 90.0
        prim_beta = 90.0
        prim_gamma = 120.0
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = HEX_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=120.0 - MIN_ANGLE),
    st.integers(min_value=0, max_value=n_order),
)
def test_RHL_standardize_cell(r1, r2, r3, conv_a, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL):
        # Prepare cell
        cell = shuffle(
            rotate(RHL(conv_a, conv_alpha, return_cell=True), r1, r2, r3), order
        )
        prim_a = conv_a
        prim_b = conv_a
        prim_c = conv_a
        prim_alpha = conv_alpha
        prim_beta = conv_alpha
        prim_gamma = conv_alpha
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = RHL_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.integers(min_value=0, max_value=n_order),
)
def test_MCL_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(
                MCL(conv_a, conv_b, conv_c, conv_alpha, return_cell=True), r1, r2, r3
            ),
            order,
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
        cell = MCL_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
        assert np.allclose(
            [a, b, c], [prim_a, prim_b, prim_c], rtol=REL_TOL, atol=ABS_TOL
        )
        assert np.allclose(alpha, prim_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(beta, prim_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)
        assert np.allclose(gamma, prim_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE)

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_LENGTH, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.integers(min_value=0, max_value=n_order),
)
def test_MCLC_standardize_cell(r1, r2, r3, conv_a, conv_b, conv_c, conv_alpha, order):
    if not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL) and (
        min(conv_a, conv_b, conv_c) / max(conv_a, conv_b, conv_c) > REL_TOL
    ):
        # Prepare cell
        cell = shuffle(
            rotate(
                MCLC(conv_a, conv_b, conv_c, conv_alpha, return_cell=True), r1, r2, r3
            ),
            order,
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
            acos(conv_b * conv_c * cos(conv_alpha * TORADIANS) / 2.0 / prim_b / prim_c)
            * TODEGREES
        )
        prim_alpha_twin = 180.0 - prim_alpha
        prim_beta = (
            acos(conv_b * conv_c * cos(conv_alpha * TORADIANS) / 2.0 / prim_a / prim_c)
            * TODEGREES
        )
        prim_beta_twin = 180.0 - prim_beta
        prim_gamma = (
            acos((conv_b**2 - conv_a**2) / 4.0 / prim_a / prim_b) * TODEGREES
        )
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = MCLC_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        # Check results
        a, b, c, alpha, beta, gamma = Cell.params(cell)
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


@given(
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0, max_value=2 * pi),
    st.floats(min_value=0.1, max_value=MAX_LENGTH),
    st.floats(min_value=0.1, max_value=MAX_LENGTH),
    st.floats(min_value=0.1, max_value=MAX_LENGTH),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.floats(min_value=MIN_ANGLE, max_value=180.0 - MIN_ANGLE),
    st.integers(min_value=0, max_value=n_order),
)
def test_TRI_standardize_cell(r1, r2, r3, a, b, c, alpha, beta, gamma, order):
    if (
        not np.allclose([r1, r2, r3], [0, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
        and (min(a, b, c) / max(a, b, c) > REL_TOL)
        and parallelepiped_check(a, b, c, alpha, beta, gamma)
    ):
        # Prepare cell
        cell = shuffle(
            rotate(TRI(a, b, c, alpha, beta, gamma, return_cell=True), r1, r2, r3),
            order,
        )
        prev_cell = cell
        prev_rcell = Cell.reciprocal(cell)

        k_a, k_b, k_c, k_alpha, k_beta, k_gamma = Cell.params(Cell.reciprocal(cell))
        old_det = np.linalg.det(cell)

        # Fix cell
        cell = TRI_standardize_cell(cell, rtol=REL_TOL, atol=ABS_TOL)

        s_a, s_b, s_c, s_alpha, s_beta, s_gamma = Cell.params(cell)
        s_k_a, s_k_b, s_k_c, s_k_alpha, s_k_beta, s_k_gamma = Cell.params(
            Cell.reciprocal(cell)
        )

        # Check that the cell is standardized
        if np.allclose(s_k_gamma, 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE):
            assert (
                compare_numerically(
                    s_k_alpha, "<=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_beta, "<=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                or compare_numerically(
                    s_k_alpha, ">=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_beta, ">=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
            )
        else:
            assert (
                compare_numerically(
                    s_k_alpha, "<=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_beta, "<=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_gamma, "<=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                or compare_numerically(
                    s_k_alpha, ">=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_beta, ">=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
                and compare_numerically(
                    s_k_gamma, ">=", 90, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
                )
            )

        # Check that parameters and angles are the same
        assert (
            compare_numerically(a, "==", s_a, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(a, "==", s_b, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(a, "==", s_c, rtol=REL_TOL, atol=ABS_TOL)
        )
        assert (
            compare_numerically(b, "==", s_a, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(b, "==", s_b, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(b, "==", s_c, rtol=REL_TOL, atol=ABS_TOL)
        )
        assert (
            compare_numerically(c, "==", s_a, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(c, "==", s_b, rtol=REL_TOL, atol=ABS_TOL)
            or compare_numerically(c, "==", s_c, rtol=REL_TOL, atol=ABS_TOL)
        )
        assert (
            compare_numerically(
                alpha, "==", s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                alpha, "==", s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                alpha, "==", s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                alpha, "==", 180.0 - s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                alpha, "==", 180.0 - s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                alpha, "==", 180.0 - s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
        )
        assert (
            compare_numerically(
                beta, "==", s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                beta, "==", s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                beta, "==", s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                beta, "==", 180.0 - s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                beta, "==", 180.0 - s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                beta, "==", 180.0 - s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
        )
        assert (
            compare_numerically(
                gamma, "==", s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                gamma, "==", s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                gamma, "==", s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                gamma, "==", 180.0 - s_alpha, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                gamma, "==", 180.0 - s_beta, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
            or compare_numerically(
                gamma, "==", 180.0 - s_gamma, rtol=REL_TOL_ANGLE, atol=ABS_TOL_ANGLE
            )
        )

        # Check that chirality is the same
        assert np.linalg.det(cell) * old_det > 0
