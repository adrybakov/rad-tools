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
import numpy as np
from math import sqrt, cos, sin
import pytest

from radtools.routines import toradians

a_set = [1, 1.0, 45, -1]


@pytest.mark.parametrize("a", a_set)
def test_CUB_cell(a):
    cell = CUB_cell(a)
    assert np.allclose(cell, np.eye(3) * a)


@pytest.mark.parametrize("a", a_set)
def test_FCC_cell(a):
    cell = FCC_cell(a)
    assert np.allclose(cell, (np.ones((3, 3)) - np.eye(3)) * a / 2)


@pytest.mark.parametrize("a", a_set)
def test_BCC_cell(a):
    cell = BCC_cell(a)
    assert np.allclose(cell, (np.ones((3, 3)) - 2 * np.eye(3)) * a / 2)


ac_set = [[1, 1], [4, 3], [3.2, 9], [-1, -1]]


@pytest.mark.parametrize("a, c", ac_set)
def test_TET_cell(a, c):
    cell = TET_cell(a, c)
    assert np.allclose(cell, np.diag([a, a, c]))


@pytest.mark.parametrize("a, c", ac_set)
def test_BCT_cell(a, c):
    cell = BCT_cell(a, c)
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, :2] *= a
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell)


abc_set = [[1, 1, 1], [4, 3, 2], [3.2, 9, 7], [-1, -1, 4]]


@pytest.mark.parametrize("a, b, c", abc_set)
def test_ORC_cell(a, b, c):
    cell = ORC_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    assert np.allclose(cell, np.diag([a, b, c]))


@pytest.mark.parametrize("a, b, c", abc_set)
def test_ORCF_cell(a, b, c):
    cell = ORCF_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell)


@pytest.mark.parametrize("a, b, c", abc_set)
def test_ORCI_cell(a, b, c):
    cell = ORCI_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    correct_cell = (np.ones((3, 3)) - 2 * np.eye(3)) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell)


@pytest.mark.parametrize("a, b, c", abc_set)
def test_ORCC_cell(a, b, c):
    cell = ORCC_cell(a, b, c)
    a, b, c = sorted([a, b, c])
    correct_cell = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]]) / 2
    correct_cell[:, 0] *= a
    correct_cell[:, 1] *= b
    correct_cell[:, 2] *= c
    assert np.allclose(cell, correct_cell)


ac_set = [[1, 1], [4, 3], [3.2, 9], [-1, -1]]


@pytest.mark.parametrize("a, c", ac_set)
def test_HEX_cell(a, c):
    cell = HEX_cell(a, c)
    correct_cell = np.array(
        [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
    )
    assert np.allclose(cell, correct_cell)


aalpha_set = [[1, 1], [4, 3], [3.2, 9], [-1, -1]]


@pytest.mark.parametrize("a, alpha", aalpha_set)
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
    assert np.allclose(cell, correct_cell)

    with pytest.raises(ValueError):
        RHL_cell(a, 120)


abcalpha_set = [[1, 2, 4, 90], [9, 3, 8, 7], [2.0, 4.0, 1, 78], [1, 5, 2, 145]]


@pytest.mark.parametrize("a, b, c, alpha", abcalpha_set)
def test_MCL_cell(a, b, c, alpha):
    cell = MCL_cell(a, b, c, alpha)
    if alpha > 90:
        alpha = alpha - 90
    alpha *= toradians
    b, c = sorted([b, c])
    correct_cell = np.array(
        [
            [a, 0, 0],
            [0, b, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )
    assert np.allclose(cell, correct_cell)


@pytest.mark.parametrize("a, b, c, alpha", abcalpha_set)
def test_MCLC_cell(a, b, c, alpha):
    cell = MCLC_cell(a, b, c, alpha)
    if alpha > 90:
        alpha = alpha - 90
    alpha *= toradians
    b, c = sorted([b, c])
    correct_cell = np.array(
        [
            [a / 2, b / 2, 0],
            [-a / 2, b / 2, 0],
            [0, c * cos(alpha), c * sin(alpha)],
        ]
    )
    assert np.allclose(cell, correct_cell)


# TODO
full_set = []


@pytest.mark.parametrize("a, b, c, alpha, beta, gamma, reciprocal", full_set)
def test_TRI_cell(a, b, c, alpha, beta, gamma, reciprocal):
    cell = TRI_cell(a, b, c, alpha, beta, gamma, reciprocal)
    pass


def test_CUB_fix_cell():
    for _ in range(0, 10):
        a = np.random.rand() / np.random.rand()
        assert np.allclose(np.eye(3) * a, CUB_fix_cell(np.eye(3) * a, 4))


def test_FCC_fix_cell():
    for _ in range(0, 10):
        a = np.random.rand() / np.random.rand()
        assert np.allclose(np.eye(3) * a, FCC_fix_cell(np.eye(3) * a, 4))


def test_BCC_fix_cell():
    for _ in range(0, 10):
        a = np.random.rand() / np.random.rand()
        assert np.allclose(np.eye(3) * a, BCC_fix_cell(np.eye(3) * a, 4))


cells_set = [
    ([[1, 1, 0], [0, 0, 1], [0, 1, 1]], [[0, 1, 1], [1, 1, 0], [0, 0, 1]]),
    ([[1, 1, 1], [0, 0, 1], [0, 1, 0]], [[0, 0, 1], [0, 1, 0], [1, 1, 1]]),
    ([[0, 0, 5], [3, 0, 0], [0, 3, 0]], [[3, 0, 0], [0, 3, 0], [0, 0, 5]]),
    ([[5, 0, 0], [0, 3, 0], [0, 0, 3]], [[0, 3, 0], [0, 0, 3], [5, 0, 0]]),
    ([[3, 0, 0], [0, 1, 0], [0, 0, 3]], [[0, 0, 3], [3, 0, 0], [0, 1, 0]]),
]


@pytest.mark.parametrize("cell, fixed_cell", cells_set)
def test_TET_fix_cell(cell, fixed_cell):
    assert np.allclose(fixed_cell, TET_fix_cell(cell, 1e-5))
