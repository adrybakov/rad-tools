from radtools.crystal.bravais_lattice.utils import (
    lattice_example,
)

from radtools.routines import volume

from radtools.crystal.lattice import Lattice

from radtools.crystal.constants import (
    BRAVAIS_LATTICE_VARIATIONS,
    MAX_LENGTH,
    MIN_LENGTH,
    MIN_ANGLE,
    REL_TOL,
    ABS_TOL,
)

import pytest
from hypothesis import given, example, strategies as st, settings
from hypothesis.extra.numpy import arrays as harrays
import numpy as np

from scipy.spatial.transform import Rotation


from math import pi

from radtools.routines import parallelepiped_check, cell_from_param


@pytest.mark.parametrize(
    "lattice_variation", BRAVAIS_LATTICE_VARIATIONS, ids=BRAVAIS_LATTICE_VARIATIONS
)
def test_lattice_example(lattice_variation):
    lattice_type = lattice_variation.translate(str.maketrans("", "", "12345ab"))
    assert lattice_example(lattice_variation).type() == lattice_type
    assert lattice_example(lattice_variation).variation == lattice_variation


@given(st.text())
def test_lattice_example_error(wrong_name):
    if wrong_name.lower() not in list(
        map(lambda x: x.lower(), BRAVAIS_LATTICE_VARIATIONS)
    ):
        with pytest.raises(ValueError):
            lattice_example(wrong_name)


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
