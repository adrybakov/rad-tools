from radtools.crystal.bravais_lattice.utils import (
    bravais_lattice_from_param,
    bravais_lattice_from_cell,
    lattice_example,
)

from radtools.routines import volume

from radtools.crystal.constants import BRAVAIS_LATTICE_VARIATIONS

import pytest
from hypothesis import given, example, strategies as st
from hypothesis.extra.numpy import arrays as harrays
import numpy as np


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


# TODO: test taking too long
@given(harrays(float, (3, 3)))
def test_bravais_lattice_from_cell(cell):
    print(cell)
    if np.linalg.det(cell) != 0:
        new_cell = bravais_lattice_from_cell(cell)
        assert np.allclose(volume(cell), volume(new_cell))
