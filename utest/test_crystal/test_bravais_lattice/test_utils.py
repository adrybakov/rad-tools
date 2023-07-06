from radtools.crystal.bravais_lattice.utils import (
    bravais_lattice_from_param,
    bravais_lattice_from_cell,
    lattice_example,
)

from radtools.routines import volume

from radtools.crystal.constants import BRAVAIS_LATTICE_VARIATIONS

import pytest
from hypothesis import given, example, strategies as st
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
@given(
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
    st.floats(),
)
def test_bravais_lattice_from_cell(a1, a2, a3, b1, b2, b3, c1, c2, c3):
    cell = [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]
    print(cell)
    if (
        a1**2 + a2**2 + a3**2 != 0
        and b1**2 + b2**2 + b3**2 != 0
        and c1**2 + c2**2 + c3**2 != 0
        and np.linalg.det(cell) != 0
    ):
        cell = [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]
        new_cell = bravais_lattice_from_cell(cell)
        assert np.allclose(volume(cell), volume(new_cell))
