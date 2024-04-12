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

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays as harrays

from radtools import Atom
from radtools.crystal.constants import ATOM_TYPES


@given(
    harrays(
        np.float64,
        (3,),
        elements=st.floats(
            allow_nan=False, allow_infinity=False, max_value=1e9, min_value=-1e9
        ),
    )
)
def test_Atom_spin(v):
    atom = Atom()

    assert np.allclose(atom.spin_direction, [0, 0, 1])
    with pytest.raises(ValueError):
        atom.spin
    with pytest.raises(ValueError):
        atom.spin_vector

    atom.spin = np.linalg.norm(v)
    assert np.allclose(atom.spin, np.linalg.norm(v))
    assert np.allclose(atom.spin_direction, [0, 0, 1])
    assert np.allclose(atom.spin_vector, np.linalg.norm(v) * atom.spin_direction)

    if np.linalg.norm(v) != 0:
        atom.spin_vector = v
        assert np.allclose(atom.spin_direction, v / np.linalg.norm(v))
        assert np.allclose(atom.spin_vector, v)


@given(
    st.text(
        min_size=1,
        max_size=3,
        alphabet=[i for i in "0123456789_-#$%!"],
    ),
    st.text(
        min_size=1,
        max_size=3,
        alphabet=[i for i in "0123456789_-#$%!"],
    ),
)
def test_Atom_type(prefix, suffix):
    if not prefix.startswith("__") and not suffix.endswith("__"):
        for atom_type in ATOM_TYPES:
            atom = Atom(prefix + atom_type + suffix)
            assert atom.type == atom_type
    else:
        with pytest.raises(ValueError):
            for atom_type in ATOM_TYPES:
                atom = Atom(prefix + atom_type + suffix)


def test_atom():
    atom = Atom()
    assert (atom.position == np.zeros(3)).all()
    assert atom.name == "X"
    atom.name = "Cr1"
    assert atom.type == "Cr"
    atom.name = "CR1"
    assert atom.type == "Cr"
    atom.name = "cr1"
    assert atom.type == "Cr"
    with pytest.raises(ValueError):
        a = atom.spin
    with pytest.raises(ValueError):
        a = atom.spin_vector
    with pytest.raises(ValueError):
        a = atom.magmom
    with pytest.raises(ValueError):
        atom.spin_vector = 23
    with pytest.raises(ValueError):
        atom.spin_vector = "adsfasdfs"
    with pytest.raises(ValueError):
        atom.spin_vector = (4, 5)
    with pytest.raises(ValueError):
        atom.spin_vector = (4, 5, 4, 5)
    atom.spin_vector = (0, 0, 1 / 2)
    assert atom.spin - 0.5 < 1e-10
    assert (atom.spin_vector - np.array([0, 0, 0.5]) < 1e-10).all()
    assert (atom.spin_direction - np.array([0, 0, 1]) < 1e-10).all()
    atom.spin = 2
    assert (atom.spin_vector - np.array([0, 0, 2]) < 1e-10).all()
    assert (atom.spin_direction - np.array([0, 0, 1]) < 1e-10).all()
    assert atom.spin - 2 < 1e-10
