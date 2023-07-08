from math import acos, pi, sqrt

import numpy as np
import pytest

from radtools.crystal.lattice import *
from radtools.routines import todegrees
from radtools.crystal.constants import ABS_TOL, REL_TOL


l = Lattice([1, 0, 0], [0, 2, 0], [0, 0, 3])


def test_init():
    pass


def test_extended_init():
    with pytest.raises(ValueError):
        l = Lattice(1, 2)
    l = Lattice(1, 1, 1, 90, 90, 60)
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10
    l = Lattice([1, 0, 0], [1 / 2, sqrt(3) / 2, 0], [0, 0, 1])
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10
    l = Lattice([[1, 0, 0], [1 / 2, sqrt(3) / 2, 0], [0, 0, 1]])
    assert l.a - 1 < 1e-10
    assert l.b - 1 < 1e-10
    assert l.c - 1 < 1e-10
    assert l.alpha - 90 < 1e-10
    assert l.beta - 90 < 1e-10
    assert l.gamma - 60 < 1e-10


def test_cell():
    assert (l.cell == np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])).all()
    assert (l.a1 == np.array([1, 0, 0])).all()
    assert (l.a2 == np.array([0, 2, 0])).all()
    assert (l.a3 == np.array([0, 0, 3])).all()
    assert l.unit_cell_volume == 6
    assert l.a == 1
    assert l.b == 2
    assert l.c == 3
    assert l.alpha == l.beta == l.gamma == 90


def test_reciprocal_cell():
    assert (
        l.reciprocal_cell == np.array([[2 * pi, 0, 0], [0, pi, 0], [0, 0, 2 / 3 * pi]])
    ).all()
    assert np.allclose(l.b1, [2 * pi, 0, 0], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(l.b2, [0, pi, 0], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(l.b3, [0, 0, 2 / 3 * pi], rtol=REL_TOL, atol=ABS_TOL)
    assert np.allclose(
        l.reciprocal_cell_volume, 4 / 3 * pi**3, rtol=REL_TOL, atol=ABS_TOL
    )
    assert np.allclose(
        [l.k_a, l.k_b, l.k_c], [2 * pi, pi, 2 / 3 * pi], rtol=REL_TOL, atol=ABS_TOL
    )

    assert np.allclose(
        [l.k_alpha, l.k_beta, l.k_gamma], [90.0, 90.0, 90.0], rtol=REL_TOL, atol=ABS_TOL
    )


def test_variation():
    assert l.variation == "ORC"
