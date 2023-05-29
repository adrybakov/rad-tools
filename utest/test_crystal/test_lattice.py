from math import pi, sqrt

import numpy as np
import pytest
import numpy as np
from math import acos, sqrt
from rad_tools.crystal.lattice import *
from rad_tools.routines import _todegrees


class TestLattice:
    l = Lattice([1, 0, 0], [0, 2, 0], [0, 0, 3])

    def test_init(self):
        assert self.l.points == {}
        assert self.l.path == []

    def test_extended_init(self):
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

    def test_undefined_pearson_symbol(self):
        with pytest.raises(RuntimeError):
            tmp = self.l.pearson_symbol
        with pytest.raises(RuntimeError):
            tmp = self.l.crystal_family
        with pytest.raises(RuntimeError):
            tmp = self.l.centring_type

    def test_cell(self):
        assert (self.l.cell == np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])).all()
        assert (self.l.a1 == np.array([1, 0, 0])).all()
        assert (self.l.a2 == np.array([0, 2, 0])).all()
        assert (self.l.a3 == np.array([0, 0, 3])).all()
        assert self.l.unit_cell_volume == 6
        assert self.l.a == 1
        assert self.l.b == 2
        assert self.l.c == 3
        assert self.l.alpha == self.l.beta == self.l.gamma == 90

    def test_reciprocal_cell(self):
        assert (
            self.l.reciprocal_cell
            == np.array([[2 * pi, 0, 0], [0, pi, 0], [0, 0, 2 / 3 * pi]])
        ).all()
        assert (self.l.b1 == np.array([2 * pi, 0, 0])).all()
        assert (self.l.b2 == np.array([0, pi, 0])).all()
        assert (self.l.b3 == np.array([0, 0, 2 / 3 * pi])).all()
        assert self.l.reciprocal_cell_volume == 4 / 3 * pi**3
        assert self.l.k_a == 2 * pi
        assert self.l.k_b == pi
        assert self.l.k_c == 2 / 3 * pi
        assert self.l.k_alpha == self.l.k_beta == self.l.k_gamma == 90

    def test_variation(self):
        assert self.l.variation == "Lattice"


def test_get_niggli():
    a = 3
    b = sqrt(27)
    c = 2
    alpha = acos(-5 / 2 / sqrt(27) / 2) * _todegrees
    beta = acos(-4 / 2 / 3 / 2) * _todegrees
    gamma = acos(-22 / 2 / 3 / sqrt(27)) * _todegrees
    assert (
        np.array([[4, 9, 9], [9 / 2, 3 / 2, 2]])
        == get_niggli(a, b, c, alpha, beta, gamma)
    ).all()


def test_niggli_run():
    alpha = 79.030
    beta = 64.130
    gamma = 64.150
    a = 4
    b = 4.472
    c = 4.583
    get_niggli(a, b, c, alpha, beta, gamma)
