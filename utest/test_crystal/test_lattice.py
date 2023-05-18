from math import pi

import pytest

from rad_tools.crystal.lattice import *


class TestLattice:
    l = Lattice([1, 0, 0], [0, 2, 0], [0, 0, 3])

    def test_init(self):
        assert (self.l.cell == np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])).all()
        assert (self.l.a1 == np.array([1, 0, 0])).all()
        assert (self.l.a2 == np.array([0, 2, 0])).all()
        assert (self.l.a3 == np.array([0, 0, 3])).all()
        assert self.l.points == {}
        assert self.l.path == []

    def test_undefined_pearson_symbol(self):
        with pytest.raises(RuntimeError):
            tmp = self.l.pearson_symbol
        with pytest.raises(RuntimeError):
            tmp = self.l.crystal_family
        with pytest.raises(RuntimeError):
            tmp = self.l.centring_type

    def test_unit_cell_volume(self):
        assert self.l.unit_cell_volume == 6

    def test_reciprocal_cell(self):
        assert (
            self.l.reciprocal_cell
            == np.array([[2 * pi, 0, 0], [0, pi, 0], [0, 0, 2 / 3 * pi]])
        ).all()
        assert (self.l.b1 == np.array([2 * pi, 0, 0])).all()
        assert (self.l.b2 == np.array([0, pi, 0])).all()
        assert (self.l.b3 == np.array([0, 0, 2 / 3 * pi])).all()
        assert self.l.reciprocal_cell_volume == 4 / 3 * pi**3
