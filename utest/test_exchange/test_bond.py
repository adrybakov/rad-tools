from math import sqrt

import pytest

from rad_tools.exchange.bond import *


class TestBond:

    bond = Bond([
        [1, 5, 2],
        [5, 8, 4],
        [2, 6, 3]])

    def test_init(self):
        bond = Bond(iso=23)
        assert bond.iso == 23
        assert (bond.dmi == np.zeros(3, dtype=float)).all()
        assert (bond.aniso == np.zeros((3, 3), dtype=float)).all()
        bond.dmi = (1, 1, 1)
        assert bond.iso == 23
        assert (bond.dmi == np.ones(3, dtype=float)).all()
        assert (bond.aniso == np.zeros((3, 3), dtype=float)).all()
        bond.aniso = [[1, 1, 0],
                      [1, -0.5, 0],
                      [0, 0, -0.5]]
        assert bond.iso == 23
        assert (bond.dmi == np.ones(3, dtype=float)).all()
        assert (bond.aniso == np.array([[1, 1, 0],
                                        [1, -0.5, 0],
                                        [0, 0, -0.5]])).all()
        bond = Bond(iso=23, matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        assert bond.iso == 1

    def test_matrix(self):
        bond = Bond()
        with pytest.raises(ValueError):
            bond.matrix = 1
        bond.matrix = [[1, 2, 0], [1, 1, 0], [0, 0, 1]]
        assert (bond.matrix == np.array([[1, 2, 0],
                                         [1, 1, 0],
                                         [0, 0, 1]])).all()
        assert (self.bond.matrix == np.array([[1, 5, 2],
                                              [5, 8, 4],
                                              [2, 6, 3]])).all()

    def test_symm_assym_matrix(self):
        bond = Bond(matrix=[[1, 2, 0], [1, 1, 0], [0, 0, 1]])
        assert (bond.symm_matrix == np.array([[1, 1.5, 0],
                                              [1.5, 1, 0],
                                              [0, 0, 1]])).all()

        assert (bond.asymm_matrix == np.array([[0, 0.5, 0],
                                               [-0.5, 0, 0],
                                               [0, 0, 0]])).all()

        assert (self.bond.symm_matrix == np.array([[1, 5, 2],
                                                   [5, 8, 5],
                                                   [2, 5, 3]])).all()

        assert (self.bond.asymm_matrix == np.array([[0, 0, 0],
                                                    [0, 0, -1],
                                                    [0, 1, 0]])).all()

    def test_iso(self):
        bond = Bond()
        bond.iso = 23
        assert bond.iso == 23
        bond.iso = None
        assert bond.iso == 0
        assert self.bond.iso == 4

    def test_iso_matrix(self):
        bond = Bond()
        bond.iso = 23
        assert (bond.iso_matrix == np.array([[23, 0, 0],
                                             [0, 23, 0],
                                             [0, 0, 23]])).all()
        bond.iso = None
        assert (bond.iso_matrix == np.array([[0, 0, 0],
                                             [0, 0, 0],
                                             [0, 0, 0]])).all()
        assert (self.bond.iso_matrix == np.array([[4, 0, 0],
                                                  [0, 4, 0],
                                                  [0, 0, 4]])).all()

    def test_aniso(self):
        bond = Bond()
        bond.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (bond.aniso == np.array([[1, 1, 0],
                                        [1, -0.5, 0],
                                        [0, 0, -0.5]])).all()
        bond.iso = 23
        assert (bond.aniso == np.array([[1, 1, 0],
                                        [1, -0.5, 0],
                                        [0, 0, -0.5]])).all()
        bond.aniso = None
        assert (bond.aniso == np.zeros((3, 3))).all()
        with pytest.raises(ValueError):
            bond.aniso = 1

        assert (self.bond.aniso == np.array([[-3, 5, 2],
                                             [5, 4, 5],
                                             [2, 5, -1]])).all()

    def test_aniso_diagonal(self):
        assert (self.bond.aniso_diagonal == np.array([-3, 4, -1])).all()

    def test_aniso_diagonal_matrix(self):
        assert (self.bond.aniso_diagonal_matrix == np.array([[-3, 0, 0],
                                                             [0, 4, 0],
                                                             [0, 0, -1]])).all()

    def test_dmi(self):
        bond = Bond()
        bond.dmi = (1, 2, 3)
        assert (bond.dmi == np.array([1, 2, 3])).all()
        bond.iso = 23
        bond.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (bond.dmi == np.array([1, 2, 3])).all()
        bond.dmi = None
        assert (bond.dmi == np.zeros(3)).all()
        with pytest.raises(ValueError):
            bond.dmi = 1
        assert (self.bond.dmi == np.array([-1,  0, 0])).all()

    def test_dmi_matrix(self):
        bond = Bond()
        bond.dmi = (1, 2, 3)
        assert (bond.dmi_matrix == np.array([[0, 3, -2],
                                             [-3, 0, 1],
                                             [2, -1, 0]])).all()
        bond.iso = 23
        bond.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (bond.dmi_matrix == np.array([[0, 3, -2],
                                             [-3, 0, 1],
                                             [2, -1, 0]])).all()
        bond.dmi = None
        assert (bond.dmi_matrix == np.array([[0, 0, 0],
                                             [0, 0, 0],
                                             [0, 0, 0]])).all()
        assert (self.bond.dmi_matrix == np.array([[0, 0, 0],
                                                  [0, 0, -1],
                                                  [0, 1, 0]])).all()

    def test_dmi_module(self):
        bond = Bond()
        bond.dmi = (1, 2, 3)
        assert bond.dmi_module == sqrt(14)
        bond.iso = 23
        bond.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert bond.dmi_module == sqrt(14)
        bond.dmi = None
        assert bond.dmi_module == 0
        assert self.bond.dmi_module == 1

    def test_round(self):
        bond = self.bond / 13
        bond.round(5)
        assert bond.iso == 0.30769
        bond.round(1)
        assert bond.iso == 0.3

    def test_addition_subtraction(self):
        bond1 = Bond(iso=1)
        bond2 = Bond(iso=2)
        bond3 = Bond(iso=3)
        with pytest.raises(TypeError):
            bond = bond1 + 1
        with pytest.raises(TypeError):
            bond = 1 + bond1
        bond = bond1 + bond3
        assert bond.iso == 4
        bond -= bond1
        assert bond.iso == 3

    def test_multiplication_division(self):
        bond1 = Bond(iso=1)
        bond = bond1 * 5
        assert bond.iso == 5
        bond = bond / 2
        assert bond.iso == 2.5
        bond = 5 * bond1
        assert bond.iso == 5
