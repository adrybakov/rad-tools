import os

import pytest
import numpy as np

from rad_tools.exchange.model import *
from rad_tools.exchange.template import *


class TestBond:

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
        assert bond.dis == 0

    def test_matrix(self):
        bond = Bond()
        with pytest.raises(ValueError):
            bond.matrix = 1
        bond.matrix = [[1, 2, 0], [1, 1, 0], [0, 0, 1]]
        assert (bond.matrix == np.array([[1, 2, 0],
                                         [1, 1, 0],
                                         [0, 0, 1]])).all()

    def test_symm_assym_matrix(self):
        bond = Bond(matrix=[[1, 2, 0], [1, 1, 0], [0, 0, 1]])
        assert (bond.symm_matrix == np.array([[1, 1.5, 0],
                                              [1.5, 1, 0],
                                              [0, 0, 1]])).all()

        assert (bond.asymm_matrix == np.array([[0, 0.5, 0],
                                              [-0.5, 0, 0],
                                              [0, 0, 0]])).all()

    def test_iso(self):
        bond = Bond()
        bond.iso = 23
        assert bond.iso == 23
        bond.iso = None
        assert bond.iso == 0

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

    def test_addition_substruction(self):
        bond1 = Bond(iso=1, distance=3.45)
        bond2 = Bond(iso=2, distance=4.32)
        bond3 = Bond(iso=3, distance=3.45)
        with pytest.raises(ValueError):
            bond = bond1 + bond2
        with pytest.raises(ValueError):
            bond = bond1 - bond2
        with pytest.raises(TypeError):
            bond = bond1 + 1
        with pytest.raises(TypeError):
            bond = 1 + bond1
        bond = bond1 + bond3
        assert bond.iso == 4
        bond -= bond1
        assert bond.iso == 3

    def test_multiplication_division(self):
        bond1 = Bond(iso=1, distance=3.45)
        bond = bond1 * 5
        assert bond.iso == 5
        assert bond.dis == 3.45
        bond = bond / 2
        assert bond.iso == 2.5
        assert bond.dis == 3.45
        bond = 5 * bond1
        assert bond.iso == 5
        assert bond.dis == 3.45


class TestExchangeModel:

    def test_cell(self):
        model = ExchangeModel()
        with pytest.raises(ValueError):
            model.cell = [3, 4, 5]

    def test_abc(self):
        model = ExchangeModel()
        model.a = [1, 2, 3]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.zeros(3)).all()
        assert (model.c == np.zeros(3)).all()
        assert (model.cell == np.array([[1, 2, 3],
                                        [0, 0, 0],
                                        [0, 0, 0]])).all()
        model.b = [4, 5, 6]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([4, 5, 6])).all()
        assert (model.c == np.zeros(3)).all()
        assert (model.cell == np.array([[1, 2, 3],
                                        [4, 5, 6],
                                        [0, 0, 0]])).all()
        model.c = [7, 8, 9]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([4, 5, 6])).all()
        assert (model.c == np.array([7, 8, 9])).all()
        assert (model.cell == np.array([[1, 2, 3],
                                        [4, 5, 6],
                                        [7, 8, 9]])).all()

    def test_add_bond(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 2, 5, 1)
        model.add_atom("Cr2", 4, 2, 1)
        model.add_atom("Cr3", 5, 1, 8)
        bond12 = Bond(iso=12)
        bond13 = Bond(iso=13)
        bond23 = Bond(iso=23)
        bond31 = Bond(iso=31)
        model.add_bond(bond12, "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(bond23, "Cr2", "Cr3", (0, 0, 0))
        model.add_bond(bond31, "Cr3", "Cr1", (0, 0, 0))
        assert len(model.bonds) == 3
        assert len(model.bonds["Cr1"]) == 1
        assert len(model.bonds["Cr2"]) == 1
        assert len(model.bonds["Cr3"]) == 1
        model.add_bond(bond13, "Cr1", "Cr3", (0, 0, 0))
        assert len(model.bonds) == 3
        assert len(model.bonds["Cr1"]) == 2
        assert len(model.bonds["Cr2"]) == 1
        assert len(model.bonds["Cr3"]) == 1

    def test_remove_bond(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 2, 5, 1)
        model.add_atom("Cr2", 4, 2, 1)
        model.add_atom("Cr3", 5, 1, 8)
        bond12 = Bond(iso=12)
        bond13 = Bond(iso=13)
        bond23 = Bond(iso=23)
        bond31 = Bond(iso=31)
        model.add_bond(bond12, "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(bond23, "Cr2", "Cr3", (0, 0, 0))
        model.add_bond(bond31, "Cr3", "Cr1", (0, 0, 0))
        model.add_bond(bond13, "Cr1", "Cr3", (0, 0, 0))
        assert len(model.bonds) == 3
        assert len(model.bonds["Cr1"]) == 2
        assert len(model.bonds["Cr2"]) == 1
        assert len(model.bonds["Cr3"]) == 1
        model.remove_bond("Cr1", "Cr2", (0, 0, 0))
        assert len(model.bonds) == 3
        assert len(model.bonds["Cr1"]) == 1
        assert len(model.bonds["Cr2"]) == 1
        assert len(model.bonds["Cr3"]) == 1
        model.remove_bond("Cr3", "Cr1", (0, 0, 0))
        assert len(model.bonds) == 2
        assert len(model.bonds["Cr1"]) == 1
        assert len(model.bonds["Cr2"]) == 1

    def test_bond_list(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        model.add_atom("Cr2", 1, 3, 5)
        model.add_atom("Cr3", 1, 3, 3)
        model.add_bond(Bond(iso=1), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2), "Cr1", "Cr3", (0, -1, 0))
        model.add_bond(Bond(iso=3), "Cr2", "Cr1", (0, 0, -3))
        bond_list = model.bond_list
        assert len(bond_list) == 3
        assert set(bond_list) == {("Cr1", "Cr2", (0, 0, 0)),
                                  ("Cr1", "Cr3", (0, -1, 0)),
                                  ("Cr2", "Cr1", (0, 0, -3))}

    def test_cell_list(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        model.add_atom("Cr2", 1, 3, 5)
        model.add_atom("Cr3", 1, 3, 3)
        model.add_bond(Bond(iso=1), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2), "Cr1", "Cr3", (0, -1, 0))
        model.add_bond(Bond(iso=3), "Cr2", "Cr1", (0, 0, -3))
        cells = model.cell_list
        assert len(cells) == 3
        assert set(cells) == {(0, 0, 0), (0, -1, 0), (0, 0, -3)}

    def test_add_atom(self):
        model = ExchangeModel()
        assert len(model.magnetic_atoms) == 0
        model.add_atom("Cr1", 2, 5, 1)
        assert len(model.magnetic_atoms) == 1
        model.add_atom("Cr2", 4, 2, 1)
        assert len(model.magnetic_atoms) == 2
        assert model.magnetic_atoms["Cr1"][0] == 2
        assert model.magnetic_atoms["Cr1"][1] == 5
        assert model.magnetic_atoms["Cr1"][2] == 1
        assert model.magnetic_atoms["Cr2"][0] == 4
        assert model.magnetic_atoms["Cr2"][1] == 2
        assert model.magnetic_atoms["Cr2"][2] == 1
        model.add_atom("Cr2", 4, 3, 1)
        assert model.magnetic_atoms["Cr2"][0] == 4
        assert model.magnetic_atoms["Cr2"][1] == 3
        assert model.magnetic_atoms["Cr2"][2] == 1

    def test_remove_atom(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 2, 5, 1)
        model.add_atom("Cr2", 4, 2, 1)
        model.add_atom("Cr3", 5, 1, 8)
        bond12 = Bond(iso=12)
        bond13 = Bond(iso=13)
        bond23 = Bond(iso=23)
        bond31 = Bond(iso=31)
        model.add_bond(bond12, "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(bond23, "Cr2", "Cr3", (0, 0, 0))
        model.add_bond(bond31, "Cr3", "Cr1", (0, 0, 0))
        model.add_bond(bond13, "Cr1", "Cr3", (0, 0, 0))
        assert len(model.bonds) == 3
        assert len(model.bonds["Cr1"]) == 2
        assert len(model.bonds["Cr2"]) == 1
        assert len(model.bonds["Cr3"]) == 1
        assert len(model.magnetic_atoms) == 3
        model.remove_atom("Cr1")
        assert len(model.magnetic_atoms) == 2
        assert len(model.bonds) == 1
        assert len(model.bonds["Cr2"]) == 1

    @pytest.mark.parametrize("a, b, c, A, B, C, volume", [
        ((1, 0, 0), (0, 1, 0), (0, 0, 1), 1, 1, 1, 1),
        ((3.588, 0, 0), (0, 4.807, 0), (0, 0, 23.571),
         3.588, 4.807, 23.571, 406.5412),
        ((4, 3, 0), (0, 0, 0), (0, 0, 0), 5, 0, 0, 0),
        ((0, 0, 0), (4, 3, 0), (0, 0, 0), 0, 5, 0, 0),
        ((0, 0, 0), (0, 0, 0), (4, 3, 0), 0, 0, 5, 0),
    ])
    def test_len_abc_and_volume(self, a, b, c, A, B, C, volume):
        model = ExchangeModel()
        model.cell = np.array([a, b, c], dtype=float)
        assert A == round(model.len_a, 4)
        assert B == round(model.len_b, 4)
        assert C == round(model.len_c, 4)
        assert volume == round(model.unit_cell_volume, 4)

    def test_get_atom_coordinates(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        x, y, z = model.get_atom_coordinates("Cr1")
        assert x == 1
        assert y == 6
        assert z == 2
        x, y, z = model.get_atom_coordinates("Cr1", R=[1, 0, 0])
        assert x == 11
        assert y == 6
        assert z == 2
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -1, 0])
        assert x == 1
        assert y == -4
        assert z == 2
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, 0, 2])
        assert x == 1
        assert y == 6
        assert z == 22
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -3, 2])
        assert x == 1
        assert y == -24
        assert z == 22
        x, y, z = model.get_atom_coordinates("Cr1", R=[3, -3, 2])
        assert x == 31
        assert y == -24
        assert z == 22
        model.a = [10, 10, 0]
        model.b = [0, 10, 10]
        model.c = [0, 0, 10]
        x, y, z = model.get_atom_coordinates("Cr1")
        assert x == 1
        assert y == 6
        assert z == 2
        x, y, z = model.get_atom_coordinates("Cr1", R=[1, 0, 0])
        assert x == 11
        assert y == 16
        assert z == 2
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -1, 0])
        assert x == 1
        assert y == -4
        assert z == -8
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, 0, 2])
        assert x == 1
        assert y == 6
        assert z == 22
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -3, 2])
        assert x == 1
        assert y == -24
        assert z == -8
        x, y, z = model.get_atom_coordinates("Cr1", R=[3, -3, 2])
        assert x == 31
        assert y == 6
        assert z == -8

    def test_get_bond_coordinate(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        model.add_atom("Cr2", 1, 3, 4)
        x, y, z = model.get_bond_coordinate("Cr1", "Cr2")
        assert x == 1
        assert y == 4.5
        assert z == 3
        x, y, z = model.get_bond_coordinate("Cr1", "Cr2", R=[3, 2, -5])
        assert x == 16
        assert y == 14.5
        assert z == -22

    def test_space_dimensions(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        model.add_atom("Cr2", 1, 3, 5)
        model.add_atom("Cr3", 1, 3, 3)
        model.add_bond(Bond(iso=1), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2), "Cr1", "Cr3", (0, -1, 0))
        model.add_bond(Bond(iso=3), "Cr2", "Cr1", (0, 0, -3))
        x_min, y_min, z_min, x_max, y_max, z_max = model.space_dimensions
        assert x_min == 1
        assert y_min == -7
        assert z_min == -28
        assert x_max == 1
        assert y_max == 6
        assert z_max == 5


class TestExchangeModelTB2J:

    model = ExchangeModelTB2J(os.path.join(
        'utest', 'exchange', 'resourses', 'exchange.out'
    ))
    template = [('Cr1', 'Cr1', (1, 0, 0)), ('Cr1', 'Cr1', (-1, 0, 0)),
                ('Cr2', 'Cr2', (1, 0, 0)), ('Cr2', 'Cr2', (-1, 0, 0)),
                ('Cr1', 'Cr1', (1, 0, 0)), ('Cr2', 'Cr2', (-1, 0, 0)),
                ('Cr1', 'Cr1', (-1, 0, 0)), ('Cr2', 'Cr2', (1, 0, 0)),
                ('Cr1', 'Cr2', (0, 0, 0)), ('Cr1', 'Cr2', (1, 0, 0)),
                ('Cr1', 'Cr2', (0, -1, 0)), ('Cr1', 'Cr2', (1, -1, 0)),
                ('Cr2', 'Cr1', (0, 0, 0)), ('Cr2', 'Cr1', (-1, 0, 0)),
                ('Cr2', 'Cr1', (0, 1, 0)), ('Cr2', 'Cr1', (-1, 1, 0)),
                ('Cr1', 'Cr2', (1, 0, 0)), ('Cr1', 'Cr2', (1, -1, 0)),
                ('Cr2', 'Cr1', (-1, 0, 0)), ('Cr2', 'Cr1', (-1, 1, 0)),
                ('Cr1', 'Cr2', (0, 0, 0)), ('Cr1', 'Cr2', (0, -1, 0)),
                ('Cr2', 'Cr1', (0, 0, 0)), ('Cr2', 'Cr1', (0, 1, 0)),
                ('Cr1', 'Cr1', (0, 1, 0)), ('Cr1', 'Cr1', (0, -1, 0)),
                ('Cr2', 'Cr2', (0, 1, 0)), ('Cr2', 'Cr2', (0, -1, 0))]


class TestInputFilename(TestExchangeModelTB2J):

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            model = ExchangeModelTB2J(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            model = ExchangeModelTB2J(
                "Ah, music. A magic beyond all we do here!")

    def test_correct_filename(self):
        model = ExchangeModelTB2J(os.path.join(
            'utest', 'exchange', 'resourses', 'exchange.out'))


class TestReadFunctions(TestExchangeModelTB2J):

    def test_read_cell(self):
        assert self.model.cell is not None
        cell_values = [[3.588, 0.000, 0.000],
                       [0.000,  4.807,  0.000],
                       [0.000,  0.000, 23.571]]
        for i in range(0, 3):
            for j in range(0, 3):
                assert self.model.cell[i][j] == cell_values[i][j]

    def test_read_atoms(self):
        assert self.model._atoms is not None
        atoms_value = {
            'Br1': (0.8970, 1.2018, -0.0668),
            'Cr1': (2.6910, 1.2018,  1.7371),
            'S1':  (2.6910, 3.6054,  2.2030),
            'Br2': (2.6910, 3.6054,  5.6376),
            'Cr2': (0.8970, 3.6054,  3.8336),
            'S2':  (0.8970, 1.2018,  3.3678)
        }
        for key in self.model._atoms:
            assert key in atoms_value
        for key in atoms_value:
            assert key in self.model._atoms
            assert atoms_value[key] == self.model._atoms[key]

    @ pytest.mark.parametrize("atom1, atom2, R, iso, aniso, dmi, distance",
                              [
                                  ('Cr1', 'Cr1',
                                   (-1, 0, 0),
                                   3.5006,
                                   [[0.006, 0, 0],
                                    [0, -0.016, 0],
                                    [0, 0, 0.01]],
                                   (0, -0.0163, 0),
                                   3.588),

                                  ('Cr1', 'Cr2',
                                      (0, 0, 0),
                                      3.0733,
                                      [[0.0007, 0, 0.004],
                                       [0, -0.0033, 0],
                                          [0.004, 0, 0.0027]],
                                      (-0.0002, 0.0001, 0.0001),
                                      3.659),

                                  ('Cr1', 'Cr1',
                                      (0, 1, 0),
                                      4.1389,
                                      [[-0.007, 0, 0],
                                       [0, 0.006, 0],
                                          [0, 0, 0.001]],
                                      (-0.0584, 0, 0),
                                      4.807),

                                  ('Cr1', 'Cr1',
                                      (-1, 1, 0),
                                      0.0362,
                                      [[-0.001, 0, 0],
                                       [0, 0, 0],
                                          [0, 0, 0.001]],
                                      (0.0194, -0.0170, 0),
                                      5.999),

                                  ('Cr2', 'Cr2',
                                      (0, -1, 0),
                                      4.1333,
                                      [[-0.007, 0, 0],
                                       [0, 0.006, 0],
                                          [0, 0, 0.001]],
                                      (-0.0568, 0, 0),
                                      4.807),

                                  ('Cr2', 'Cr2',
                                      (0, 2, 0),
                                      0.1206,
                                      [[-0.0007, 0, 0],
                                       [0, 0.0003, 0],
                                          [0, 0, 0.0003]],
                                      (0.0363, 0, 0),
                                      9.614),

                                  ('Cr2', 'Cr1',
                                      (1, -1, 0),
                                      0.0028,
                                      [[0, 0, 0],
                                       [0, 0, 0],
                                          [0, 0, 0]],
                                      (0, 0, 0),
                                      9.239),

                                  ('Cr1', 'Cr2',
                                      (1, -2, 0),
                                      0.549,
                                      [[0.0003, 0, 0],
                                       [0, -0.0007, 0],
                                          [0, 0, 0.0003]],
                                      (0.0001, -0.0001, 0),
                                      7.721),

                                  ('Cr2', 'Cr1',
                                      (-1, 1, 0),
                                      3.0733,
                                      [[0.0007, 0, -0.004],
                                       [0, -0.0033, 0],
                                          [-0.004, 0, 0.0027]],
                                      (-0.0002, 0.0001, -0.0001),
                                      3.659),
                              ])
    def test_read_exchange_examples(self, atom1, atom2, R, iso, aniso, dmi, distance):
        assert round(self.model.bonds[atom1][atom2][R].iso, 4) == iso
        assert round(self.model.bonds[atom1][atom2][R].dis, 4) == distance
        for i in range(0, 3):
            assert round(self.model.bonds[atom1][atom2][R].dmi[i], 4) == dmi[i]
            for j in range(0, 3):
                assert round(
                    self.model.bonds[atom1][atom2][R].aniso[i][j], 4) == aniso[i][j]

    def test_magnetic_atoms(self):
        assert len(self.model.magnetic_atoms) == 2
        assert (self.model.magnetic_atoms['Cr1'] == (
            2.6910, 1.2018,  1.7371)).all()
        assert (self.model.magnetic_atoms['Cr2'] == (
            0.8970, 3.6054,  3.8336)).all()


class TestFilter(TestExchangeModelTB2J):

    def count_entries(self, dictionary):
        i = 0
        if dictionary is not None:
            for atom1 in dictionary:
                for atom2 in dictionary[atom1]:
                    for R in dictionary[atom1][atom2]:
                        i += 1
        return i

    def test_instance_data_type(self):
        filtered_model = self.model.filter(max_distance=5)
        assert type(self.model) == type(filtered_model)

    @ pytest.mark.parametrize("max_distance, elements_number", [
        (4.807, 16),
        (6, 24),
        (0, 0),
        (5, 16),
        (4.0, 12),
        (1000, 72)
    ])
    def test_filter_by_max_distance(self, max_distance, elements_number):
        filtered_model = self.model.filter(max_distance=max_distance)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("min_distance, elements_number", [
        (0, 72),
        (6, 48),
        (9.6, 4),
        (10, 0)
    ])
    def test_filter_by_min_distance(self, min_distance, elements_number):
        filtered_model = self.model.filter(min_distance=min_distance)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("R_vector, elements_number", [
        ((0, 0, 0), 2),
        ([(0, 0, 0), (1, 0, 0)], 6),
        ([(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 10),
        ([], 0),
        ([(0, 0, 0), (1, 0, 0), (1, 0, 0)], 6)
    ])
    def test_filter_by_R_vector(self, R_vector, elements_number):
        filtered_model = self.model.filter(R_vector=R_vector)
        assert self.count_entries(filtered_model.bonds) == elements_number

    def test_filter_by_template(self):
        filtered_model = self.model.filter(
            template=self.template)
        assert self.count_entries(filtered_model.bonds) == 16

    @ pytest.mark.parametrize("min_distance, max_distance, elements_number", [
        (4.807, 4.807, 4),
        (3.6, 5, 12),
        (0, 4.807, 16),
        (0, 10, 72),
        (4, 10, 60)
    ])
    def test_filter_by_max_distance_and_min_distance(self,
                                                     min_distance,
                                                     max_distance,
                                                     elements_number):
        filtered_model = self.model.filter(max_distance=max_distance,
                                           min_distance=min_distance)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("max_distance, elements_number", [
        (4.807, 16),
        (6, 16),
        (0, 0),
        (5, 16),
        (4.0, 12),
        (1000, 16),
        (8, 16)
    ])
    def test_filter_by_max_distance_and_template(self,
                                                 max_distance,
                                                 elements_number):
        filtered_model = self.model.filter(max_distance=max_distance,
                                           template=self.template)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("min_distance, elements_number", [
        (0, 16),
        (3.6, 12),
        (6, 0),
        (9.6, 0),
        (10, 0)
    ])
    def test_filter_by_min_distance_and_template(self,
                                                 min_distance,
                                                 elements_number):
        filtered_model = self.model.filter(min_distance=min_distance,
                                           template=self.template)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("R_vector, elements_number", [
        ((0, 0, 0), 2),
        ([(0, 0, 0), (1, 0, 0)], 5),
        ([(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 8),
        ([], 0),
        ([(0, 0, 0), (1, 0, 0), (1, 0, 0)], 5)
    ])
    def test_filter_by_R_vector_and_template(self, R_vector, elements_number):
        filtered_model = self.model.filter(
            R_vector=R_vector, template=self.template)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("max_distance, R_vector, elements_number", [
        (0, (0, 0, 0), 0),
        (10, (0, 0, 0), 2),
        (5, [(0, 0, 0), (1, 0, 0)], 5),
        (10, [(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 10),
        (5, [(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 8),
        (100, [], 0)
    ])
    def test_filter_by_everything(self, max_distance, R_vector, elements_number):
        filtered_model = self.model.filter(
            R_vector=R_vector, max_distance=max_distance)
        assert self.count_entries(filtered_model.bonds) == elements_number
