import pytest
import numpy as np

from rad_tools.exchange.model import *


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
        assert (model.b == np.array([0, 1, 0])).all()
        assert (model.c == np.array([0, 0, 1])).all()
        assert (model.cell == np.array([[1, 2, 3],
                                        [0, 1, 0],
                                        [0, 0, 1]])).all()
        model.b = [4, 5, 6]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([4, 5, 6])).all()
        assert (model.c == np.array([0, 0, 1])).all()
        assert (model.cell == np.array([[1, 2, 3],
                                        [4, 5, 6],
                                        [0, 0, 1]])).all()
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

    @ pytest.mark.parametrize("a, b, c, A, B, C, volume", [
        ((1, 0, 0), (0, 1, 0), (0, 0, 1), 1, 1, 1, 1),
        ((3.588, 0, 0), (0, 4.807, 0), (0, 0, 23.571),
         3.588, 4.807, 23.571, 406.5412),
        ((4, 3, 0), (0, 1, 0), (0, 0, 1), 5, 1, 1, 4),
        ((1, 0, 0), (4, 3, 0), (0, 0, 1), 1, 5, 1, 3),
        ((1, 0, 0), (0, 1, 0), (0, 4, 3), 1, 1, 5, 3),
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

    def test_get_bond_coordinates(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 1, 6, 2)
        model.add_atom("Cr2", 1, 3, 4)
        x, y, z = model.get_bond_coordinates("Cr1", "Cr2")
        assert x == 1
        assert y == 4.5
        assert z == 3
        x, y, z = model.get_bond_coordinates("Cr1", "Cr2", R=[3, 2, -5])
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

    def test_remove_double_bonds(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (1, 0, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (-1, 0, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr1", "Cr2", (1, 0, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (-1, 1, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (1, -1, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr2", "Cr1", (1, -1, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 1, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, -1, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr1", "Cr2", (0, 1, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 0, -2))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, 0, 2))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr2", "Cr1", (0, 0, 2))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, 0, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr2", "Cr1", (0, 0, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.3, 0.2)
        model.add_atom("Cr2", 0.1, 0.6, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, 0, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr1", "Cr2", (0, 0, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.2, 0.3, 0.2)
        model.add_atom("Cr2", 0.1, 0.6, 0.5)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, 0, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr2", "Cr1", (0, 0, 0))

        model = ExchangeModel()
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.6, 0.1)
        model.add_bond(Bond(iso=1, dmi=(1, 2, 3)), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2, dmi=(-1, -2, -3)), "Cr2", "Cr1", (0, 0, 0))
        model.remove_double_bonds()
        assert len(model.bond_list) == 1
        assert model.bond_list[0] == ("Cr2", "Cr1", (0, 0, 0))

    def test_filter(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        model.add_bond(Bond(iso=12, distance=3), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=12, distance=3), "Cr2", "Cr1", (0, 0, 0))
        model.add_bond(Bond(iso=12, distance=4), "Cr1", "Cr1", (1, 0, 0))
        model.add_bond(Bond(iso=12, distance=4), "Cr1", "Cr1", (-1, 0, 0))
        model.add_bond(Bond(iso=12, distance=4), "Cr2", "Cr2", (1, 0, 0))
        model.add_bond(Bond(iso=12, distance=4), "Cr2", "Cr2", (-1, 0, 0))
        model.add_bond(Bond(iso=12, distance=5), "Cr1", "Cr1", (0, 1, 0))
        model.add_bond(Bond(iso=12, distance=5), "Cr1", "Cr1", (0, -1, 0))
        model.add_bond(Bond(iso=12, distance=5), "Cr2", "Cr2", (0, 1, 0))
        model.add_bond(Bond(iso=12, distance=5), "Cr2", "Cr2", (0, -1, 0))
        model.add_bond(Bond(iso=12, distance=6), "Cr2", "Cr1", (1, 1, 0))
        model.add_bond(Bond(iso=12, distance=6), "Cr1", "Cr2", (-1, -1, 0))
        assert len(model.bond_list) == 12
        filtered_model = model.filtered(max_distance=4)
        assert len(filtered_model.bond_list) == 6
        filtered_model = model.filtered(min_distance=4)
        assert len(filtered_model.bond_list) == 10
        filtered_model = model.filtered(min_distance=4, max_distance=5)
        assert len(filtered_model.bond_list) == 8
        filtered_model = model.filtered(R_vector=(0, 0, 0))
        assert len(filtered_model.bond_list) == 2
        filtered_model = model.filtered(R_vector=[(0, 0, 0), (1, 0, 0)])
        assert len(filtered_model.bond_list) == 4
        filtered_model = model.filtered(template=[("Cr1", "Cr2", (0, 0, 0))])
        assert len(filtered_model.bond_list) == 1
        filtered_model = model.filtered(template=[("Cr1", "Cr2", (0, 0, 0))],
                                        R_vector=[(0, 0, 0), (1, 0, 0)])
        assert len(filtered_model.bond_list) == 1
