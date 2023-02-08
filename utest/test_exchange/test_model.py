from math import sqrt, pi

import pytest
import numpy as np

from rad_tools.exchange.model import ExchangeModel
from rad_tools.exchange.bond import Bond
from rad_tools.exchange.template import ExchangeTemplate


class TestExchangeModel:

    def test_iteration(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        bonds = [(12, "Cr1", "Cr2", (0, 0, 0)),
                 (12, "Cr2", "Cr1", (0, 0, 0)),
                 (12, "Cr1", "Cr1", (1, 0, 0)),
                 (12, "Cr1", "Cr1", (-1, 0, 0)),
                 (12, "Cr2", "Cr2", (1, 0, 0)),
                 (12, "Cr2", "Cr2", (-1, 0, 0)),
                 (12, "Cr1", "Cr1", (0, 2, 0)),
                 (12, "Cr1", "Cr1", (0, -2, 0)),
                 (12, "Cr2", "Cr2", (0, 2, 0)),
                 (12, "Cr2", "Cr2", (0, -2, 0)),
                 (12, "Cr2", "Cr1", (2, 2, 0)),
                 (12, "Cr1", "Cr2", (-2, -2, 0))]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(Bond(iso=iso), atom1, atom2, R)
        for atom1, atom2, R in model:
            assert isinstance(atom1, str)
            assert isinstance(atom2, str)
            assert isinstance(R, tuple)

    def test_contains(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        bonds = [(12, "Cr1", "Cr2", (0, 0, 0)),
                 (12, "Cr2", "Cr1", (0, 0, 0)),
                 (12, "Cr1", "Cr1", (1, 0, 0)),
                 (12, "Cr1", "Cr1", (-1, 0, 0)),
                 (12, "Cr2", "Cr2", (1, 0, 0)),
                 (12, "Cr2", "Cr2", (-1, 0, 0)),
                 (12, "Cr1", "Cr1", (0, 2, 0)),
                 (12, "Cr1", "Cr1", (0, -2, 0)),
                 (12, "Cr2", "Cr2", (0, 2, 0)),
                 (12, "Cr2", "Cr2", (0, -2, 0)),
                 (12, "Cr2", "Cr1", (2, 2, 0)),
                 (12, "Cr1", "Cr2", (-2, -2, 0))]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(Bond(iso=iso), atom1, atom2, R)
        assert ("Cr1", "Cr2", (0, 0, 0)) in model
        assert ("Cr2", "Cr2", (1, 0, 0)) in model
        assert ("Cr1", "Cr2", (4, 0, 0)) not in model
        assert ("Cr1", "Cr2", (0, 8, 0)) not in model

    def test_getitem(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        bonds = [(12, "Cr1", "Cr2", (0, 0, 0)),
                 (12, "Cr2", "Cr1", (0, 0, 0)),
                 (12, "Cr1", "Cr1", (1, 0, 0)),
                 (12, "Cr1", "Cr1", (-1, 0, 0)),
                 (12, "Cr2", "Cr2", (1, 0, 0)),
                 (12, "Cr2", "Cr2", (-1, 0, 0)),
                 (12, "Cr1", "Cr1", (0, 2, 0)),
                 (12, "Cr1", "Cr1", (0, -2, 0)),
                 (12, "Cr2", "Cr2", (0, 2, 0)),
                 (12, "Cr2", "Cr2", (0, -2, 0)),
                 (12, "Cr2", "Cr1", (2, 2, 0)),
                 (12, "Cr1", "Cr2", (-2, -2, 0))]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(Bond(iso=iso), atom1, atom2, R)
        assert model[("Cr1", "Cr2", (0, 0, 0))].iso == 12
        assert model[("Cr2", "Cr2", (1, 0, 0))].iso == 12
        assert model[("Cr1", "Cr2", (-2, -2, 0))].iso == 12
        with pytest.raises(KeyError):
            assert model[("Cr1", "Cr2", (0, 8, 0))].iso == 12

    def test_cell_error(self):
        model = ExchangeModel()
        with pytest.raises(ValueError):
            model.cell = [3, 4, 5]
        # with pytest.raises(ValueError):
        #     model.cell = [[3, 4, 5],
        #                   [9, 8, 10],
        #                   [1, 2, 3]]

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
        model.cell = [[2, 5, 3],
                      [7, 3, 1],
                      [8, 4, 9]]
        assert (model.cell == np.array([[2, 5, 3],
                                        [7, 3, 1],
                                        [8, 4, 9]])).all()

    def test_len_abc(self):
        model = ExchangeModel()
        model.a = [1, 2, 3]
        model.b = [4, 5, 6]
        model.c = [7, 8, 9]
        assert model.len_a == sqrt(14)
        assert model.len_b == sqrt(77)
        assert model.len_c == sqrt(194)

    def test_unit_cell_volume(self):
        model = ExchangeModel()
        model.a = [1, 2, 3]
        model.b = [4, 5, 6]
        model.c = [7, 8, 8]
        assert model.unit_cell_volume == np.dot(
            np.array([1, 2, 3]), np.cross(np.array([4, 5, 6]), np.array([7, 8, 8])))

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

    def test_b123(self):
        model = ExchangeModel()
        a = [1, 2, 3]
        b = [4, 5, 6]
        c = [7, 8, 8]
        model.a = a
        model.b = b
        model.c = c
        assert (model.b1 == 2 * pi / 3 * np.cross(b, c)).all()
        assert (model.b2 == 2 * pi / 3 * np.cross(c, a)).all()
        assert (model.b3 == 2 * pi / 3 * np.cross(a, b)).all()

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
        cells = model.cell_list.tolist()
        cells = [tuple(i) for i in cells]
        assert len(cells) == 3
        assert (0, 0, 0) in cells
        assert (0, -1, 0) in cells
        assert (0, 0, -3) in cells

    def test_number_spins_in_unit_cell(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 1, 6, 2)
        assert model.number_spins_in_unit_cell == 1
        model.add_atom("Cr2", 1, 3, 5)
        assert model.number_spins_in_unit_cell == 2
        model.add_atom("Cr3", 1, 3, 3)
        assert model.number_spins_in_unit_cell == 3
        model.add_atom("Cr3", 4, 9, 1)
        assert model.number_spins_in_unit_cell == 3
        model.remove_atom("Cr1")
        assert model.number_spins_in_unit_cell == 2
        model.remove_atom("Cr2")
        assert model.number_spins_in_unit_cell == 1
        model.remove_atom("Cr3")
        assert model.number_spins_in_unit_cell == 0

    def test_space_dimensions(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_atom("Cr3", 0.1, 0.3, 0.3)
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

    def test_round(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.5)
        model.add_atom("Cr3", 0.1, 0.3, 0.3)
        model.add_bond(Bond(iso=1/3), "Cr1", "Cr2", (0, 0, 0))
        model.add_bond(Bond(iso=2/7), "Cr1", "Cr3", (0, -1, 0))
        model.add_bond(Bond(iso=3/11), "Cr2", "Cr1", (0, 0, -3))
        model.round(4)
        assert model[("Cr1", "Cr2", (0, 0, 0))].iso == 0.3333
        assert model[("Cr1", "Cr3", (0, -1, 0))].iso == 0.2857
        assert model[("Cr2", "Cr1", (0, 0, -3))].iso == 0.2727
        model.round(1)
        assert model[("Cr1", "Cr2", (0, 0, 0))].iso == 0.3
        assert model[("Cr1", "Cr3", (0, -1, 0))].iso == 0.3
        assert model[("Cr2", "Cr1", (0, 0, -3))].iso == 0.3

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
        assert ("Cr1", "Cr2", (0, 0, 0)) in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) in model
        model.add_bond(bond13, "Cr1", "Cr3", (0, 0, 0))
        assert len(model.bonds) == 4
        assert ("Cr1", "Cr2", (0, 0, 0)) in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) in model
        assert ("Cr1", "Cr3", (0, 0, 0)) in model

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
        assert len(model.bonds) == 4
        assert ("Cr1", "Cr2", (0, 0, 0)) in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) in model
        assert ("Cr1", "Cr3", (0, 0, 0)) in model
        model.remove_bond("Cr1", "Cr2", (0, 0, 0))
        assert len(model.bonds) == 3
        assert ("Cr1", "Cr2", (0, 0, 0)) not in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) in model
        assert ("Cr1", "Cr3", (0, 0, 0)) in model
        model.remove_bond("Cr3", "Cr1", (0, 0, 0))
        assert len(model.bonds) == 2
        assert ("Cr1", "Cr2", (0, 0, 0)) not in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) not in model
        assert ("Cr1", "Cr3", (0, 0, 0)) in model

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
        assert len(model.bonds) == 4
        assert ("Cr1", "Cr2", (0, 0, 0)) in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) in model
        assert ("Cr1", "Cr3", (0, 0, 0)) in model
        assert len(model.magnetic_atoms) == 3
        model.remove_atom("Cr1")
        assert len(model.magnetic_atoms) == 2
        assert len(model.bonds) == 1
        assert ("Cr1", "Cr2", (0, 0, 0)) not in model
        assert ("Cr2", "Cr3", (0, 0, 0)) in model
        assert ("Cr3", "Cr1", (0, 0, 0)) not in model
        assert ("Cr1", "Cr3", (0, 0, 0)) not in model

    def test_get_atom_coordinates(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
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
        assert y == 7
        assert z == 8
        x, y, z = model.get_atom_coordinates("Cr1", R=[1, 0, 0])
        assert x == 11
        assert y == 17
        assert z == 8
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -1, 0])
        assert x == 1
        assert y == -3
        assert z == -2
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, 0, 2])
        assert x == 1
        assert y == 7
        assert z == 28
        x, y, z = model.get_atom_coordinates("Cr1", R=[0, -3, 2])
        assert x == 1
        assert y == -23
        assert z == -2
        x, y, z = model.get_atom_coordinates("Cr1", R=[3, -3, 2])
        assert x == 31
        assert y == 7
        assert z == -2

    def test_get_bond_centre_coordinates(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.4)
        x, y, z = model.get_bond_centre_coordinates("Cr1", "Cr2")
        assert x == 1
        assert y == 4.5
        assert z == 3
        x, y, z = model.get_bond_centre_coordinates("Cr1", "Cr2", R=[3, 2, -5])
        assert x == 16
        assert y == 14.5
        assert z == -22

    def test_get_bond_vector(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.4)
        vector = model.get_bond_vector("Cr1", "Cr2")
        assert (vector == np.array([0, -3, 2])).all()
        vector = model.get_bond_vector("Cr1", "Cr2", R=[3, 2, -5])
        assert (vector == np.array([30, 17, -48])).all()

    def test_get_distance(self):
        model = ExchangeModel()
        model.a = [10, 0, 0]
        model.b = [0, 10, 0]
        model.c = [0, 0, 10]
        model.add_atom("Cr1", 0.1, 0.6, 0.2)
        model.add_atom("Cr2", 0.1, 0.3, 0.4)
        d = model.get_distance("Cr1", "Cr2")
        assert d == sqrt(13)
        d = model.get_distance("Cr1", "Cr2", R=[3, 2, -5])
        assert d == sqrt(900 + 17**2 + 48**2)

    def test_filter(self):
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        bonds = [(12, "Cr1", "Cr2", (0, 0, 0)),
                 (12, "Cr2", "Cr1", (0, 0, 0)),
                 (12, "Cr1", "Cr1", (1, 0, 0)),
                 (12, "Cr1", "Cr1", (-1, 0, 0)),
                 (12, "Cr2", "Cr2", (1, 0, 0)),
                 (12, "Cr2", "Cr2", (-1, 0, 0)),
                 (12, "Cr1", "Cr1", (0, 2, 0)),
                 (12, "Cr1", "Cr1", (0, -2, 0)),
                 (12, "Cr2", "Cr2", (0, 2, 0)),
                 (12, "Cr2", "Cr2", (0, -2, 0)),
                 (12, "Cr2", "Cr1", (2, 2, 0)),
                 (12, "Cr1", "Cr2", (-2, -2, 0))]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(Bond(iso=iso), atom1, atom2, R)

        assert len(model.bonds) == 12
        filtered_model = model.filtered(max_distance=1)
        assert len(filtered_model.bonds) == 6
        filtered_model = model.filtered(min_distance=1)
        assert len(filtered_model.bonds) == 10
        filtered_model = model.filtered(min_distance=1, max_distance=2)
        assert len(filtered_model.bonds) == 8
        filtered_model = model.filtered(R_vector=(0, 0, 0))
        assert len(filtered_model.bonds) == 2
        filtered_model = model.filtered(R_vector=[(0, 0, 0), (1, 0, 0)])
        assert len(filtered_model.bonds) == 4
        filtered_model = model.filtered(template=[("Cr1", "Cr2", (0, 0, 0))])
        assert len(filtered_model.bonds) == 1
        filtered_model = model.filtered(template=[("Cr1", "Cr2", (0, 0, 0))],
                                        R_vector=[(0, 0, 0), (1, 0, 0)])
        assert len(filtered_model.bonds) == 1

    def test_force_symmetry(self):
        template1 = ExchangeTemplate()
        template2 = ExchangeTemplate()
        template3 = ExchangeTemplate()
        model = ExchangeModel()
        model.add_atom("Cr1", 0.25, 0.25, 0)
        model.add_atom("Cr2", 0.75, 0.75, 0)
        bonds = [([[4, 2, 3],
                   [4, 8, 6],
                   [7, 8, 12]], "Cr1", "Cr2", (0, 0, 0)),
                 ([[1, 4, 7],
                   [2, 5, 8],
                   [3, 6, 9]], "Cr2", "Cr1", (0, 0, 0)),

                 ([[1, 4, 0],
                   [4, 1, -1],
                   [0, 1, 0]], "Cr1", "Cr1", (1, 0, 0)),
                 ([[1, 4, 0],
                   [4, 1, 1],
                   [0, -1, 0]], "Cr1", "Cr1", (-1, 0, 0)),
                 ([[1, 2, 0],
                   [2, 1, -2],
                   [0, 2, 0]], "Cr2", "Cr2", (1, 0, 0)),
                 ([[1, 2, 0],
                   [2, 1, 2],
                   [0, -2, 0]], "Cr2", "Cr2", (-1, 0, 0))]
        for matrix, atom1, atom2, R in bonds:
            model.add_bond(Bond(matrix=matrix), atom1, atom2, R)

        template1.names = {"J1": [("Cr1", "Cr2", (0, 0, 0)),
                                  ("Cr2", "Cr1", (0, 0, 0))],
                           "J2": [("Cr1", "Cr1", (1, 0, 0)),
                                  ("Cr1", "Cr1", (-1, 0, 0)),
                                  ("Cr2", "Cr2", (1, 0, 0)),
                                  ("Cr2", "Cr2", (-1, 0, 0))]}
        template2.names = {"J1": [("Cr1", "Cr2", (0, 0, 0)),
                                  ("Cr2", "Cr1", (0, 0, 0))],
                           "J2": [("Cr1", "Cr1", (1, 0, 0)),
                                  ("Cr1", "Cr1", (-1, 0, 0)),
                                  ("Cr2", "Cr2", (1, 0, 0))]}
        template3.names = {"J1": [("Cr1", "Cr2", (0, 0, 0)),
                                  ("Cr2", "Cr1", (0, 0, 0))],
                           "J2": [("Cr1", "Cr1", (1, 0, 0)),
                                  ("Cr1", "Cr1", (-1, 0, 0))]}

        model.force_symmetry(template=template1)
        assert len(model.bonds) == 6
        assert (model[("Cr1", "Cr2", (0, 0, 0))].matrix == np.array([[2.5, 2, 3],
                                                                     [4, 6.5, 6],
                                                                     [7, 8, 10.5]])).all()
        assert (model[("Cr2", "Cr1", (0, 0, 0))].matrix == np.array([[2.5, 4, 7],
                                                                     [2, 6.5, 8],
                                                                     [3, 6, 10.5]])).all()

        assert (model[("Cr1", "Cr1", (1, 0, 0))].matrix == np.array([[1, 3, 0],
                                                                     [3, 1, -1.5],
                                                                     [0, 1.5, 0]])).all()
        assert (model[("Cr1", "Cr1", (-1, 0, 0))].matrix == np.array([[1, 3, 0],
                                                                      [3, 1,
                                                                       1.5],
                                                                      [0, -1.5, 0]])).all()
        assert (model[("Cr2", "Cr2", (1, 0, 0))].matrix == np.array([[1, 3, 0],
                                                                     [3, 1, -1.5],
                                                                     [0, 1.5, 0]])).all()
        assert (model[("Cr2", "Cr2", (-1, 0, 0))].matrix == np.array([[1, 3, 0],
                                                                      [3, 1,
                                                                       1.5],
                                                                      [0, -1.5, 0]])).all()
        for matrix, atom1, atom2, R in bonds:
            model.add_bond(Bond(matrix=matrix), atom1, atom2, R)
        model.force_symmetry(template=template2)
        assert len(model.bonds) == 5
        assert (model[("Cr1", "Cr2", (0, 0, 0))].matrix == np.array([[2.5, 2, 3],
                                                                     [4, 6.5, 6],
                                                                     [7, 8, 10.5]])).all()
        assert (model[("Cr2", "Cr1", (0, 0, 0))].matrix == np.array([[2.5, 4, 7],
                                                                     [2, 6.5, 8],
                                                                     [3, 6, 10.5]])).all()

        assert (model[("Cr1", "Cr1", (1, 0, 0))].matrix == np.array([[1, 10/3, 0],
                                                                     [10/3,
                                                                      1, -4/3],
                                                                     [0, 4/3, 0]])).all()
        assert (model[("Cr1", "Cr1", (-1, 0, 0))].matrix == np.array([[1, 10/3, 0],
                                                                      [10/3, 1,
                                                                       4/3],
                                                                      [0, -4/3, 0]])).all()
        assert (model[("Cr2", "Cr2", (1, 0, 0))].matrix == np.array([[1, 10/3, 0],
                                                                     [10/3,
                                                                      1, -4/3],
                                                                     [0, 4/3, 0]])).all()
        assert ("Cr2", "Cr2", (-1, 0, 0)) not in model

        for matrix, atom1, atom2, R in bonds:
            model.add_bond(Bond(matrix=matrix), atom1, atom2, R)
        model.force_symmetry(template=template3)
        assert len(model.bonds) == 4
        assert (model[("Cr1", "Cr2", (0, 0, 0))].matrix == np.array([[2.5, 2, 3],
                                                                     [4, 6.5, 6],
                                                                     [7, 8, 10.5]])).all()
        assert (model[("Cr2", "Cr1", (0, 0, 0))].matrix == np.array([[2.5, 4, 7],
                                                                     [2, 6.5, 8],
                                                                     [3, 6, 10.5]])).all()

        assert (model[("Cr1", "Cr1", (1, 0, 0))].matrix == np.array([[1, 4, 0],
                                                                     [4,
                                                                      1, -1],
                                                                     [0, 1, 0]])).all()
        assert (model[("Cr1", "Cr1", (-1, 0, 0))].matrix == np.array([[1, 4, 0],
                                                                      [4, 1,
                                                                       1],
                                                                      [0, -1, 0]])).all()

        assert ("Cr2", "Cr2", (11, 0, 0)) not in model
        assert ("Cr2", "Cr2", (-1, 0, 0)) not in model
