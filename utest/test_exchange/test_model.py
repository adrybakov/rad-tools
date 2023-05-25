from math import pi, sqrt

import numpy as np
import pytest

from rad_tools.crystal.atom import Atom
from rad_tools.exchange.parameter import ExchangeParameter
from rad_tools.exchange.model import ExchangeModel, NotationError
from rad_tools.exchange.template import ExchangeTemplate


class TestExchangeModel:
    def test_iteration(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (0.25, 0.25, 0))
        Cr2 = Atom("Cr2", (0.75, 0.75, 0))
        model.add_atom(Cr1)
        model.add_atom(Cr2)
        bonds = [
            (12, Cr1, Cr2, (0, 0, 0)),
            (2, Cr2, Cr1, (0, 0, 0)),
            (6, Cr1, Cr1, (1, 0, 0)),
            (3.425, Cr1, Cr1, (-1, 0, 0)),
            (5.3, Cr2, Cr2, (1, 0, 0)),
            (7.34, Cr2, Cr2, (-1, 0, 0)),
            (12.4, Cr1, Cr1, (0, 2, 0)),
            (34, Cr1, Cr1, (0, -2, 0)),
            (1.098, Cr2, Cr2, (0, 2, 0)),
            (0.0054, Cr2, Cr2, (0, -2, 0)),
            (0.35, Cr2, Cr1, (2, 2, 0)),
            (-2.35, Cr1, Cr2, (-2, -2, 0)),
        ]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(iso=iso), atom1, atom2, R)
        for i, (atom1, atom2, R, J) in enumerate(model):
            assert isinstance(atom1, Atom)
            assert isinstance(atom2, Atom)
            assert isinstance(R, tuple)
            assert J == ExchangeParameter(iso=bonds[i][0])

    def test_contains(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (0.25, 0.25, 0))
        Cr2 = Atom("Cr2", (0.75, 0.75, 0))
        model.add_atom(Cr1)
        model.add_atom(Cr2)
        bonds = [
            (12, Cr1, Cr2, (0, 0, 0)),
            (12, Cr2, Cr1, (0, 0, 0)),
            (12, Cr1, Cr1, (1, 0, 0)),
            (12, Cr1, Cr1, (-1, 0, 0)),
            (12, Cr2, Cr2, (1, 0, 0)),
            (12, Cr2, Cr2, (-1, 0, 0)),
            (12, Cr1, Cr1, (0, 2, 0)),
            (12, Cr1, Cr1, (0, -2, 0)),
            (12, Cr2, Cr2, (0, 2, 0)),
            (12, Cr2, Cr2, (0, -2, 0)),
            (12, Cr2, Cr1, (2, 2, 0)),
            (12, Cr1, Cr2, (-2, -2, 0)),
        ]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(iso=iso), atom1, atom2, R)
        assert (Cr1, Cr2, (0, 0, 0)) in model
        assert (Cr2, Cr2, (1, 0, 0)) in model
        assert (Cr1, Cr2, (4, 0, 0)) not in model
        assert (Cr1, Cr2, (0, 8, 0)) not in model

    def test_getitem(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (0.25, 0.25, 0))
        Cr2 = Atom("Cr2", (0.75, 0.75, 0))
        model.add_atom(Cr1)
        model.add_atom(Cr2)
        bonds = [
            (12, Cr1, Cr2, (0, 0, 0)),
            (12, Cr2, Cr1, (0, 0, 0)),
            (12, Cr1, Cr1, (1, 0, 0)),
            (12, Cr1, Cr1, (-1, 0, 0)),
            (12, Cr2, Cr2, (1, 0, 0)),
            (12, Cr2, Cr2, (-1, 0, 0)),
            (12, Cr1, Cr1, (0, 2, 0)),
            (12, Cr1, Cr1, (0, -2, 0)),
            (12, Cr2, Cr2, (0, 2, 0)),
            (12, Cr2, Cr2, (0, -2, 0)),
            (12, Cr2, Cr1, (2, 2, 0)),
            (12, Cr1, Cr2, (-2, -2, 0)),
        ]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(iso=iso), atom1, atom2, R)
        assert model[(Cr1, Cr2, (0, 0, 0))].iso == 12
        assert model[(Cr2, Cr2, (1, 0, 0))].iso == 12
        assert model[(Cr1, Cr2, (-2, -2, 0))].iso == 12
        with pytest.raises(KeyError):
            assert model[(Cr1, Cr2, (0, 8, 0))].iso == 12

    def test_cell_error(self):
        model = ExchangeModel()
        with pytest.raises(ValueError):
            model.cell = [3, 4, 5]

    def test_abc(self):
        model = ExchangeModel()
        model.cell = [[1, 2, 3], [0, 1, 0], [0, 0, 1]]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([0, 1, 0])).all()
        assert (model.c == np.array([0, 0, 1])).all()
        assert (model.cell == np.array([[1, 2, 3], [0, 1, 0], [0, 0, 1]])).all()
        model.cell = [[1, 2, 3], [4, 5, 6], [0, 0, 1]]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([4, 5, 6])).all()
        assert (model.c == np.array([0, 0, 1])).all()
        assert (model.cell == np.array([[1, 2, 3], [4, 5, 6], [0, 0, 1]])).all()
        model.cell = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        assert (model.a == np.array([1, 2, 3])).all()
        assert (model.b == np.array([4, 5, 6])).all()
        assert (model.c == np.array([7, 8, 9])).all()
        assert (model.cell == np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])).all()
        model.cell = [[2, 5, 3], [7, 3, 1], [8, 4, 9]]
        assert (model.cell == np.array([[2, 5, 3], [7, 3, 1], [8, 4, 9]])).all()

    def test_len_abc(self):
        model = ExchangeModel()
        model.cell = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        assert model.len_a == sqrt(14)
        assert model.len_b == sqrt(77)
        assert model.len_c == sqrt(194)

    def test_unit_cell_volume(self):
        model = ExchangeModel()
        model.cell = [[1, 2, 3], [4, 5, 6], [7, 8, 8]]
        assert model.unit_cell_volume == np.dot(
            np.array([1, 2, 3]), np.cross(np.array([4, 5, 6]), np.array([7, 8, 8]))
        )

    @pytest.mark.parametrize(
        "a, b, c, A, B, C, volume",
        [
            ((1, 0, 0), (0, 1, 0), (0, 0, 1), 1, 1, 1, 1),
            (
                (3.588, 0, 0),
                (0, 4.807, 0),
                (0, 0, 23.571),
                3.588,
                4.807,
                23.571,
                406.5412,
            ),
            ((4, 3, 0), (0, 1, 0), (0, 0, 1), 5, 1, 1, 4),
            ((1, 0, 0), (4, 3, 0), (0, 0, 1), 1, 5, 1, 3),
            ((1, 0, 0), (0, 1, 0), (0, 4, 3), 1, 1, 5, 3),
        ],
    )
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
        model.cell = [[1, 2, 3], [4, 5, 6], [7, 8, 8]]
        assert (model.b1 == 2 * pi / 3 * np.cross(b, c)).all()
        assert (model.b2 == 2 * pi / 3 * np.cross(c, a)).all()
        assert (model.b3 == 2 * pi / 3 * np.cross(a, b)).all()

    def test_cell_list(self):
        model = ExchangeModel()
        model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        Cr1 = Atom("Cr1", (1, 6, 2))
        Cr2 = Atom("Cr2", (1, 3, 2))
        Cr3 = Atom("Cr3", (1, 6, 2))
        model.add_bond(ExchangeParameter(iso=1), Cr1, Cr2, (0, 0, 0))
        model.add_bond(ExchangeParameter(iso=2), Cr1, Cr3, (0, -1, 0))
        model.add_bond(ExchangeParameter(iso=3), Cr2, Cr1, (0, 0, -3))
        cells = model.cell_list.tolist()
        cells = [tuple(i) for i in cells]
        assert len(cells) == 3
        assert (0, 0, 0) in cells
        assert (0, -1, 0) in cells
        assert (0, 0, -3) in cells

    def test_number_spins_in_unit_cell(self):
        Cr1 = Atom("Cr1", (1, 6, 2))
        Cr2 = Atom("Cr2", (1, 3, 5))
        Cr3 = Atom("Cr3", (1, 3, 3))
        model = ExchangeModel()
        model.add_atom(Cr1)
        assert model.number_spins_in_unit_cell == 1
        model.add_atom(Cr2)
        assert model.number_spins_in_unit_cell == 2
        model.add_atom(Cr3)
        assert model.number_spins_in_unit_cell == 3
        model.remove_atom(Cr3)
        Cr3 = Atom("Cr3", (4, 9, 1))
        model.add_atom(Cr3)
        assert model.number_spins_in_unit_cell == 3
        model.remove_atom(Cr1)
        assert model.number_spins_in_unit_cell == 2
        model.remove_atom(Cr2)
        assert model.number_spins_in_unit_cell == 1
        model.remove_atom(Cr3)
        assert model.number_spins_in_unit_cell == 0

    def test_space_dimensions(self):
        model = ExchangeModel()
        model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        Cr1 = Atom("Cr1", (1, 6, 2))
        Cr2 = Atom("Cr2", (1, 3, 5))
        Cr3 = Atom("Cr3", (1, 3, 3))
        model.add_bond(ExchangeParameter(iso=1), Cr1, Cr2, (0, 0, 0))
        model.add_bond(ExchangeParameter(iso=2), Cr1, Cr3, (0, -1, 0))
        model.add_bond(ExchangeParameter(iso=3), Cr2, Cr1, (0, 0, -3))
        x_min, y_min, z_min, x_max, y_max, z_max = model.space_dimensions
        assert x_min == 1
        assert y_min == -7
        assert z_min == -28
        assert x_max == 1
        assert y_max == 6
        assert z_max == 5

    def test_add_bond(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (2, 5, 1))
        Cr2 = Atom("Cr2", (4, 2, 1))
        Cr3 = Atom("Cr3", (5, 1, 8))
        bond12 = ExchangeParameter(iso=12)
        bond13 = ExchangeParameter(iso=13)
        bond23 = ExchangeParameter(iso=23)
        bond31 = ExchangeParameter(iso=31)
        model.add_bond(bond12, Cr1, Cr2, (0, 0, 0))
        model.add_bond(bond23, Cr2, Cr3, (0, 0, 0))
        model.add_bond(bond31, Cr3, Cr1, (0, 0, 0))
        assert len(model.bonds) == 3
        assert (Cr1, Cr2, (0, 0, 0)) in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) in model
        model.add_bond(bond13, Cr1, Cr3, (0, 0, 0))
        assert len(model.bonds) == 4
        assert (Cr1, Cr2, (0, 0, 0)) in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) in model
        assert (Cr1, Cr3, (0, 0, 0)) in model

    def test_remove_bond(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (2, 5, 1))
        Cr2 = Atom("Cr2", (4, 2, 1))
        Cr3 = Atom("Cr3", (5, 1, 8))
        bond12 = ExchangeParameter(iso=12)
        bond13 = ExchangeParameter(iso=13)
        bond23 = ExchangeParameter(iso=23)
        bond31 = ExchangeParameter(iso=31)
        model.add_bond(bond12, Cr1, Cr2, (0, 0, 0))
        model.add_bond(bond23, Cr2, Cr3, (0, 0, 0))
        model.add_bond(bond31, Cr3, Cr1, (0, 0, 0))
        model.add_bond(bond13, Cr1, Cr3, (0, 0, 0))
        assert len(model.bonds) == 4
        assert (Cr1, Cr2, (0, 0, 0)) in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) in model
        assert (Cr1, Cr3, (0, 0, 0)) in model
        model.remove_bond(Cr1, Cr2, (0, 0, 0))
        assert len(model.bonds) == 3
        assert (Cr1, Cr2, (0, 0, 0)) not in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) in model
        assert (Cr1, Cr3, (0, 0, 0)) in model
        model.remove_bond(Cr3, Cr1, (0, 0, 0))
        assert len(model.bonds) == 2
        assert (Cr1, Cr2, (0, 0, 0)) not in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) not in model
        assert (Cr1, Cr3, (0, 0, 0)) in model

    def test_add_atom(self):
        model = ExchangeModel()
        assert len(model.magnetic_atoms) == 0
        Cr1 = Atom("Cr1", (2, 5, 1))
        model.add_atom(Cr1)
        assert len(model.magnetic_atoms) == 1
        Cr2 = Atom("Cr2", (4, 2, 1))
        model.add_atom(Cr2)
        assert len(model.magnetic_atoms) == 2
        assert model.magnetic_atoms[0].position[0] == 2
        assert model.magnetic_atoms[0].position[1] == 5
        assert model.magnetic_atoms[0].position[2] == 1
        assert model.magnetic_atoms[1].position[0] == 4
        assert model.magnetic_atoms[1].position[1] == 2
        assert model.magnetic_atoms[1].position[2] == 1
        model.remove_atom(Cr2)
        Cr2 = Atom("Cr2", (4, 3, 1))
        model.add_atom(Cr2)
        assert model.magnetic_atoms[1].position[0] == 4
        assert model.magnetic_atoms[1].position[1] == 3
        assert model.magnetic_atoms[1].position[2] == 1

    def test_remove_atom(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (2, 5, 1))
        Cr2 = Atom("Cr2", (4, 2, 1))
        Cr3 = Atom("Cr3", (5, 1, 8))
        bond12 = ExchangeParameter(iso=12)
        bond13 = ExchangeParameter(iso=13)
        bond23 = ExchangeParameter(iso=23)
        bond31 = ExchangeParameter(iso=31)
        model.add_bond(bond12, Cr1, Cr2, (0, 0, 0))
        model.add_bond(bond23, Cr2, Cr3, (0, 0, 0))
        model.add_bond(bond31, Cr3, Cr1, (0, 0, 0))
        model.add_bond(bond13, Cr1, Cr3, (0, 0, 0))
        assert len(model.bonds) == 4
        assert (Cr1, Cr2, (0, 0, 0)) in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) in model
        assert (Cr1, Cr3, (0, 0, 0)) in model
        assert len(model.magnetic_atoms) == 3
        model.remove_atom(Cr1)
        assert len(model.magnetic_atoms) == 2
        assert len(model.bonds) == 1
        assert (Cr1, Cr2, (0, 0, 0)) not in model
        assert (Cr2, Cr3, (0, 0, 0)) in model
        assert (Cr3, Cr1, (0, 0, 0)) not in model
        assert (Cr1, Cr3, (0, 0, 0)) not in model

    def test_get_atom_coordinates(self):
        model = ExchangeModel()
        model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        Cr1 = Atom("Cr1", (1, 6, 2))
        model.add_atom(Cr1)
        x, y, z = model.get_atom_coordinates(Cr1)
        assert x == 1
        assert y == 6
        assert z == 2
        x, y, z = model.get_atom_coordinates(Cr1, R=[1, 0, 0])
        assert x == 11
        assert y == 6
        assert z == 2
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, -1, 0])
        assert x == 1
        assert y == -4
        assert z == 2
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, 0, 2])
        assert x == 1
        assert y == 6
        assert z == 22
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, -3, 2])
        assert x == 1
        assert y == -24
        assert z == 22
        x, y, z = model.get_atom_coordinates(Cr1, R=[3, -3, 2])
        assert x == 31
        assert y == -24
        assert z == 22
        model.cell = [[10, 10, 0], [0, 10, 10], [0, 0, 10]]
        x, y, z = model.get_atom_coordinates(Cr1)
        assert x == 1
        assert y == 7
        assert z == 8
        x, y, z = model.get_atom_coordinates(Cr1, R=[1, 0, 0])
        assert x == 11
        assert y == 17
        assert z == 8
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, -1, 0])
        assert x == 1
        assert y == -3
        assert z == -2
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, 0, 2])
        assert x == 1
        assert y == 7
        assert z == 28
        x, y, z = model.get_atom_coordinates(Cr1, R=[0, -3, 2])
        assert x == 1
        assert y == -23
        assert z == -2
        x, y, z = model.get_atom_coordinates(Cr1, R=[3, -3, 2])
        assert x == 31
        assert y == 7
        assert z == -2

    def test_get_bond_vector(self):
        model = ExchangeModel()
        model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        Cr1 = Atom("Cr1", (1, 6, 2))
        Cr2 = Atom("Cr2", (1, 3, 4))
        model.add_atom(Cr1)
        model.add_atom(Cr2)
        vector = model.get_bond_vector(Cr1, Cr2)
        assert (vector == np.array([0, -3, 2])).all()
        vector = model.get_bond_vector(Cr1, Cr2, R=[3, 2, -5])
        assert (vector == np.array([30, 17, -48])).all()

    def test_get_distance(self):
        model = ExchangeModel()
        model.cell = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        Cr1 = Atom("Cr1", (1, 6, 2))
        Cr2 = Atom("Cr2", (1, 3, 4))
        model.add_atom(Cr1)
        model.add_atom(Cr2)
        d = model.get_distance(Cr1, Cr2)
        assert d == sqrt(13)
        d = model.get_distance(Cr1, Cr2, R=[3, 2, -5])
        assert d == sqrt(900 + 17**2 + 48**2)

    def test_filter(self):
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (0.25, 0.25, 0))
        Cr2 = Atom("Cr2", (0.75, 0.75, 0))
        bonds = [
            (12, Cr1, Cr2, (0, 0, 0)),
            (12, Cr2, Cr1, (0, 0, 0)),
            (12, Cr1, Cr1, (1, 0, 0)),
            (12, Cr1, Cr1, (-1, 0, 0)),
            (12, Cr2, Cr2, (1, 0, 0)),
            (12, Cr2, Cr2, (-1, 0, 0)),
            (12, Cr1, Cr1, (0, 2, 0)),
            (12, Cr1, Cr1, (0, -2, 0)),
            (12, Cr2, Cr2, (0, 2, 0)),
            (12, Cr2, Cr2, (0, -2, 0)),
            (12, Cr2, Cr1, (2, 2, 0)),
            (12, Cr1, Cr2, (-2, -2, 0)),
        ]
        for iso, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(iso=iso), atom1, atom2, R)

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
        # In the template there are literals, not objects, because in general
        # template only has information about literals and R.
        filtered_model = model.filtered(template=[("Cr1", "Cr2", (0, 0, 0))])
        assert len(filtered_model.bonds) == 1
        filtered_model = model.filtered(
            template=[("Cr1", "Cr2", (0, 0, 0))], R_vector=[(0, 0, 0), (1, 0, 0)]
        )
        assert len(filtered_model.bonds) == 1

    def test_force_symmetry(self):
        template1 = ExchangeTemplate()
        template2 = ExchangeTemplate()
        template3 = ExchangeTemplate()
        model = ExchangeModel()
        Cr1 = Atom("Cr1", (0.25, 0.25, 0))
        Cr2 = Atom("Cr2", (0.75, 0.75, 0))
        bonds = [
            ([[4, 2, 3], [4, 8, 6], [7, 8, 12]], Cr1, Cr2, (0, 0, 0)),
            ([[1, 4, 7], [2, 5, 8], [3, 6, 9]], Cr2, Cr1, (0, 0, 0)),
            ([[1, 4, 0], [4, 1, -1], [0, 1, 0]], Cr1, Cr1, (1, 0, 0)),
            ([[1, 4, 0], [4, 1, 1], [0, -1, 0]], Cr1, Cr1, (-1, 0, 0)),
            ([[1, 2, 0], [2, 1, -2], [0, 2, 0]], Cr2, Cr2, (1, 0, 0)),
            ([[1, 2, 0], [2, 1, 2], [0, -2, 0]], Cr2, Cr2, (-1, 0, 0)),
        ]
        for matrix, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(matrix=matrix), atom1, atom2, R)

        template1.names = {
            "J1": [("Cr1", "Cr2", (0, 0, 0)), ("Cr2", "Cr1", (0, 0, 0))],
            "J2": [
                ("Cr1", "Cr1", (1, 0, 0)),
                ("Cr1", "Cr1", (-1, 0, 0)),
                ("Cr2", "Cr2", (1, 0, 0)),
                ("Cr2", "Cr2", (-1, 0, 0)),
            ],
        }
        template2.names = {
            "J1": [("Cr1", "Cr2", (0, 0, 0)), ("Cr2", "Cr1", (0, 0, 0))],
            "J2": [
                ("Cr1", "Cr1", (1, 0, 0)),
                ("Cr1", "Cr1", (-1, 0, 0)),
                ("Cr2", "Cr2", (1, 0, 0)),
            ],
        }
        template3.names = {
            "J1": [("Cr1", "Cr2", (0, 0, 0)), ("Cr2", "Cr1", (0, 0, 0))],
            "J2": [("Cr1", "Cr1", (1, 0, 0)), ("Cr1", "Cr1", (-1, 0, 0))],
        }

        model.force_symmetry(template=template1)
        assert len(model.bonds) == 6
        assert (
            model[(Cr1, Cr2, (0, 0, 0))].matrix
            == np.array([[2.5, 2, 3], [4, 6.5, 6], [7, 8, 10.5]])
        ).all()
        assert (
            model[(Cr2, Cr1, (0, 0, 0))].matrix
            == np.array([[2.5, 4, 7], [2, 6.5, 8], [3, 6, 10.5]])
        ).all()

        assert (
            model[(Cr1, Cr1, (1, 0, 0))].matrix
            == np.array([[1, 3, 0], [3, 1, -1.5], [0, 1.5, 0]])
        ).all()
        assert (
            model[(Cr1, Cr1, (-1, 0, 0))].matrix
            == np.array([[1, 3, 0], [3, 1, 1.5], [0, -1.5, 0]])
        ).all()
        assert (
            model[(Cr2, Cr2, (1, 0, 0))].matrix
            == np.array([[1, 3, 0], [3, 1, -1.5], [0, 1.5, 0]])
        ).all()
        assert (
            model[(Cr2, Cr2, (-1, 0, 0))].matrix
            == np.array([[1, 3, 0], [3, 1, 1.5], [0, -1.5, 0]])
        ).all()
        for matrix, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(matrix=matrix), atom1, atom2, R)
        model.force_symmetry(template=template2)
        assert len(model.bonds) == 5
        assert (
            model[(Cr1, Cr2, (0, 0, 0))].matrix
            == np.array([[2.5, 2, 3], [4, 6.5, 6], [7, 8, 10.5]])
        ).all()
        assert (
            model[(Cr2, Cr1, (0, 0, 0))].matrix
            == np.array([[2.5, 4, 7], [2, 6.5, 8], [3, 6, 10.5]])
        ).all()

        assert (
            model[(Cr1, Cr1, (1, 0, 0))].matrix
            == np.array([[1, 10 / 3, 0], [10 / 3, 1, -4 / 3], [0, 4 / 3, 0]])
        ).all()
        assert (
            model[(Cr1, Cr1, (-1, 0, 0))].matrix
            == np.array([[1, 10 / 3, 0], [10 / 3, 1, 4 / 3], [0, -4 / 3, 0]])
        ).all()
        assert (
            model[(Cr2, Cr2, (1, 0, 0))].matrix
            == np.array([[1, 10 / 3, 0], [10 / 3, 1, -4 / 3], [0, 4 / 3, 0]])
        ).all()
        assert (Cr2, Cr2, (-1, 0, 0)) not in model

        for matrix, atom1, atom2, R in bonds:
            model.add_bond(ExchangeParameter(matrix=matrix), atom1, atom2, R)
        model.force_symmetry(template=template3)
        assert len(model.bonds) == 4
        assert (
            model[(Cr1, Cr2, (0, 0, 0))].matrix
            == np.array([[2.5, 2, 3], [4, 6.5, 6], [7, 8, 10.5]])
        ).all()
        assert (
            model[(Cr2, Cr1, (0, 0, 0))].matrix
            == np.array([[2.5, 4, 7], [2, 6.5, 8], [3, 6, 10.5]])
        ).all()

        assert (
            model[(Cr1, Cr1, (1, 0, 0))].matrix
            == np.array([[1, 4, 0], [4, 1, -1], [0, 1, 0]])
        ).all()
        assert (
            model[(Cr1, Cr1, (-1, 0, 0))].matrix
            == np.array([[1, 4, 0], [4, 1, 1], [0, -1, 0]])
        ).all()

        assert (Cr2, Cr2, (11, 0, 0)) not in model
        assert (Cr2, Cr2, (-1, 0, 0)) not in model

    def test_notation_manipulation(self):
        model = ExchangeModel()
        Cr = Atom("Cr", (0, 0, 0), spin=3 / 2)
        model[Cr, Cr, (1, 0, 0)] = ExchangeParameter(iso=1)
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        with pytest.raises(NotationError):
            model.double_counting
        with pytest.raises(NotationError):
            model.spin_normalized
        with pytest.raises(NotationError):
            model.factor_one_half
        with pytest.raises(NotationError):
            model.factor_two
        with pytest.raises(NotationError):
            model.minus_sign

        model.notation = "standard"
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        model.notation = "standard"
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        assert model.double_counting
        model.double_counting = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 2
        model.double_counting = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        assert not model.spin_normalized
        model.spin_normalized = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 9 / 4
        model.spin_normalized = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        assert not model.factor_one_half
        assert not model.factor_two
        model.factor_one_half = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 2
        model.factor_one_half = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        model.factor_two = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 0.5
        model.factor_two = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        model.factor_one_half = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 2
        model.factor_two = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        model.factor_two = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        model.factor_one_half = False
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        assert model.minus_sign
        model.minus_sign = False
        assert model[Cr, Cr, (1, 0, 0)].iso == -1
        model.minus_sign = True
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        model.notation = "SpinW"
        assert model[Cr, Cr, (1, 0, 0)].iso == -1
        model.notation = "TB2J"
        assert model[Cr, Cr, (1, 0, 0)].iso == 9 / 4

    def test_predefined_notations(self):
        model = ExchangeModel()
        Cr = Atom("Cr", (0, 0, 0), spin=3 / 2)
        model[Cr, Cr, (1, 0, 0)] = ExchangeParameter(iso=1)
        assert model[Cr, Cr, (1, 0, 0)].iso == 1

        with pytest.raises(NotationError):
            model.double_counting
        with pytest.raises(NotationError):
            model.spin_normalized
        with pytest.raises(NotationError):
            model.factor_one_half
        with pytest.raises(NotationError):
            model.factor_two
        with pytest.raises(NotationError):
            model.minus_sign

        model.notation = "standard"
        assert model[Cr, Cr, (1, 0, 0)].iso == 1
        assert (
            model.double_counting
            and not model.spin_normalized
            and not model.factor_one_half
            and not model.factor_two
            and model.minus_sign
        )

        model.notation = "TB2J"
        assert model[Cr, Cr, (1, 0, 0)].iso == 9 / 4
        assert (
            model.double_counting
            and model.spin_normalized
            and not model.factor_one_half
            and not model.factor_two
            and model.minus_sign
        )

        model.notation = "SpinW"
        assert model[Cr, Cr, (1, 0, 0)].iso == -1
        assert (
            model.double_counting
            and not model.spin_normalized
            and not model.factor_one_half
            and not model.factor_two
            and not model.minus_sign
        )

    def test_ferromagnetic_energy(self):
        model = ExchangeModel()
        model.cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        Cr1 = Atom("Cr1", (1, 1, 1))
        Cr2 = Atom("Cr2", (1, 1, 1))
        Cr3 = Atom("Cr3", (1, 1, 1))
        model.add_bond(ExchangeParameter(iso=1), Cr1, Cr2, (0, 0, 0))
        model.add_bond(ExchangeParameter(iso=2), Cr1, Cr3, (0, -1, 0))
        model.add_bond(ExchangeParameter(iso=3), Cr2, Cr1, (0, 0, -3))
        assert model.ferromagnetic_energy() == -6
        assert model.ferromagnetic_energy(theta=23, phi=234) == -6
        model.add_bond(
            ExchangeParameter(iso=3, aniso=[[1, 0, 0], [0, 2, 0], [0, 0, -3]]),
            Cr2,
            Cr1,
            (0, 0, -1),
        )
        assert model.ferromagnetic_energy() == -6
        assert model.ferromagnetic_energy(theta=90) == -10
        assert model.ferromagnetic_energy(theta=90, phi=90) == -11
        assert model.ferromagnetic_energy(theta=90, phi=45) == -10.5
