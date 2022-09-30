import os
from re import template

import pytest
import numpy as np

from rad_tools.tb2j_tools.file_logic import *
from rad_tools.tb2j_tools.template_logic import *


class TestExchangeModelTB2J:

    model = ExchangeModelTB2J(os.path.join(
        'utest', 'tb2j_tools', 'resourses', 'exchange.out'
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
            'utest', 'tb2j_tools', 'resourses', 'exchange.out'))


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

    @pytest.mark.parametrize("atom1, atom2, R, iso, aniso, dmi, distance",
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
        assert self.model.magnetic_atoms['Cr1'] == (2.6910, 1.2018,  1.7371)
        assert self.model.magnetic_atoms['Cr2'] == (0.8970, 3.6054,  3.8336)


class TestFilter(TestExchangeModelTB2J):

    def count_entries(self, dictionary):
        i = 0
        if dictionary is not None:
            for atom1 in dictionary:
                for atom2 in dictionary[atom1]:
                    for R in dictionary[atom1][atom2]:
                        i += 1
        return i

    @pytest.mark.parametrize("max_distance, elements_number", [
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

    @pytest.mark.parametrize("min_distance, elements_number", [
        (0, 72),
        (6, 48),
        (9.6, 4),
        (10, 0)
    ])
    def test_filter_by_min_distance(self, min_distance, elements_number):
        filtered_model = self.model.filter(min_distance=min_distance)
        assert self.count_entries(filtered_model.bonds) == elements_number

    @pytest.mark.parametrize("R_vector, elements_number", [
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

    @pytest.mark.parametrize("min_distance, max_distance, elements_number", [
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

    @pytest.mark.parametrize("max_distance, elements_number", [
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

    @pytest.mark.parametrize("min_distance, elements_number", [
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

    @pytest.mark.parametrize("R_vector, elements_number", [
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

    @pytest.mark.parametrize("max_distance, R_vector, elements_number", [
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


class TestExchangeModel:

    @pytest.mark.parametrize("a, b, c, A, B, C", [
        ((1, 0, 0), (0, 1, 0), (0, 0, 1), 1, 1, 1),
        ((3.588, 0, 0), (0, 4.807, 0), (0, 0, 23.571), 3.588, 4.807, 23.571),
        ((4, 3, 0), (0, 0, 0), (0, 0, 0), 5, 0, 0),
        ((0, 0, 0), (4, 3, 0), (0, 0, 0), 0, 5, 0),
        ((0, 0, 0), (0, 0, 0), (4, 3, 0), 0, 0, 5),
    ])
    def test_get_lattice_vectors_length(self, a, b, c, A, B, C):
        model = ExchangeModel()
        model.cell = np.array([a, b, c], dtype=float)
        assert A == round(model.get_lattice_vectors_length()[0], 4)
        assert B == round(model.get_lattice_vectors_length()[1], 4)
        assert C == round(model.get_lattice_vectors_length()[2], 4)
