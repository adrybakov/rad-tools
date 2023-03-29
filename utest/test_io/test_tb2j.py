import os

import numpy as np
import pytest

from rad_tools.io import read_template
from rad_tools.io.tb2j import *


class TestReadExchangeModel:

    model = read_exchange_model(os.path.join('utest',
                                             'test_io',
                                             'resources',
                                             'exchange.out'))

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            model = read_exchange_model(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            model = read_exchange_model("Ah, music. " +
                                        "A magic beyond all we do here!")

    def test_read_cell(self):
        cell_values = [[3.588, 0.000, 0.000],
                       [0.000,  4.807,  0.000],
                       [0.000,  0.000, 23.571]]
        assert (self.model.cell == cell_values).all()

    def test_nonmagnetic_atoms(self):
        assert len(self.model.nonmagnetic_atoms) == 4
        assert (np.around(self.model.nonmagnetic_atoms['Br1'] - np.array(
            [0.25, 0.2500104, -0.00283399]), 4) == np.zeros(3)).all()
        assert (np.around(self.model.nonmagnetic_atoms['Br2'] - np.array(
            [0.75, 0.7500312, 0.23917526]), 4) == np.zeros(3)).all()
        assert (np.around(self.model.nonmagnetic_atoms['S1'] - np.array(
            [0.75, 0.7500312, 0.09346231]), 4) == np.zeros(3)).all()
        assert (np.around(self.model.nonmagnetic_atoms['S2'] - np.array(
            [0.25, 0.2500104, 0.14287896]), 4) == np.zeros(3)).all()

    def test_magnetic_atoms(self):
        assert len(self.model.magnetic_atoms) == 2
        assert (np.around(self.model.magnetic_atoms['Cr1'] - np.array(
            [0.75, 0.2500104, 0.07369649]), 4) == np.zeros(3)).all()
        assert (np.around(self.model.magnetic_atoms['Cr2'] - np.array(
            [0.25, 0.7500312, 0.16264053]), 4) == np.zeros(3)).all()

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
        assert round(self.model[(atom1, atom2, R)].iso, 4) == iso
        assert round(self.model.get_distance(
            atom1, atom2, R), 2) == round(distance, 2)
        for i in range(0, 3):
            assert round(
                self.model[(atom1, atom2, R)].dmi[i], 4) == dmi[i]
            for j in range(0, 3):
                assert round(
                    self.model[(atom1, atom2, R)].aniso[i][j], 4) == aniso[i][j]

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

    def test_instance_data_type(self):
        filtered_model = self.model.filtered(max_distance=5)
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
        filtered_model = self.model.filtered(max_distance=max_distance)
        assert len(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("min_distance, elements_number", [
        (0, 72),
        (6, 48),
        (9.6, 4),
        (10, 0)
    ])
    def test_filter_by_min_distance(self, min_distance, elements_number):
        filtered_model = self.model.filtered(min_distance=min_distance)
        assert len(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("R_vector, elements_number", [
        ((0, 0, 0), 2),
        ([(0, 0, 0), (1, 0, 0)], 6),
        ([(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 10),
        ([], 0),
        ([(0, 0, 0), (1, 0, 0), (1, 0, 0)], 6)
    ])
    def test_filter_by_R_vector(self, R_vector, elements_number):
        filtered_model = self.model.filtered(R_vector=R_vector)
        assert len(filtered_model.bonds) == elements_number

    def test_filter_by_template(self):
        filtered_model = self.model.filtered(
            template=self.template)
        assert len(filtered_model.bonds) == 16
        template = read_template(os.path.join('utest',
                                              'test_io',
                                              'resources',
                                              'template.txt'))
        filtered_model = self.model.filtered(
            template=template)
        assert len(filtered_model.bonds) == 16

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
        filtered_model = self.model.filtered(max_distance=max_distance,
                                             min_distance=min_distance)
        assert len(filtered_model.bonds) == elements_number

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
        filtered_model = self.model.filtered(max_distance=max_distance,
                                             template=self.template)
        assert len(filtered_model.bonds) == elements_number

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
        filtered_model = self.model.filtered(min_distance=min_distance,
                                             template=self.template)
        assert len(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("R_vector, elements_number", [
        ((0, 0, 0), 2),
        ([(0, 0, 0), (1, 0, 0)], 5),
        ([(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 8),
        ([], 0),
        ([(0, 0, 0), (1, 0, 0), (1, 0, 0)], 5)
    ])
    def test_filter_by_R_vector_and_template(self, R_vector, elements_number):
        filtered_model = self.model.filtered(
            R_vector=R_vector, template=self.template)
        assert len(filtered_model.bonds) == elements_number

    @ pytest.mark.parametrize("max_distance, R_vector, elements_number", [
        (0, (0, 0, 0), 0),
        (10, (0, 0, 0), 2),
        (5, [(0, 0, 0), (1, 0, 0)], 5),
        (10, [(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 10),
        (5, [(0, 0, 0), (1, 0, 0), (-1, 0, 0)], 8),
        (100, [], 0)
    ])
    def test_filter_by_everything(self, max_distance, R_vector, elements_number):
        filtered_model = self.model.filtered(
            R_vector=R_vector, max_distance=max_distance)
        assert len(filtered_model.bonds) == elements_number
