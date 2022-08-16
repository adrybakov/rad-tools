import os

import pytest
import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModel


class TestExchangeModel:

    tmp_model = ExchangeModel(os.path.join(
        'rad_tools', 'utest', 'tb2j_tools', 'resourses', 'exchange.out'
    ))


class TestInputFilename(TestExchangeModel):

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            tmp_model = ExchangeModel(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            tmp_model = ExchangeModel(
                "Ah, music. A magic beyond all we do here!")

    def test_correct_filename(self):
        tmp_model = ExchangeModel(os.path.join(
            'rad_tools', 'utest', 'tb2j_tools', 'resourses', 'exchange.out'))


class TestReadFunctions(TestExchangeModel):

    def test_read_cell(self):
        assert self.tmp_model.cell is not None
        cell_values = [[3.588, 0.000, 0.000],
                       [0.000,  4.807,  0.000],
                       [0.000,  0.000, 23.571]]
        for i in range(0, 3):
            for j in range(0, 3):
                assert self.tmp_model.cell[i][j] == cell_values[i][j]

    def test_read_atoms(self):
        assert self.tmp_model.atoms is not None
        atoms_value = {
            'Br1': [0.8970, 1.2018, -0.0668, 5.6114,  0.0000, -0.0000, -0.0535],
            'Cr1': [2.6910, 1.2018,  1.7371, 4.3510, -0.0000,  0.0000,  3.2857],
            'S1':  [2.6910, 3.6054,  2.2030, 5.0369, -0.0000,  0.0000, -0.2315],
            'Br2': [2.6910, 3.6054,  5.6376, 5.6111,  0.0000,  0.0000, -0.0537],
            'Cr2': [0.8970, 3.6054,  3.8336, 4.3519, -0.0000,  0.0000,  3.2863],
            'S2':  [0.8970, 1.2018,  3.3678, 5.0369, -0.0000, -0.0000, -0.2316]
        }
        for key in self.tmp_model.atoms:
            assert key in atoms_value
        for key in atoms_value:
            assert key in self.tmp_model.atoms
            assert atoms_value[key] == self.tmp_model.atoms[key]

    def test_read_orb_decomposition(self):
        assert self.tmp_model.orb_for_decomposition is not None
        orb_decomposition_value = {
            'Cr1': ['orb_1', 'orb_2', 'orb_3', 'orb_4', 'orb_5',
                    'orb_6', 'orb_7', 'orb_8', 'orb_9', 'orb_10'],
            'Cr2': ['orb_1', 'orb_2', 'orb_3', 'orb_4', 'orb_5',
                    'orb_6', 'orb_7', 'orb_8', 'orb_9', 'orb_10']
        }
        for key in self.tmp_model.orb_for_decomposition:
            assert key in orb_decomposition_value
        for key in orb_decomposition_value:
            assert key in self.tmp_model.orb_for_decomposition
            assert orb_decomposition_value[key] ==\
                self.tmp_model.orb_for_decomposition[key]

    def test_read_exchange_basic(self):
        assert self.tmp_model.iso is not None
        assert self.tmp_model.aniso is not None
        assert self.tmp_model.dmi is not None
        assert self.tmp_model.distance is not None

    @pytest.mark.parametrize("atom_1,atom_2,R,iso,aniso,dmi,distance",
                             [
                                 ('Cr1', 'Cr1',
                                  (-1, 0, 0),
                                  3.5386,
                                  [[-0.032, 0, 0],
                                   [0, -0.054, 0],
                                   [0, 0, -0.028]],
                                  (0, -0.0163, 0),
                                  3.588),

                                 ('Cr1', 'Cr2',
                                  (0, 0, 0),
                                  3.0830,
                                  [[-0.009, 0, 0.004],
                                   [0, -0.013, 0],
                                   [0.004, 0, -0.007]],
                                  (-0.0002, 0.0001, 0.0001),
                                  3.659),

                                 ('Cr1', 'Cr1',
                                  (0, 1, 0),
                                  4.1479,
                                  [[-0.016, 0, 0],
                                   [0, -0.003, 0],
                                   [0, 0, -0.008]],
                                  (-0.0584, 0, 0),
                                  4.807),

                                 ('Cr1', 'Cr1',
                                  (-1, 1, 0),
                                  0.0422,
                                  [[-0.007, 0, 0],
                                   [0, -0.006, 0],
                                   [0, 0, -0.005]],
                                  (0.0194, -0.0170, 0),
                                  5.999),

                                 ('Cr2', 'Cr2',
                                  (0, -1, 0),
                                  4.1423,
                                  [[-0.016, 0, 0],
                                   [0, -0.003, 0],
                                   [0, 0, -0.008]],
                                  (-0.0568, 0, 0),
                                  4.807),

                                 ('Cr2', 'Cr2',
                                  (0, 2, 0),
                                  0.1209,
                                  [[-0.001, 0, 0],
                                   [0, 0, 0],
                                   [0, 0, 0]],
                                  (0.0363, 0, 0),
                                  9.614),

                                 ('Cr2', 'Cr1',
                                  (1, -1, 0),
                                  0.0038,
                                  [[-0.001, 0, 0],
                                   [0, -0.001, 0],
                                   [0, 0, -0.001]],
                                  (0, 0, 0),
                                  9.239),

                                 ('Cr1', 'Cr2',
                                  (1, -2, 0),
                                  0.5503,
                                  [[-0.001, 0, 0],
                                   [0, -0.002, 0],
                                   [0, 0, -0.001]],
                                  (0.0001, -0.0001, 0),
                                  7.721),

                                 ('Cr2', 'Cr1',
                                  (-1, 1, 0),
                                  3.0830,
                                  [[-0.009, 0, -0.004],
                                   [0, -0.013, 0],
                                   [-0.004, 0, -0.007]],
                                  (-0.0002, 0.0001, -0.0001),
                                  3.659),
                             ])
    def test_read_exchange_examples(self, atom_1, atom_2, R, iso, aniso, dmi, distance):
        assert self.tmp_model.iso[atom_1][atom_2][R] == iso
        assert self.tmp_model.dmi[atom_1][atom_2][R] == dmi
        assert self.tmp_model.distance[atom_1][atom_2][R] == distance
        for i in range(0, 3):
            for j in range(0, 3):
                assert self.tmp_model.aniso[atom_1][atom_2][R][i][j] == aniso[i][j]

    def test_magnetic_atoms(self):
        assert self.tmp_model.magnetic_atoms == ['Cr1', 'Cr2']


class TestFilter(TestExchangeModel):

    def test_wrong_entry(self):
        with pytest.raises(ValueError):
            self.tmp_model.filter(distance=4, template='some_template')

    def count_entries(self, dictionary):
        i = 0
        for atom_1 in dictionary:
            for atom_2 in dictionary[atom_1]:
                for R in dictionary[atom_1][atom_2]:
                    i += 1
        return i

    @pytest.mark.parametrize("distance, elements_number", (
        (4.807, 16),
        (6, 24),
        (0, 0),
        (5, 16),
        (4, 12)
    ))
    def test_filter_by_distance(self, distance, elements_number):
        self.tmp_model.filter(distance=distance)
        assert self.count_entries(self.tmp_model.iso) == elements_number
        assert self.count_entries(self.tmp_model.aniso) == elements_number
        assert self.count_entries(self.tmp_model.dmi) == elements_number
        assert self.count_entries(self.tmp_model.distance) == elements_number
