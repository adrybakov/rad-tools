import os

import pytest
import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModel


class TestInputFilename:

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


class TestReadFunctions:

    def test_read_cell(self):
        tmp_model = ExchangeModel(os.path.join(
            'rad_tools', 'utest', 'tb2j_tools', 'resourses', 'exchange.out'))
        assert tmp_model.cell is not None
        cell_values = [[3.588, 0.000, 0.000],
                       [0.000,  4.807,  0.000],
                       [0.000,  0.000, 23.571]]
        for i in range(0, 3):
            for j in range(0, 3):
                assert tmp_model.cell[i][j] == cell_values[i][j]
