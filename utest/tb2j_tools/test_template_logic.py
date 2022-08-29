import os

import pytest

from rad_tools.tb2j_tools.template_logic import ExchangeTemplate


class TestInputFilename():

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            tmp_model = ExchangeTemplate(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            tmp_model = ExchangeTemplate(
                "Ah, music. A magic beyond all we do here!")

    def test_correct_filename(self):
        tmp_model = ExchangeTemplate(os.path.join(
            'utest', 'tb2j_tools', 'resourses', 'exchange.out'))


class TestExchangeTemplate:

    tmp_model = ExchangeTemplate(os.path.join(
        'utest', 'tb2j_tools', 'resourses', 'template.txt'
    ))


class TestReadFunctions(TestExchangeTemplate):

    def count_entries(self, dictionary):
        i = 0
        for neighbor in dictionary:
            for atom_1 in dictionary[neighbor]:
                for atom_2 in dictionary[neighbor][atom_1]:
                    for R in dictionary[neighbor][atom_1][atom_2]:
                        i += 1
        return i

    def test_read_neighbors(self):
        assert self.count_entries(self.tmp_model.template) == 28
        assert self.tmp_model.names_for_plot == {'J1': '$J_1$',
                                                 'J1(1)': '$J_1(1)$',
                                                 'J1(2)': '$J_1(2)$',
                                                 'J2': 'J2',
                                                 'J2(1)': '$J_2(1)$',
                                                 'J2(2)': 'J2(2)',
                                                 'J3': '$J_3$'}

    def test_template(self):
        assert self.tmp_model.template == {'J1': {
            'Cr1': {'Cr1': {(1, 0, 0): True, (-1, 0, 0): True}},
            'Cr2': {'Cr2': {(1, 0, 0): True, (-1, 0, 0): True}}
        },

            'J1(1)': {
            'Cr1': {'Cr1': {(1, 0, 0): True}},
            'Cr2': {'Cr2': {(-1, 0, 0): True}}
        },

            'J1(2)': {
            'Cr1': {'Cr1': {(-1, 0, 0): True}},
            'Cr2': {'Cr2': {(1, 0, 0): True}}
        },
            'J2': {
            'Cr1': {'Cr2': {(0, 0, 0): True,
                            (1, 0, 0): True,
                            (0, -1, 0): True,
                            (1, -1, 0): True}},
            'Cr2': {'Cr1': {(0, 0, 0): True,
                            (-1, 0, 0): True,
                            (0, 1, 0): True,
                            (-1, 1, 0): True}}
        },
            'J2(1)': {
            'Cr1': {'Cr2': {(1, 0, 0): True,
                            (1, -1, 0): True}},
            'Cr2': {'Cr1': {(-1, 0, 0): True,
                            (-1, 1, 0): True}}
        },
            'J2(2)': {
            'Cr1': {'Cr2': {(0, 0, 0): True,
                            (0, -1, 0): True}},
            'Cr2': {'Cr1': {(0, 0, 0): True,
                            (0, 1, 0): True}}
        },
            'J3': {
            'Cr1': {'Cr1': {(0, 1, 0): True,
                            (0, -1, 0): True}},
            'Cr2': {'Cr2': {(0, 1, 0): True,
                            (0, -1, 0): True}}
        }
        }
        assert len(self.tmp_model.plained_template) == 28
