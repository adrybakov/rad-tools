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
            'tb2j_tools', 'resourses', 'exchange.out'))


class TestExchangeTemplate:

    tmp_model = ExchangeTemplate(os.path.join(
        'tb2j_tools', 'resourses', 'template.txt'
    ))


class TestReadFunctions(TestExchangeTemplate):

    def count_entries(self, dictionary):
        i = 0
        for name in dictionary:
            i += len(dictionary[name])
        return i

    def test_read_neighbors(self):
        assert self.count_entries(self.tmp_model.names) == 28
        assert self.tmp_model.latex_names == {'J1': '$J_1$',
                                              'J1(1)': '$J_1(1)$',
                                              'J1(2)': '$J_1(2)$',
                                              'J2': 'J2',
                                              'J2(1)': '$J_2(1)$',
                                              'J2(2)': '$J_2(2)$',
                                              'J3': '$J_3$'}

    def test_template(self):
        assert self.tmp_model.names == {'J1': [
            ('Cr1', 'Cr1', (1, 0, 0)),
            ('Cr1', 'Cr1', (-1, 0, 0)),
            ('Cr2', 'Cr2', (1, 0, 0)),
            ('Cr2', 'Cr2', (-1, 0, 0))
        ],

            'J1(1)': [
            ('Cr1', 'Cr1', (1, 0, 0)),
            ('Cr2', 'Cr2', (-1, 0, 0))
        ],

            'J1(2)': [
            ('Cr1', 'Cr1', (-1, 0, 0)),
            ('Cr2', 'Cr2', (1, 0, 0))
        ],
            'J2': [
            ('Cr1', 'Cr2', (0, 0, 0)),
            ('Cr1', 'Cr2', (1, 0, 0)),
            ('Cr1', 'Cr2', (0, -1, 0)),
            ('Cr1', 'Cr2', (1, -1, 0)),
            ('Cr2', 'Cr1', (0, 0, 0)),
            ('Cr2', 'Cr1', (-1, 0, 0)),
            ('Cr2', 'Cr1', (0, 1, 0)),
            ('Cr2', 'Cr1', (-1, 1, 0))
        ],
            'J2(1)': [
            ('Cr1', 'Cr2', (1, 0, 0)),
            ('Cr1', 'Cr2', (1, -1, 0)),
            ('Cr2', 'Cr1', (-1, 0, 0)),
            ('Cr2', 'Cr1', (-1, 1, 0))
        ],
            'J2(2)': [
            ('Cr1', 'Cr2', (0, 0, 0)),
            ('Cr1', 'Cr2', (0, -1, 0)),
            ('Cr2', 'Cr1', (0, 0, 0)),
            ('Cr2', 'Cr1', (0, 1, 0))
        ],
            'J3': [
            ('Cr1', 'Cr1', (0, 1, 0)),
            ('Cr1', 'Cr1', (0, -1, 0)),
            ('Cr2', 'Cr2', (0, 1, 0)),
            ('Cr2', 'Cr2', (0, -1, 0))
        ]
        }
        assert len(self.tmp_model.get_list()) == 28
