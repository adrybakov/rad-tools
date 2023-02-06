import os

import pytest

from rad_tools.io.internal import *


class ReadTemplate():

    template = read_template(os.path.join('utest',
                                          'test_io',
                                          'resources',
                                          'exchange.out'))

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            template = read_template(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            template = read_template("Ah, music. " +
                                     "A magic beyond all we do here!")

    def count_entries(self, dictionary):
        i = 0
        for name in dictionary:
            i += len(dictionary[name])
        return i

    def test_read_neighbors(self):
        assert self.count_entries(self.template.names) == 28
        assert self.template.latex_names == {'J1': '$J_1$',
                                             'J1(1)': '$J_1(1)$',
                                             'J1(2)': '$J_1(2)$',
                                             'J2': 'J2',
                                             'J2(1)': '$J_2(1)$',
                                             'J2(2)': '$J_2(2)$',
                                             'J3': '$J_3$'}

    def test_template(self):
        assert self.template.names == {'J1': [
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
        assert len(self.template.get_list()) == 28
