from copy import deepcopy

import numpy as np


class ExchangeModel:

    def __init__(self, filename) -> None:
        self.headers_to_functions = {
            'Cell (Angstrom):\n': self.read_cell,
            'Atoms:  \n': self.read_atoms,
            'Orbitals used in decomposition: \n': self.read_orb_decomposition,
            'Exchange: \n': self.read_exchange
        }
        self.major_sep = '=' * 90 + '\n'
        self.minor_sep = '-' * 88 + '\n'
        self.cell = None
        self.atoms = None
        self.magnetic_atoms = []
        self.iso = None
        self.aniso = None
        self.dmi = None
        self.distance = None
        self.resolved_data_attributes = {
            'iso': self.iso,
            'aniso': self.aniso,
            'dmi': self.dmi,
            'distance': self.distance
        }
        self.resolved_data = [
            None,
            None,
            None,
            None
        ]
        self.resolved_data_nofilter = [
            None,
            None,
            None,
            None
        ]
        self.resolved_data_keys = {
            'iso': 0,
            'aniso': 1,
            'dmi': 2,
            'distance': 3
        }
        self.resolved_data_names_to_read = {
            'iso': 'J_iso:',
            'aniso': 'J_ani:',
            'dmi': 'DMI:',
            'distance': None
        }
        self.resolved_data_format_functions = {
            'iso': self.format_iso,
            'aniso': self.format_aniso,
            'dmi': self.format_dmi,
            'distance': None
        }
        with open(filename, 'r') as file:
            self.data = file.readlines()
        i = 0
        while i < len(self.data):
            if self.data[i] in self.headers_to_functions:
                i = self.headers_to_functions[self.data[i]](i)
            else:
                i += 1

    def update_attributes(self):
        for

    def read_cell(self, i):
        i += 1
        self.cell = np.array(
            [
                list(map(float, self.data[i].split())),
                list(map(float, self.data[i + 1].split())),
                list(map(float, self.data[i + 2].split()))
            ],
            dtype=float
        )
        return i + 3

    def read_atoms(self, i):
        self.atoms = {}
        i += 3
        while 'Total' not in self.data[i]:
            self.atoms[self.data[i].split()[0]] = list(
                map(float, self.data[i].split()[1:]))
            i += 1
        return i + 1

    def read_orb_decomposition(self, i):
        self.orb_for_decomposition = {}
        i += 2
        while ':' in self.data[i]:
            self.orb_for_decomposition[self.data[i].split()[0]] = \
                self.data[i].replace('[', '').replace(']', '')\
                .replace(',', '').replace('\'', '').split()[2:]
            i += 1
        return i

    def prepare_dicts(self, atom_1, atom_2):
        for d in range(0, len(self.resolved_data)):
            if self.resolved_data[d] is None:
                self.resolved_data[d] = {}
            if atom_1 not in self.resolved_data[d]:
                self.resolved_data[d][atom_1] = {}
            if atom_2 not in self.resolved_data[d][atom_1]:
                self.resolved_data[d][atom_1][atom_2] = {}

    def format_iso(self, i):
        return float(self.data[i].split()[1])

    def format_aniso(self, i):
        tmp = []
        for j in range(0, 3):
            i += 1
            tmp.append(list(map(float, self.data[i]
                                .replace('[', '').replace(']', '').split())))
        return tmp

    def format_dmi(self, i):
        return tuple(map(float, self.data[i].replace('(', '').replace(')', '')
                         .replace('[Testing!]', '').split()[1:]))

    def read_exchange(self, i):
        i += 2
        while i < len(self.data):
            while i < len(self.data) and self.minor_sep not in self.data[i]:
                for name in self.resolved_data_names_to_read:
                    if self.resolved_data_names_to_read[name]\
                        and self.resolved_data_names_to_read[name]\
                            in self.data[i]:
                        self.resolved_data[
                            self.resolved_data_keys[name]
                        ][atom_1][atom_2][R] =\
                            self.resolved_data_format_functions[name](i)
                i += 1
            else:
                if i < len(self.data):
                    i += 1
                    tmp = self.data[i].replace(
                        '(', '').replace(')', '').replace(',', '').split()
                    atom_1 = tmp[0]
                    atom_2 = tmp[1]
                    R = tuple(map(int, tmp[2:5]))
                    distance = float(tmp[9])
                    if atom_1 not in self.magnetic_atoms:
                        self.magnetic_atoms.append(atom_1)
                    if atom_2 not in self.magnetic_atoms:
                        self.magnetic_atoms.append(atom_2)
                    self.prepare_dicts(atom_1, atom_2)
                    self.resolved_data[
                        self.resolved_data_keys['distance']
                    ][atom_1][atom_2][R] = distance
            i += 1
        return i

    def filter(self, distance=None, template=None):
        if distance is not None and template is not None:
            raise ValueError("Do not try to filter by distance and\
                             by template at the same time, please")
        if distance is not None:
            for d in range(0, len(self.resolved_data)):
                if self.resolved_data_nofilter[d] is None:
                    self.resolved_data_nofilter[d] =\
                        deepcopy(self.resolved_data[d])
                else:
                    self.resolved_data[d] =\
                        deepcopy(self.resolved_data_nofilter[d])

            for atom_1 in self.resolved_data_nofilter[d]:
                for atom_2 in self.resolved_data_nofilter[d][atom_1]:
                    for R in self.resolved_data_nofilter[d][atom_1][atom_2]:
                        del self.resolved_data[d][atom_1][atom_2][R]
                    if not self.resolved_data[d][atom_1][atom_2]:
                        del self.resolved_data[d][atom_1][atom_2]
                if not self.resolved_data[d][atom_1]:
                    del self.resolved_data[d][atom_1]
