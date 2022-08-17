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
        self.atom_1 = None
        self.atom_2 = None
        self.R = None
        self.cell = None
        self.atoms = None
        self.magnetic_atoms = []
        self.list_for_sort = []
        self.data = [None, None, None, None]
        self.keys = {'iso': 0,
                     'aniso': 1,
                     'dmi': 2,
                     'distance': 3}
        self.names_to_read = {'iso': 'J_iso:',
                              'aniso': 'J_ani:',
                              'dmi': 'DMI:',
                              'distance': self.minor_sep}
        self.format_functions = {'iso': self.format_iso,
                                 'aniso': self.format_aniso,
                                 'dmi': self.format_dmi,
                                 'distance': self.format_distance}
        with open(filename, 'r') as file:
            self.file_data = file.readlines()
        i = 0
        while i < len(self.file_data):
            if self.file_data[i] in self.headers_to_functions:
                i = self.headers_to_functions[self.file_data[i]](i)
            else:
                i += 1
        self.update_attributes()
        self.data_nofilter = deepcopy(self.data)

    def update_attributes(self):
        self.iso = self.data[self.keys['iso']]
        self.aniso = self.data[self.keys['aniso']]
        self.dmi = self.data[self.keys['dmi']]
        self.distance = self.data[self.keys['distance']]

    def reset_data(self):
        for d in range(0, len(self.data)):
            self.data[d] = None

    def read_cell(self, i):
        i += 1
        self.cell = np.array(
            [
                list(map(float, self.file_data[i].split())),
                list(map(float, self.file_data[i + 1].split())),
                list(map(float, self.file_data[i + 2].split()))
            ],
            dtype=float
        )
        return i + 3

    def read_atoms(self, i):
        self.atoms = {}
        i += 3
        while 'Total' not in self.file_data[i]:
            self.atoms[self.file_data[i].split()[0]] = list(
                map(float, self.file_data[i].split()[1:]))
            i += 1
        return i + 1

    def read_orb_decomposition(self, i):
        self.orb_for_decomposition = {}
        i += 2
        while ':' in self.file_data[i]:
            self.orb_for_decomposition[self.file_data[i].split()[0]] = \
                self.file_data[i].replace('[', '').replace(']', '')\
                .replace(',', '').replace('\'', '').split()[2:]
            i += 1
        return i

    def prepare_dicts(self, atom_1, atom_2, key):
        if self.data[key] is None:
            self.data[key] = {}
        if atom_1 not in self.data[key]:
            self.data[key][atom_1] = {}
        if atom_2 not in self.data[key][atom_1]:
            self.data[key][atom_1][atom_2] = {}

    def format_iso(self, i):
        return float(self.file_data[i].split()[1])

    def format_aniso(self, i):
        tmp = []
        for _ in range(0, 3):
            i += 1
            tmp.append(list(map(float, self.file_data[i]
                                .replace('[', '').replace(']', '').split())))
        return tmp

    def format_dmi(self, i):
        return tuple(map(float, self.file_data[i].replace('(', '').replace(')', '')
                         .replace('[Testing!]', '').split()[1:]))

    def format_distance(self, i):
        i += 1
        tmp = self.file_data[i].replace(
            '(', '').replace(')', '').replace(',', '').split()
        self.atom_1 = tmp[0]
        self.atom_2 = tmp[1]
        self.R = tuple(map(int, tmp[2:5]))
        distance = float(tmp[9])
        if self.atom_1 not in self.magnetic_atoms:
            self.magnetic_atoms.append(self.atom_1)
        if self.atom_2 not in self.magnetic_atoms:
            self.magnetic_atoms.append(self.atom_2)
        self.list_for_sort.append((self.atom_1, self.atom_2, self.R, distance))
        self.prepare_dicts(self.atom_1, self.atom_2,
                           key=self.keys['distance'])
        return distance

    def read_exchange(self, i):
        i += 2
        while i < len(self.file_data):
            for name in self.names_to_read:
                if (self.names_to_read[name] and
                        self.names_to_read[name] in self.file_data[i]):
                    if name != 'distance':
                        self.prepare_dicts(self.atom_1,
                                           self.atom_2,
                                           key=self.keys[name])
                    self.data[self.keys[name]][self.atom_1][self.atom_2][self.R]\
                        = self.format_functions[name](i)
            i += 1
        return i

    # TODO write sorting with template after the template will be done
    def filter(self,
               distance=None,
               number=None,
               template=None,
               from_scratch=False):
        if from_scratch:
            self.data = deepcopy(self.data_nofilter)

        # TODO
        if template is not None:
            pass
        if number is not None:
            tmp = deepcopy(self.data)
            self.reset_data()
            for d in range(0, len(tmp)):
                if tmp[d] is not None:
                    for i in range(0, min(number, len(self.list_for_sort))):
                        atom_1 = self.list_for_sort[i][0]
                        atom_2 = self.list_for_sort[i][1]
                        R = self.list_for_sort[i][2]
                        self.prepare_dicts(atom_1,
                                           atom_2,
                                           key=d)
                        self.data[d][atom_1][atom_2][R] =\
                            tmp[d][atom_1][atom_2][R]

        if distance is not None:
            tmp = deepcopy(self.data)
            for d in range(0, len(tmp)):
                if tmp[d] is not None:
                    for atom_1 in tmp[d]:
                        for atom_2 in tmp[d][atom_1]:
                            for R in tmp[d][atom_1][atom_2]:
                                if tmp[
                                    self.keys['distance']
                                ][atom_1][atom_2][R] > distance:
                                    del self.data[d][atom_1][atom_2][R]
                            if not self.data[d][atom_1][atom_2]:
                                del self.data[d][atom_1][atom_2]
                        if not self.data[d][atom_1]:
                            del self.data[d][atom_1]
                    if not self.data[d]:
                        self.data[d] = None

        self.update_attributes()
