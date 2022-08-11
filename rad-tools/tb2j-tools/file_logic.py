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
        self.orb_for_decomposition = None
        with open(filename, 'r') as file:
            self.data = file.readlines()
        i = 0
        while i < len(self.data):
            if self.data[i] in self.headers_to_functions:
                i = self.headers_to_functions(self.data[i])(i + 1)
            else:
                i += 1

    def read_cell(self, i):
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
        i += 2
        while 'Total' not in self.data[i]:
            self.atoms[self.data.split()[0]] = list(
                map(float, self.data.split()[1:]))
            i += 1
        return i + 1

    def read_orb_decomposition(self, i):
        i += 1
        while ':' in self.data[i]:
            self.orb_for_decomposition[self.data.split()[0]] = \
                self.data[i].replace('[', '').replace(']', '')\
                .replace(',', '').replace('\'', '').split()
            i += 1
        return i

    def read_exchange(self, i):
        i += 1
        while i < len(self.data):
            while self.minor_sep not in self.data[i]:
                pass
            else:
                pass
        return i
