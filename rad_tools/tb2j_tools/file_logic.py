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
        self.iso = None
        self.aniso = None
        self.dmi = None
        self.distance = None
        with open(filename, 'r') as file:
            self.data = file.readlines()
        i = 0
        print(self.data[2])
        while i < len(self.data):
            if self.data[i] in self.headers_to_functions:
                i = self.headers_to_functions[self.data[i]](i)
            else:
                i += 1

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
        if self.iso is None:
            self.iso = {}
        if self.aniso is None:
            self.aniso = {}
        if self.dmi is None:
            self.dmi = {}
        if self.distance is None:
            self.distance = {}
        if atom_1 not in self.iso:
            self.iso[atom_1] = {}
            self.aniso[atom_1] = {}
            self.dmi[atom_1] = {}
            self.distance[atom_1] = {}
        if atom_2 not in self.iso[atom_1]:
            self.iso[atom_1][atom_2] = {}
            self.aniso[atom_1][atom_2] = {}
            self.dmi[atom_1][atom_2] = {}
            self.distance[atom_1][atom_2] = {}

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
                if 'J_iso:' in self.data[i]:
                    self.iso[atom_1][atom_2][R] = self.format_iso(i)
                if 'J_ani:' in self.data[i]:
                    self.aniso[atom_1][atom_2][R] = self.format_aniso(i)
                if 'DMI:' in self.data[i]:
                    self.dmi[atom_1][atom_2][R] = self.format_dmi(i)
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
                    self.prepare_dicts(atom_1, atom_2)
                    self.distance[atom_1][atom_2][R] = distance
            i += 1
        return i
