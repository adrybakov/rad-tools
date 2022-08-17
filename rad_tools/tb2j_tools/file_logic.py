from copy import deepcopy


class ExchangeModel:
    """
    Class containing the whole information extracted from TB2J *.out

    Keeping resolved data from the file
    and provide a number of functions to process it.

    Parameters
    ----------
    filename : str
        Path to the TB2J *.out file.

    Attributes
    ----------



    Methods
    -------
    filter(distance=None, number=None, template=None, from_scratch=False)
        Filter all present exchange parameters
        with respect to the provided conditions

    """

    def __init__(self, filename: str) -> None:

        self.cell = []
        self.atoms = {}
        self.magnetic_atoms = []
        self.data = []

        self.__major_sep = '=' * 90 + '\n'
        self.__minor_sep = '-' * 88 + '\n'
        self.__garbage = str.maketrans({'(': None,
                                        ')': None,
                                        '[': None,
                                        ']': None,
                                        ',': None,
                                        '\'': None})

        self.headers_to_functions = {
            'Cell (Angstrom):\n': self._read_cell,
            'Atoms:  \n': self._read_atoms,
            'Orbitals used in decomposition: \n': self._read_orb_decomposition,
            'Exchange: \n': self._read_exchange
        }

        self.__list_for_sort = []
        self.__index = {}
        self.__names_to_read = {'iso': 'J_iso:',
                                'aniso': 'J_ani:',
                                'dmi': 'DMI:',
                                'distance': self.__minor_sep}
        self.__format_functions = {'iso': self._format_iso,
                                   'aniso': self._format_aniso,
                                   'dmi': self._format_dmi,
                                   'distance': self._format_distance}
        i = 0
        for key in self.__names_to_read:
            self.__index[key] = i
            self.data.append({})
            i += 1
        with open(filename, 'r') as file:
            self._file_data = file.readlines()
        i = 0
        while i < len(self._file_data):
            if self._file_data[i] in self.headers_to_functions:
                i = self.headers_to_functions[self._file_data[i]](i)
            else:
                i += 1
        self._update_attributes()
        self.__data_nofilter = deepcopy(self.data)

    def _update_attributes(self):
        self.iso = self.data[self.__index['iso']]
        self.aniso = self.data[self.__index['aniso']]
        self.dmi = self.data[self.__index['dmi']]
        self.distance = self.data[self.__index['distance']]

    def _reset_data(self):
        for d in range(0, len(self.data)):
            self.data[d] = {}

    def _read_cell(self, i):
        i += 1
        self.cell = [
            list(map(float, self._file_data[i].split())),
            list(map(float, self._file_data[i + 1].split())),
            list(map(float, self._file_data[i + 2].split()))
        ]
        return i + 3

    def _read_atoms(self, i):
        i += 3
        while 'Total' not in self._file_data[i]:
            self.atoms[self._file_data[i].split()[0]] = list(
                map(float, self._file_data[i].split()[1:]))
            i += 1
        return i + 1

    def _read_orb_decomposition(self, i):
        self.orb_for_decomposition = {}
        i += 2
        while ':' in self._file_data[i]:
            self.orb_for_decomposition[self._file_data[i].split()[0]] =\
                self._file_data[i].translate(self.__garbage).split()[2:]
            i += 1
        return i

    def _prepare_dicts(self, atom_1, atom_2, key):
        if atom_1 not in self.data[key]:
            self.data[key][atom_1] = {}
        if atom_2 not in self.data[key][atom_1]:
            self.data[key][atom_1][atom_2] = {}

    def _format_iso(self, i):
        return float(self._file_data[i].split()[1])

    def _format_aniso(self, i):
        tmp = []
        for _ in range(0, 3):
            i += 1
            tmp.append(list(map(float, self._file_data[i]
                                .translate(self.__garbage).split())))
        return tmp

    def _format_dmi(self, i):
        return tuple(map(float, self._file_data[i].translate(self.__garbage)
                         .replace('Testing!', '').split()[1:]))

    def _format_distance(self, i):
        i += 1

        tmp = self._file_data[i].translate(self.__garbage).split()
        atom_1 = tmp[0]
        atom_2 = tmp[1]
        R = tuple(map(int, tmp[2:5]))
        distance = float(tmp[9])

        if atom_1 not in self.magnetic_atoms:
            self.magnetic_atoms.append(atom_1)
        if atom_2 not in self.magnetic_atoms:
            self.magnetic_atoms.append(atom_2)

        self.__list_for_sort.append((atom_1, atom_2, R, distance))
        return atom_1, atom_2, R, distance

    def _read_exchange(self, i):
        i += 2
        while i < len(self._file_data):
            for name in self.__names_to_read:
                if (self.__names_to_read[name] and
                        self.__names_to_read[name] in self._file_data[i]):

                    if name == 'distance':
                        atom_1, atom_2, R, tmp =\
                            self.__format_functions[name](i)
                    else:
                        tmp = self.__format_functions[name](i)

                    self._prepare_dicts(atom_1,
                                        atom_2,
                                        key=self.__index[name])
                    self.data[self.__index[name]][atom_1][atom_2][R] = tmp
            i += 1
        return i

    # TODO write sorting with template after the template will be done
    def filter(self,
               distance=None,
               number=None,
               template=None,
               from_scratch=False):
        if from_scratch:
            self.data = deepcopy(self.__data_nofilter)

        # TODO
        if template is not None:
            pass

        if number is not None:
            tmp = deepcopy(self.data)
            self._reset_data()
            for d in range(0, len(tmp)):
                if tmp[d] is not None:
                    for i in range(0, min(number, len(self.__list_for_sort))):
                        atom_1 = self.__list_for_sort[i][0]
                        atom_2 = self.__list_for_sort[i][1]
                        R = self.__list_for_sort[i][2]
                        self._prepare_dicts(atom_1,
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
                                    self.__index['distance']
                                ][atom_1][atom_2][R] > distance:
                                    del self.data[d][atom_1][atom_2][R]
                            if not self.data[d][atom_1][atom_2]:
                                del self.data[d][atom_1][atom_2]
                        if not self.data[d][atom_1]:
                            del self.data[d][atom_1]
                    if not self.data[d]:
                        self.data[d] = None

        self._update_attributes()
