from copy import deepcopy
from typing import Union, Tuple
import numpy as np

from rad_tools.tb2j_tools.template_logic import ExchangeTemplate


class Bond:
    """
    Class with implemented logic for one bond

    Parameters
    ----------
    iso : float
        Value of isotropic exchange parameter in meV. If `iso` is 
        not specified then it will be 0.
        J
        Matrix from:

             J | 0 | 0 
            ---|---|---
             0 | J | 0 
            ---|---|---
             0 | 0 | J 
    aniso : 3 x 3 np.ndarray of floats
        3 x 3 matrix of symmetric anisotropic exchange in meV. If `aniso` 
        is not specified then it will be filled with zeros.
            J_xx | J_xy | J_xz
            -----|------|-----
            J_xy | J_yy | J_yz
            -----|------|-----
            J_xz | J_yz | J_zz
    dmi : 3 x 1 np.ndarray of floats
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV. If `dmi` 
        is not specified then it will be filled with zeros.
        [D_x, D_y, D_z]
        Matrix form:
               0  |  D_z | -D_y
            ------|------|------
             -D_z |   0  | D_x
            ------|------|------
              D_y | -D_x |  0
    matrix : 3 x 3 np.ndarray of floats
        Exchange matrix in meV. If `matrix` is specified then `iso`, 
        `aniso` and `dmi` will be ignored and derived from `matrix`.
        If `matrix` is not specified then it will be derived from 
        `iso`, `aniso` and `dmi`.
            J_xx | J_xy | J_xz
            -----|------|-----
            J_yx | J_yy | J_yz
            -----|------|-----
            J_zx | J_zy | J_zz
    """

    def __init__(self,
                 iso=None,
                 aniso=None,
                 dmi=None,
                 matrix=None) -> None:
        self.iso = 0.
        self.aniso = np.zeros((3, 3), dtype=float)
        self.dmi = np.zeros(3, dtype=float)
        self.matrix = np.zeros((3, 3), dtype=float)
        if matrix is not None:
            self.matrix = np.array(matrix, dtype=float)
            self.from_matrix()
        else:
            if iso is not None:
                self.iso = float(iso)
            if aniso is not None:
                self.aniso = np.array(aniso, dtype=float)
            if dmi is not None:
                self.dmi = np.array(dmi, dtype=float)
            self.to_matrix()

    def to_matrix(self):
        """
        Combine isotropic, anisotropic and dmi exchange into exchange matrix.

        Parameters
        ----------
        isotropic : float
            Value of isotropic exchange parameter in meV.
        anisotropic : 3 x 3 np.ndarray of floats
            Matrix of symmetric anisotropic exchange in meV.
        dmi : 3 x 1 np.ndarray of floats
            Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV.
        """
        self.matrix = np.zeros((3, 3), dtype=float)
        self.matrix += self.aniso
        self.matrix += self.iso * np.identity(3, dtype=float)
        self.matrix += np.array([[0, self.dmi[2], -self.dmi[1]],
                                 [-self.dmi[2], 0, self.dmi[0]],
                                 [self.dmi[1], -self.dmi[0], 0]],
                                dtype=float)

    def from_matrix(self):
        """
        Decompose matrix into isotropic, anisotropic and dmi exchange.

        Parameters
        ----------
        matrix : 3 x 3 np.ndarray of floats
            Exchange matrix in meV.
        """
        symm = (self.matrix + self.matrix.T) / 2
        assym = (self.matrix - self.matrix.T) / 2
        self.dmi = np.array([assym[1][2], assym[2][0], assym[0][1]],
                            dtype=float)
        self.iso = np.trace(symm) / 3
        self.aniso = symm - self.iso * np.identity(3, dtype=float)


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
    cell : list of list of float
        3x3 matrix of lattice vectors in Angstrom.
        Directly read from "Cell (Angstrom):" section 
        of TB2J file.
        [
            [a_x, a_y, a_z],
            [b_x, b_y, b_z],
            [c_x, x_y, c_z]
        ]

    atoms : dict 
        Info read from "Atoms:" section  of TB2J file.
        Keys are the marks of atoms (str), values are the 
        list of float with the data
        {
            "Atom_1" : [x, y, z, ...],
            ...
        }

    magnetic_atoms : list of str
        list of atom's marks for which the entry in "Exchange:" 
        section of TB2J file is present.

    iso : dict
        Isotropic Exchange parameters read from 
        "Exchange:" section of TB2J file. 
        atom_1: str - mark of the atom in central unit cell.
        atom_2: str - mark of the atom in the other unit cell.
        R_v - radius vector between the central unit cell and the other unit cell.
        R_v = (a_step: int, b_step: int, c_step: int). 
        J_iso: float -  Isotropic exchange parameter 
        in the notation of TB2J in meV.
        {
            "atom_1" : 
            {
                "atom_2": 
                {
                    R_v : J_iso,
                    ...
                },
                ...
            },
            ...
        }

    aniso : dict
        Anisotropic symmetric Exchange parameters read 
        from "Exchange:" section of TB2J file. 
        Structure is the same as in iso, but with J_aniso
        instead of J_iso.
        J_aniso : list of list of float
        [
            [J_xx, J_xy, J_xz],
            [J_yx, J_yy, J_yz],
            [J_zx, J_zy, J_zz]
        ]
        Note: J_ij = J_ji.

    dmi : dict
        Dzyaloshinskii-Moriya interaction parameters read 
        from "Exchange:" section of TB2J file. 
        Structure is the same as in iso, but with D
        instead of J_iso.
        D: tuple of float - DM vector (D_x, D_y, D_z)

    distance : dict
        Distance between pairs of magnetic atoms. 
        Structure is the same as in iso, but with d instead of J_iso.
        d : float - distance between two magnetic atoms in Angstrom.

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
        self.iso = {}
        self.aniso = {}
        self.dmi = {}
        self.distance = {}

        self.__major_sep = '=' * 90 + '\n'
        self.__minor_sep = '-' * 88 + '\n'
        self.__garbage = str.maketrans({'(': None,
                                        ')': None,
                                        '[': None,
                                        ']': None,
                                        ',': None,
                                        '\'': None})

        self.__headers_to_functions = {
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

        self.__data = []
        i = 0
        for key in self.__names_to_read:
            self.__index[key] = i
            self.__data.append({})
            i += 1
        with open(filename, 'r') as file:
            self.__file_data = file.readlines()
        i = 0
        while i < len(self.__file_data):
            if self.__file_data[i] in self.__headers_to_functions:
                i = self.__headers_to_functions[self.__file_data[i]](i)
            else:
                i += 1
        self._update_attributes()
        self.__data_nofilter = deepcopy(self.__data)

    def _update_attributes(self):
        """
        Updates the values of public attributes based on the 
        state of the main <<working horse>> of this class `.__data`.
        """

        self.iso = self.__data[self.__index['iso']]
        self.aniso = self.__data[self.__index['aniso']]
        self.dmi = self.__data[self.__index['dmi']]
        self.distance = self.__data[self.__index['distance']]

    def _reset_data(self):
        """
        Reset the data in `.__data` to the empty dict.
        """

        for d in range(0, len(self.__data)):
            self.__data[d] = {}

    def _read_cell(self, i: int):
        """
        Read and format the information about cell vectors.

        `.cell` attribute will be filled with content here.

        Parameters
        ----------
        i : int
            Index of the line with Cell section header 
            (see `.__headers_to_functions`) in the `.__file_data`.

        Returns
        -------
        i : int
            Index of next line after the last line with the cell data 
            in the `.__file_data`.
        """

        i += 1
        self.cell = [
            list(map(float, self.__file_data[i].split())),
            list(map(float, self.__file_data[i + 1].split())),
            list(map(float, self.__file_data[i + 2].split()))
        ]
        return i + 3

    def _read_atoms(self, i: int):
        """
        Read and format the information about the atoms.

        `.atoms` attribute will be filled with content here.

        Parameters
        ----------
        i : int
            Index of the line with Atoms section header 
            (see `.__headers_to_functions`) in the `.__file_data`.

        Returns
        -------
        i : int
            Index of next line after the last line with the atoms data 
            in the `.__file_data` 
            (next line after the one that contains "Total" in it).
        """

        i += 3
        while 'Total' not in self.__file_data[i]:
            self.atoms[self.__file_data[i].split()[0]] = list(
                map(float, self.__file_data[i].split()[1:]))
            i += 1
        return i + 1

    def _read_orb_decomposition(self, i: int):
        """
        Read and format the information about the names 
        of orbitals for orbital decomposition.

        `.orb_for_decomposition` attribute will be filled with content here.
        May be useless since the names of orbital are meaningless in TB2J, 
        but who knows.

        Parameters
        ----------
        i : int
            Index of the line with Orbitals used in decomposition section 
            header (see `.__headers_to_functions`) in the `.__file_data`.

        Returns
        -------
        i : int
            Index of next line after the last line with the orbital 
            decomposition data in the `.__file_data`.
        """

        self.orb_for_decomposition = {}
        i += 2
        while ':' in self.__file_data[i]:
            self.orb_for_decomposition[self.__file_data[i].split()[0]] =\
                self.__file_data[i].translate(self.__garbage).split()[2:]
            i += 1
        return i

    def _prepare_dicts(self, atom_1: str, atom_2: str, key: int):
        """
        Since the `.__data` attribute is an array of nested dicts
        the dicts have to be prepared. This function does it.

        Parameters
        ----------
        atom_1 : str
            mark of the first-level atom (see `iso`). 
        atom_2 : str
            mark of the second-level atom (see `iso`). 
        key : int
            Index of the dict in `.__data` to prepare. 
            It is recommended to use `.__index` attribute with desired keyword 
            (i.e. "iso", "aniso", ...) for passing the `key` to this function.
        """

        if atom_1 not in self.__data[key]:
            self.__data[key][atom_1] = {}
        if atom_2 not in self.__data[key][atom_1]:
            self.__data[key][atom_1][atom_2] = {}

    def _format_iso(self, i: int):
        """
        Format the isotropic exchange parameter from the line of TB2J file 
        to float.

        Parameters
        ----------
        i : int
            Index of the line with isotropic exchange parameter header 
            (see `.__names_to_read`) in `.__file_data`.

        Returns
        -------
        J_iso : float
            Isotropic exchange parameter in meV. In the notation of TB2J.
        """

        return float(self.__file_data[i].split()[1])

    def _format_aniso(self, i: int):
        """
        Format the anisotropic exchange parameter from the lines of TB2J file.

        Parameters
        ----------
        i : int
            Index of the line with anisotropic exchange parameter header 
            (see `.__names_to_read`) in `.__file_data`.

        Returns
        -------
        J_aniso : list of list of float
            Symmetric anisotropic exchange 3x3 matrix in meV. 
            In the notation of TB2J.
        """

        J_aniso = []
        for _ in range(0, 3):
            i += 1
            J_aniso.append(list(map(float, self.__file_data[i]
                                .translate(self.__garbage).split())))
        return J_aniso

    def _format_dmi(self, i: int):
        """
        Format the Dzyaloshinskii-Moriya interaction parameter 
        from the lines of TB2J file.

        Parameters
        ----------
        i : int
            Index of the line with Dzyaloshinskii-Moriya interaction parameter 
            header (see `.__names_to_read`) in `.__file_data`.

        Returns
        -------
        D : tuple of float
            Dzyaloshinskii-Moriya interaction vector in meV. 
            In the notation of TB2J.
        """

        return tuple(map(float, self.__file_data[i].translate(self.__garbage)
                         .replace('Testing!', '').split()[1:]))

    def _format_distance(self, i: int):
        """
        Format the first line of TB2J file`s small exchange block.

        This function is more important then the other formatters. 
        It sets the `atom_1`, `atom_2` and `R_v` which are valid for current 
        exchange block and fill `.magnetic_atoms` and `.__list_for_sort` with
        content.

        Parameters
        ----------
        i : int
            Index of the line previous to the first line of exchange block 
            (see `.__names_to_read`) in `.__file_data`.

        Returns
        -------
        atom_1 : str
            mark of the first-level atom (see `iso`). 
        atom_2 : str
            mark of the second-level atom (see `iso`). 
        R_v : tuple of float
            Radius vector between the central unit cell and the other unit cell.
        distance : float
            Distance between pair (`atom_1`, `atom_2`) in Angstrom. 
        """

        i += 1

        tmp = self.__file_data[i].translate(self.__garbage).split()
        atom_1 = tmp[0]
        atom_2 = tmp[1]
        R_v = tuple(map(int, tmp[2:5]))
        distance = float(tmp[9])

        if atom_1 not in self.magnetic_atoms:
            self.magnetic_atoms.append(atom_1)
        if atom_2 not in self.magnetic_atoms:
            self.magnetic_atoms.append(atom_2)

        self.__list_for_sort.append((atom_1, atom_2, R_v, distance))
        return atom_1, atom_2, R_v, distance

    def _read_exchange(self, i: int):
        """
        Read and format the Exchange section.

        `.__data` attribute and corresponding public attributes will be 
        filled with content here.

        Parameters
        ----------
        i : int
            Index of the line with Exchange section header 
            (see `.__headers_to_functions`) in the `.__file_data`.

        Returns
        -------
        i : int
            Index of next line after the last line of Exchange section in the 
            `.__file_data`. This is effectivly the len(`self.__file_data`)
        """

        i += 2
        while i < len(self.__file_data):
            for name in self.__names_to_read:
                if (self.__names_to_read[name] and
                        self.__names_to_read[name] in self.__file_data[i]):

                    if name == 'distance':
                        atom_1, atom_2, R_v, tmp =\
                            self.__format_functions[name](i)
                    else:
                        tmp = self.__format_functions[name](i)

                    self._prepare_dicts(atom_1,
                                        atom_2,
                                        key=self.__index[name])
                    self.__data[self.__index[name]][atom_1][atom_2][R_v] = tmp
            i += 1
        return i

    def filter(self,
               distance: Union[float, int] = None,
               number: int = None,
               template: ExchangeTemplate = None,
               from_scratch: bool = False):
        """
        Filter the Exchange entries based on the given conditions.

        The result will be defined by logical conjugate of the conditions. 
        Saying so the filtering will be performed for each given condition 
        one by one.

        `.__data` attribute and corresponding public attributes will be 
        updated here.

        Parameters
        ----------
        distance : float | int, optional
            Distance for sorting, the condition is the less or equal.
        number : int, optional
            The exact amount of interaction pairs to be kept. The pairs are 
            sorted by distance and basically have exactly the same order as in 
            TB2J file (from top to down). It will be interpreted as the amount 
            of exchange blocks in TB2J file to be kept.
        template : ExchangeTemplate
            template instance for filtering neighbors
        from_scratch : bool, default=False
            If True then the data will be restored from the original file 
            (through `.__data_nofilter` attribute) before sorting.
        """

        if from_scratch:
            self.__data = deepcopy(self.__data_nofilter)

        tmp_list = []
        if template is not None and number is not None:
            for i in self.__list_for_sort[:number]:
                if (i not in tmp_list and
                        (i[0], i[1], i[2]) in template.plained_template):
                    tmp_list.append(i)
        elif template is not None:
            for i in template.plained_template:
                if i not in tmp_list:
                    tmp_list.append(i)
        elif number is not None:
            for i in self.__list_for_sort[:number]:
                if i not in tmp_list:
                    tmp_list.append(i)

        if template is not None or number is not None:
            tmp = deepcopy(self.__data)
            self._reset_data()
            for d in range(0, len(tmp)):
                if tmp[d]:
                    for i in range(0, len(tmp_list)):
                        atom_1 = tmp_list[i][0]
                        atom_2 = tmp_list[i][1]
                        R_v = tmp_list[i][2]
                        self._prepare_dicts(atom_1,
                                            atom_2,
                                            key=d)
                        self.__data[d][atom_1][atom_2][R_v] =\
                            tmp[d][atom_1][atom_2][R_v]

        if distance is not None:
            tmp = deepcopy(self.__data)
            for d in range(0, len(tmp)):
                if tmp[d]:
                    for atom_1 in tmp[d]:
                        for atom_2 in tmp[d][atom_1]:
                            for R_v in tmp[d][atom_1][atom_2]:
                                if tmp[
                                    self.__index['distance']
                                ][atom_1][atom_2][R_v] > distance:
                                    del self.__data[d][atom_1][atom_2][R_v]
                            if not self.__data[d][atom_1][atom_2]:
                                del self.__data[d][atom_1][atom_2]
                        if not self.__data[d][atom_1]:
                            del self.__data[d][atom_1]
                    if not self.__data[d]:
                        self.__data[d] = {}

        self._update_attributes()
