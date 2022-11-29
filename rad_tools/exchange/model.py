from copy import deepcopy

import numpy as np

from rad_tools.routines import exchange_from_matrix, exchange_to_matrix


class Bond:
    r"""
    Class with implemented logic for one bond.

    Parameters
    ----------

    iso : int or float, default None
        Value of isotropic exchange parameter in meV. 

    aniso : 3 x 3 np.ndarray, None
        3 x 3 matrix of symmetric anisotropic exchange in meV. 

    dmi : 3 x 1 np.ndarray, None
        Dzyaroshinsky-Moria interaction vector :math:`(D_x, D_y, D_z)` in meV. 

    matrix : 3 x 3 np.ndarray, None
        Exchange matrix in meV. If ``matrix`` is specified then ``iso`` ,
        ``aniso`` and ``dmi`` will be ignored and derived from ``matrix`` .
        If ``matrix`` is not specified then it will be derived from
        ``iso`` , ``aniso`` and ``dmi`` .

    distance : float, default 0
        Length of the bond.

    Attributes
    ----------

    iso : float
        Value of isotropic exchange parameter in meV. If ``iso`` is
        not specified then it will be 0.

        Matrix form: ::

            [[J, 0, 0],
             [0, J, 0],
             [0, 0, J]]

    aniso : 3 x 3 np.ndarray of floats
        3 x 3 matrix of symmetric anisotropic exchange in meV. If ``aniso``
        is not specified then it will be filled with zeros.

        Matrix form: ::

            [[J_xx, J_xy, J_xz],
             [J_xy, J_yy, J_yz],
             [J_xz, J_yz, J_zz]]

    dmi : 3 x 1 np.ndarray of floats
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV. If ``dmi``
        is not specified then it will be filled with zeros.

        Vector form: ::

            [D_x, D_y, D_z]

        Matrix form: ::

            [[0, D_z, -D_y],
             [-D_z, 0, D_x],
             [D_y, -D_x, 0]]

    matrix : 3 x 3 np.ndarray of floats
        Exchange matrix in meV. If ``matrix`` is specified then ``iso`` ,
        ``aniso`` and ``dmi`` will be ignored and derived from ``matrix`` .
        If ``matrix`` is not specified then it will be derived from
        ``iso`` , ``aniso`` and ``dmi`` .

        Matrix form ::

            [[J_xx, J_xy, J_xz],
             [J_yx, J_yy, J_yz],
             [J_zx, J_zy, J_zz]]

    distance : float
        Length of the bond.

    """

    distance_tolerance = 1E-05

    def __init__(self,
                 iso=None,
                 aniso=None,
                 dmi=None,
                 matrix=None,
                 distance=None) -> None:
        self.iso = 0.
        self.aniso = np.zeros((3, 3), dtype=float)
        self.dmi = np.zeros(3, dtype=float)
        self.matrix = np.zeros((3, 3), dtype=float)
        self.dis = 0.

        if matrix is not None:
            self.matrix = np.array(matrix, dtype=float)
        if (self.matrix == np.zeros((3, 3))).all():
            if iso is not None:
                self.iso = float(iso)
            if aniso is not None:
                self.aniso = np.array(aniso, dtype=float)
            if dmi is not None:
                self.dmi = np.array(dmi, dtype=float)
            self.matrix = exchange_to_matrix(self.iso, self.aniso, self.dmi)
            if distance is not None:
                self.dis = float(distance)
        # To ensure the correct decomposition
        self.iso, self.aniso, self.dmi = exchange_from_matrix(self.matrix)

    def __add__(self, other):
        iso = (self.iso + other.iso)
        aniso = (self.aniso + other.aniso)
        dmi = (self.dmi + other.dmi)
        if self.dis - other.dis < self.distance_tolerance:
            dis = self.dis
        else:
            raise ValueError('Two bonds have different distance, could hot add. '
                             f'(Tolerance: {self.distance_tolerance})')
        return Bond(iso=iso, aniso=aniso, dmi=dmi, distance=dis)

    def __sub__(self, other):
        iso = (self.iso - other.iso)
        aniso = (self.aniso - other.aniso)
        dmi = (self.dmi - other.dmi)
        if self.dis - other.dis < self.distance_tolerance:
            dis = self.dis
        else:
            raise ValueError('Two bonds have different distance, '
                             'could hot subtraction. '
                             f'(Tolerance: {self.distance_tolerance})')
        return Bond(iso=iso, aniso=aniso, dmi=dmi, distance=dis)

    def __mul__(self, number):
        iso = self.iso * number
        aniso = self.aniso * number
        dmi = self.dmi * number
        return Bond(iso=iso, aniso=aniso, dmi=dmi, distance=self.dis)

    def __rmul__(self, number):
        iso = number * self.iso
        aniso = number * self.aniso
        dmi = number * self.dmi
        return Bond(iso=iso, aniso=aniso, dmi=dmi, distance=self.dis)


class ExchangeModel:
    r"""
    Class with the logic for the exchange model.

    Attributes
    ----------

    cell : 3 x 3 np.ndarray of floats
        3 x 3 matrix of lattice vectors in Angstrom. ::

            [[a_x  a_y  a_z],
             [b_x  b_y  b_z],
             [c_x  x_y  c_z]]

    magnetic_atoms : dict
       Dictionary with keys : str - marks of atoms and value : 1 x 3 np.ndarray
       - coordinate of the atom in Angstroms.

    bonds : dict
        Dictionary of bonds. ::

            {Atom_1: {Atom_2: {R: bond, ...}, ...}, ...}
    """

    def __init__(self) -> None:
        self.cell = np.zeros((3, 3), dtype=float)
        self.magnetic_atoms = {}
        self.bonds = {}

    def add_atom(self, name, x, y, z):
        r"""
        Add magnetic atom to the model.

        Parameters
        ----------
        name : str
            Mark for the atom. Note: if an atom with the same mark already
            exists in ``magnetic_atoms`` then it will be rewritten.

        x : int or float
            x coordinate of the atom, in Angstroms.

        y : int or float
            y coordinate of the atom, in Angstroms.

        z : int or float
            z coordinate of the atom, in Angstroms.
        """
        self.magnetic_atoms[name] = (float(x), float(y), float(z))

    def remove_atom(self, name):
        r"""
        Remove magnetic atom from the model

        Note: this method will remove atom from ``magnetic_atoms`` and all the
        bonds, which starts or ends in this atom.

        Parameters
        ----------
        name : str
            Mark for the atom.
        """
        if name in self.magnetic_atoms:
            del self.magnetic_atoms[name]
        if name in self.bonds:
            del self.bonds[name]
        for atom1 in self.bonds:
            if name in self.bonds[atom1]:
                del self.bonds[atom1][name]

    def add_bond(self, bond, atom1, atom2, R):
        r"""
        Add one bond to the model.

        Parameters
        ----------
        bond : :py:class:`Bond`
            An instance of :py:class:`Bond` class with the information about
            exchange parameters.
        atom1 : str
            Name of atom1 in (0, 0, 0) unit cell.
        atom2 : str
            Name of atom2 in R unit cell.
        R : tuple of ints
            Radius vector of the unit cell for atom2.
        """
        if atom1 not in self.bonds:
            self.bonds[atom1] = {}
        if atom2 not in self.bonds[atom1]:
            self.bonds[atom1][atom2] = {}
        self.bonds[atom1][atom2][R] = bond

    def remove_bond(self, atom1, atom2, R):
        r"""
        Remove one bond from the model.

        Parameters
        ----------
        atom1 : str
            Name of atom1 in (0, 0, 0) unit cell.

        atom2 : str
            Name of atom2 in R unit cell.

        R : tuple of ints
            Radius vector of the unit cell for atom2.
        """
        if atom1 in self.bonds and \
            atom2 in self.bonds[atom1] and \
                R in self.bonds[atom1][atom2]:
            del self.bonds[atom1][atom2][R]
            if not self.bonds[atom1][atom2]:
                del self.bonds[atom1][atom2]
            if not self.bonds[atom1]:
                del self.bonds[atom1]

    def get_lattice_vectors_length(self):
        r"""
        Getter for the lattice vectors length.

        Returns
        -------
        a : float
            Length of lattice vector a.

        b : float
            Length of lattice vector b.

        c : float
            Length of lattice vector c.
        """
        a = np.sqrt(np.sum(self.cell[0]**2))
        b = np.sqrt(np.sum(self.cell[1]**2))
        c = np.sqrt(np.sum(self.cell[2]**2))
        return a, b, c

    def get_atom_coordinates(self, atom1, atom2, R):
        r"""
        Getter for the pair of atom coordinates.

        Parameters
        ----------
        atom1 : str
            Name of atom1 in (0, 0, 0) unit cell.

        atom2 : str
            Name of atom2 in R unit cell.

        R : tuple of ints
            Radius vector of the unit cell for atom2.

        Returns
        -------
        x1 : float
            x coordinate of atom1 in the cell (0 0 0).

        y1 : float
            y coordinate of atom1 in the cell (0 0 0).

        z1 : float
            z coordinate of atom1 in the cell (0 0 0).

        x2 : float
            x coordinate of atom2 in the cell R.

        y2 : float
            y coordinate of atom2 in the cell R.

        z2 : float
            z coordinate of atom2 in the cell R.
        """
        x1 = self.magnetic_atoms[atom1][0]
        y1 = self.magnetic_atoms[atom1][1]
        z1 = self.magnetic_atoms[atom1][2]

        x2 = np.sum((R * self.cell.T).T, axis=0)[0] +\
            self.magnetic_atoms[atom2][0]
        y2 = np.sum((R * self.cell.T).T, axis=0)[1] +\
            self.magnetic_atoms[atom2][1]
        z2 = np.sum((R * self.cell.T).T, axis=0)[2] +\
            self.magnetic_atoms[atom2][2]
        return x1, y1, z1, x2, y2, z2

    def get_bond_coordinate(self, atom1, atom2, R):
        r"""
        Getter for the middle point of the bond.

        Parameters
        ----------
        atom1 : str
            Name of atom1 in (0, 0, 0) unit cell.

        atom2 : str
            Name of atom2 in R unit cell.

        R : tuple of ints
            Radius vector of the unit cell for atom2.

        Returns
        -------
        x : float
            x coordinate of the point in the middle of bond.

        y : float
            y coordinate of the point in the middle of bond.

        z : float
            z coordinate of the point in the middle of bond.

        """
        x1, y1, z1, x2, y2, z2 = self.get_atom_coordinates(atom1, atom2, R)

        return (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2

    def get_space_dimensions(self):
        r"""
        Getter for the sample size.

        Return the maximum coordinate for all 3 axes.
        Maximumm is taken for all directions.

        Returns
        -------
        X : float
            Maximum x coordinate.

        Y : float
            Maximum y coordinate.

        Z : float
            Maximum z coordinate.
        """
        X, Y, Z = 0, 0, 0
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:
                    x1, y1, z1, x2, y2, z2 = self.get_atom_coordinates(atom1,
                                                                       atom2,
                                                                       R)
                    X = max(abs(x1), abs(x2), X)
                    Y = max(abs(y1), abs(y2), Y)
                    Z = max(abs(z1), abs(z2), Z)
        return X, Y, Z

    def get_cells(self):
        r"""
        Getter for the list of cells.

        Return the list of R, specifying all cell that are present in the model.

        Returns
        -------
        cells : list of tuples of ints
            List of R.
        """
        cells = set()
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:
                    cells.add(R)
        return list(cells)

    def get_bond_list(self):
        r"""
        Getter for the list of bonds from the model.

        Returns
        -------
        bond_list : list
            List with all the bonds from the model.::

                [(atom1, atom2, R), ...]
        """
        bond_list = []
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:
                    bond_list.append((atom1, atom2, R))
        return bond_list

    # TODO It is ugly, redo in a beautifull way.
    def remove_double_bonds(self):
        r"""
        Remove double bonds.

        If for atom pair atom1, atom2 exist bond 1-2 and 2-1 they will be merged.
        All values will be calcuated as (1-2 + 2-1) / 2.
        Note: this method is not modifying the instance at which it is called.
        It will create a new instance with sorted ``bonds`` and all the other
        attributes will be copied (through deepcopy).

        Returns
        -------
        undoubled_model : ExchangeModel
            Exchange model after reoving double bonds.
        """
        unbounded_model = deepcopy(self)
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                if atom1 in unbounded_model.bonds and\
                   atom2 in unbounded_model.bonds[atom1] and\
                   (0, 0, 0) in unbounded_model.bonds[atom1][atom2] and\
                   atom2 in unbounded_model.bonds and\
                   atom1 in unbounded_model.bonds[atom2] and\
                   (0, 0, 0) in unbounded_model.bonds[atom2][atom1]:
                    unbounded_model.bonds[atom1][atom2][(0, 0, 0)] = (unbounded_model.bonds[atom1][atom2][(0, 0, 0)] +
                                                                      unbounded_model.bonds[atom2][atom1][(0, 0, 0)]) * 0.5
                    unbounded_model.remove_bond(atom2, atom1, (0, 0, 0))
        return unbounded_model

    def filter(self,
               max_distance=None,
               min_distance=None,
               template=None,
               R_vector=None):
        r"""
        Filter the exchange entries based on the given conditions.

        The result will be defined by logical conjugate of the conditions.
        Saying so the filtering will be performed for each given condition
        one by one.
        Note: this method is not modifying the instance at which it is called.
        It will create a new instance with sorted :py:attr:`bonds` and all the other
        attributes will be copied (through :py:func:`deepcopy`).

        Parameters
        ----------

        max_distance : float or int, optional
            Distance for sorting, the condition is <<less or equal>>.

        min_distance : float or int, optional
            Distance for sorting, the condition is <<more or equal>>.

        template : list
            List of pairs, which will remain in the model. ::

                [(atom1, atom2, R), ...]

        R_vector : tuple of ints or list of tuples of ints
            Tuple of 3 integers or list of tuples, specifying the R vectors,
            which will be kept after filtering.

        Returns
        -------

        filtered_model : ExchangeModel
            Exchange model after filtering.
        """
        filtered_model = deepcopy(self)

        if template is not None:
            filtered_model.bonds = {}
            for atom1, atom2, R in template:
                filtered_model.add_bond(self.bonds[atom1][atom2][R],
                                        atom1, atom2, R)
        bonds_for_removal = set()
        if R_vector is not None:
            if type(R_vector) == tuple:
                R_vector = {R_vector}
            elif type(R_vector) == list:
                R_vector = set(R_vector)
        if min_distance is not None or\
           max_distance is not None or\
           R_vector is not None:
            for atom1 in filtered_model.bonds:
                for atom2 in filtered_model.bonds[atom1]:
                    for R in filtered_model.bonds[atom1][atom2]:

                        dis = filtered_model.bonds[atom1][atom2][R].dis

                        if max_distance is not None and dis > max_distance:
                            bonds_for_removal.add((atom1, atom2, R))

                        if min_distance is not None and dis < min_distance:
                            bonds_for_removal.add((atom1, atom2, R))

                        if R_vector is not None and not R in R_vector:
                            bonds_for_removal.add((atom1, atom2, R))

        for atom1, atom2, R in bonds_for_removal:
            filtered_model.remove_bond(atom1, atom2, R)
        return filtered_model


class ExchangeModelTB2J(ExchangeModel):
    r"""
    Class containing the exchange model extracted from TB2J output file.

    Parameters
    ----------

    filename : str
        Path to the TB2J output file.
    """

    # Flags for reading from TB2J file
    _major_sep = '=' * 90
    _minor_sep = '-' * 88
    _garbage = str.maketrans({'(': None,
                              ')': None,
                             '[': None,
                              ']': None,
                              ',': None,
                              '\'': None})
    _cell_flag = 'Cell (Angstrom):'
    _atoms_flag = 'Atoms:'
    _atom_end_flag = 'Total'
    _exchange_flag = 'Exchange:'
    _iso_flag = 'J_iso:'
    _aniso_flag = 'J_ani:'
    _dmi_flag = 'DMI:'

    # Store info about all atoms, not only magnetic ones
    # {mark: (x, y, z), ...}
    _atoms = {}

    def __init__(self, filename) -> None:
        super().__init__()
        file = open(filename, 'r')
        line = True
        self.file_order = []

        # Read everything before exchange
        while line:
            line = file.readline()

            # Read cell
            if line and self._cell_flag in line:
                self.cell = np.array([
                    list(map(float, file.readline().split())),
                    list(map(float, file.readline().split())),
                    list(map(float, file.readline().split()))
                ])

            # Read atoms
            if line and self._atoms_flag in line:
                line = file.readline()
                line = file.readline()
                line = file.readline().split()
                while line and self._atom_end_flag not in line:
                    self._atoms[line[0]] = tuple(map(float, line[1:4]))
                    line = file.readline().split()

            # Check if the exchange section is reached
            if line and self._exchange_flag in line:
                break

        # Read exchange
        while line:
            while line and self._minor_sep not in line:
                line = file.readline()
            line = file.readline().translate(self._garbage).split()
            atom1 = line[0]
            atom2 = line[1]
            R = tuple(map(int, line[2:5]))
            self.file_order.append((atom1, atom2, R))
            distance = float(line[-1])
            iso = None,
            aniso = None
            dmi = None
            while line and self._minor_sep not in line:
                line = file.readline()

                # Read isotropic exchange
                if line and self._iso_flag in line:
                    iso = float(line.split()[-1])

                # Read anisotropic exchange
                if line and self._aniso_flag in line:
                    aniso = np.array([
                        list(map(float, file.readline().translate(
                            self._garbage).split())),
                        list(map(float, file.readline().translate(
                            self._garbage).split())),
                        list(map(float, file.readline().translate(
                            self._garbage).split()))
                    ])

                # Read DMI
                if line and self._dmi_flag in line:
                    dmi = tuple(map(float, line.translate(
                        self._garbage).split()[-3:]))

            # Adding info from the exchange block to the ExchangeModel structure
            if atom1 not in self.magnetic_atoms:
                self.add_atom(atom1, *self._atoms[atom1])
            if atom2 not in self.magnetic_atoms:
                self.add_atom(atom2, *self._atoms[atom2])
            bond = Bond(iso=iso, aniso=aniso, dmi=dmi, distance=distance)
            self.add_bond(bond, atom1, atom2, R)

    # TODO Think about the class type problem
    def filter(self,
               max_distance=None,
               min_distance=None,
               template=None,
               R_vector=None):
        r"""
        Filter the exchange entries based on the given conditions.

        Call :py:meth:`ExchangeModel.filter` method from parent class and
        update :py:attr:`file_order`.
        """
        filtered_model = deepcopy(self)
        result = super().filter(max_distance, min_distance, template, R_vector)
        filtered_model.file_order = []
        filtered_model.bonds = result.bonds
        filtered_model.cell = result.cell
        filtered_model.magnetic_atoms = result.magnetic_atoms
        bond_list = set(filtered_model.get_bond_list())
        for key in self.file_order:
            if key in bond_list:
                filtered_model.file_order.append(key)
        return filtered_model
