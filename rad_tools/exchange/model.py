r"""
Exchange model.
"""

from copy import deepcopy
from math import sqrt

import numpy as np

from rad_tools.exchange.template import ExchangeTemplate


class Bond:
    r"""
    Exchange bond.

    If ``matrix`` is specified then ``iso``, ``aniso`` and ``dmi`` will 
    be ignored and derived from ``matrix``. If ``matrix`` is not specified 
    then it will be derived from ``iso``, ``aniso`` and ``dmi``.

    Parameters
    ----------
    iso : int or float, default None
        Value of isotropic exchange parameter in meV.
    aniso : 3 x 3 array, None
        3 x 3 matrix of symmetric anisotropic exchange in meV.
    dmi : 3 x 1 array, None
        Dzyaroshinsky-Moria interaction vector :math:`(D_x, D_y, D_z)` in meV.
    matrix : 3 x 3 array, None
        Exchange matrix in meV. If ``matrix`` is specified then ``iso`` ,
        ``aniso`` and ``dmi`` will be ignored and derived from ``matrix`` .
        If ``matrix`` is not specified then it will be derived from
        ``iso`` , ``aniso`` and ``dmi`` .
    distance : float, default 0
        Length of the bond.

    Attributes
    ----------
    matrix : 3 x 3 array of floats
        Exchange matrix in meV. 

        Matrix form 

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_yx, J_yy, J_yz],
             [J_zx, J_zy, J_zz]]

    iso : float
    aniso : 3 x 3 array of floats
    dmi : 3 x 1 array of floats
    dis : float
        Length of the bond.
    symm_matrix : 3 x 3 array of floats
    asymm_matrix : 3 x 3 array of floats
    dmi_module : float
    dmi_vs_iso : float
    """

    distance_tolerance = 1E-10

    def __init__(self,
                 iso=None,
                 aniso=None,
                 dmi=None,
                 matrix=None,
                 distance=None) -> None:
        self._matrix = np.zeros((3, 3), dtype=float)
        self._iso = 0.
        self._aniso = np.zeros((3, 3), dtype=float)
        self._dmi = np.zeros(3, dtype=float)
        self.dis = 0.

        if matrix is None:
            self.iso = iso
            self.dmi = dmi
            self.aniso = aniso
        else:
            self.matrix = matrix

        if distance is not None:
            self.dis = float(distance)

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self, new_matrix):
        if new_matrix is None:
            new_matrix = np.zeros((3, 3), dtype=float)
        new_matrix = np.array(new_matrix)
        if new_matrix.shape != (3, 3):
            raise ValueError("Matrix shape have to be equal to (3, 3)")
        self._matrix = new_matrix

    @property
    def symm_matrix(self):
        r"""
        Symmetric part of exchange matrix.
        """

        return (self.matrix + self.matrix.T) / 2

    @property
    def asymm_matrix(self):
        """
        Asymmetric part of exchange matrix.
        """

        return (self.matrix - self.matrix.T) / 2

    @property
    def iso(self):
        r"""
        Value of isotropic exchange parameter in meV. 
        If ``iso`` is not specified then it will be 0.

        Matrix form: 

        .. code-block:: python

            [[J, 0, 0],
             [0, J, 0],
             [0, 0, J]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as 

        .. math::
            J_{iso} = \dfrac{1}{3}Tr(\mathbf{J})
        """
        return np.trace(self.symm_matrix) / 3

    @iso.setter
    def iso(self, new_iso):
        if new_iso is not None:
            new_iso = new_iso - self.iso
        else:
            new_iso = -self.iso
        self.matrix = (self.matrix +
                       (new_iso) *
                       np.identity(3, dtype=float))

    @property
    def aniso(self):
        r"""
        3 x 3 matrix of symmetric anisotropic exchange in meV. 

        Matrix form: 

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_xy, J_yy, J_yz],
             [J_xz, J_yz, J_zz]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as 

        .. math::
            J_{aniso} = \mathbf{J}_{symm} - \dfrac{1}{3}Tr(\mathbf{J})
        """
        return self.symm_matrix - self.iso * np.identity(3, dtype=float)

    @aniso.setter
    def aniso(self, new_aniso):
        if new_aniso is not None:
            new_aniso = np.array(new_aniso, dtype=float)
            if new_aniso.shape != (3, 3):
                raise ValueError("Aniso matrix shape have " +
                                 "to be equal to (3, 3)")
            new_aniso = new_aniso - self.aniso
        else:
            new_aniso = -self.aniso
        self.matrix = self.matrix + new_aniso

    @property
    def dmi(self):
        r"""
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV. 

        Vector form: 

        .. code-block:: python

            [D_x, D_y, D_z]

        Matrix form: 

        .. code-block:: python

            [[0, D_z, -D_y],
             [-D_z, 0, D_x],
             [D_y, -D_x, 0]]

        Derived from antisymmetric part of exchange matrix (:math:`\mathbf{J}`).
        """
        return np.array([self.asymm_matrix[1][2],
                         self.asymm_matrix[2][0],
                         self.asymm_matrix[0][1]],
                        dtype=float)

    @dmi.setter
    def dmi(self, new_dmi):
        if new_dmi is not None:
            new_dmi = np.array(new_dmi, dtype=float)
            if new_dmi.shape != (3,):
                raise ValueError(
                    f"DMI have to be a 3 component vector. {new_dmi}")
            new_dmi = new_dmi - self.dmi
        else:
            new_dmi = -self.dmi
        dmi_matrix = np.array([[0, new_dmi[2], -new_dmi[1]],
                               [-new_dmi[2], 0, new_dmi[0]],
                               [new_dmi[1], -new_dmi[0], 0]],
                              dtype=float)
        self.matrix = self.matrix + dmi_matrix

    @property
    def dmi_module(self):
        r"""
        Length of the DMI vector in th e units of exchange interaction.
        """

        return sqrt(np.sum(self.dmi**2))

    @property
    def dmi_vs_iso(self):
        r"""
        Relative strength of DMI.

        .. math::

            \dfrac{\vert\vec{D}\vert}{\vert J_{iso}\vert}
        """
        return abs(self.dmi_module/self.iso)

    def __add__(self, other):
        if isinstance(other, Bond):
            if abs(self.dis - other.dis) < self.distance_tolerance:
                dis = self.dis
            else:
                raise ValueError('Two bonds have different distance, could hot add. '
                                 f'(Tolerance: {self.distance_tolerance})')
            matrix = self.matrix + other.matrix
            return Bond(matrix=matrix, distance=dis)
        else:
            raise TypeError

    def __sub__(self, other):
        if isinstance(other, Bond):
            if abs(self.dis - other.dis) < self.distance_tolerance:
                dis = self.dis
            else:
                raise ValueError('Two bonds have different distance, could hot add. '
                                 f'(Tolerance: {self.distance_tolerance})')
            matrix = self.matrix - other.matrix
            return Bond(matrix=matrix, distance=dis)
        else:
            raise TypeError

    def __mul__(self, number):
        matrix = self.matrix * number
        return Bond(matrix=matrix, distance=self.dis)

    def __rmul__(self, number):
        return self.__mul__(number)

    def __truediv__(self, number):
        matrix = self.matrix / number
        return Bond(matrix=matrix, distance=self.dis)


class ExchangeModel:
    r"""
    Class with the logic for the exchange model.

    Attributes
    ----------
    bonds : dict
        Dictionary of bonds.

        .. code-block:: python

            {Atom_1: {Atom_2: {R: bond, ...}, ...}, ...}
    magnetic_atoms : dict
        Dictionary with keys : str - marks of atoms and value : 1 x 3 np.ndarray
        - coordinate of the atom in Angstroms.
    noomagnetic_atoms : dist
        Dictionary with keys : str - marks of atoms and value : 1 x 3 np.ndarray
        - coordinate of the atom in Angstroms. Only non-magnetic atoms.
    cell : 3 x 3 array of floats
    a : array
    b : array
    c : array
    len_a : float
    len_b : float
    len_c : float
    unit_cell_volume : float
    bond_list : list
    cell_list : list
    space_dimensions : tuple
    """

    def __init__(self) -> None:
        self._cell = np.identity(3, dtype=float)
        self.magnetic_atoms = {}
        self.nonmagnetic_atoms = {}
        self._bond_list = []
        self.bonds = {}

    @property
    def cell(self):
        r"""
        3 x 3 matrix of lattice vectors in Angstrom.

        .. code-block:: python

            [[a_x  a_y  a_z],
             [b_x  b_y  b_z],
             [c_x  x_y  c_z]]

        Default cell is orthonormal.
        """
        return self._cell

    @cell.setter
    def cell(self, new_cell):
        new_cell = np.array(new_cell)
        if new_cell.shape != (3, 3):
            raise ValueError("Cell matrix shape have to be equal (3, 3)")
        if np.linalg.det(new_cell) == 0:
            raise ValueError("Lattice vectors are complanar.")
        self._cell = new_cell

    @property
    def a(self):
        r"""
        First lattice vector :math:`\vec{a} = (a_x, a_y, a_z)`.
        """

        return self.cell[0]

    @a.setter
    def a(self, new_a):
        new_a = np.array(new_a)
        if new_a.shape != (3,):
            raise ValueError("new_a is not a 3 component vector")
        self._cell[0] = new_a

    @property
    def b(self):
        r"""
        Second lattice vector :math:`\vec{b} = (b_x, b_y, b_z)`.
        """

        return self.cell[1]

    @b.setter
    def b(self, new_b):
        new_b = np.array(new_b)
        if new_b.shape != (3,):
            raise ValueError("new_b is not a 3 component vector")
        self._cell[1] = new_b

    @property
    def c(self):
        r"""
        Third lattice vector :math:`\vec{b} = (b_x, b_y, b_z)`.
        """

        return self.cell[2]

    @c.setter
    def c(self, new_c):
        new_c = np.array(new_c)
        if new_c.shape != (3,):
            raise ValueError("new_c is not a 3 component vector")
        self._cell[2] = new_c

    @property
    def len_a(self):
        r"""
        Length of lattice vector :math:`\vec{a}`.
        """
        return sqrt(np.sum(self.cell[0]**2))

    @property
    def len_b(self):
        r"""
        Length of lattice vector :math:`\vec{b}`.
        """
        return sqrt(np.sum(self.cell[1]**2))

    @property
    def len_c(self):
        r"""
        Length of lattice vector :math:`\vec{c}`.
        """
        return sqrt(np.sum(self.cell[2]**2))

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        .. math::

            V = \delta_{\vec{A}}\cdot(\delta_{\vec{B}}\times
            \delta_{\vec{C}})
        """

        return np.dot(self.a, np.cross(self.b, self.c))

    @property
    def bond_list(self):
        r"""
        List of bonds from the model.

        Returns
        -------
        bond_list : list
            List with all the bonds from the model.

            .. code-block:: python

                [(atom1, atom2, R), ...]
        """
        bond_list = []
        for atom1, atom2, R in self._bond_list:
            if (atom1 in self.bonds
                and atom2 in self.bonds[atom1]
                    and R in self.bonds[atom1][atom2]):
                bond_list.append((atom1, atom2, R))
        return bond_list

    @property
    def cell_list(self):
        r"""
        List of cells from the model.

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

    def add_bond(self, bond: Bond, atom1, atom2, R):
        r"""
        Add one bond to the model.

        It is important to update ``bond_list`` here, 
        since it alows to track the order of the bonds.

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
        self._bond_list.append((atom1, atom2, R))

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

    def add_atom(self, name, x, y, z):
        r"""
        Add magnetic atom to the model.

        Parameters
        ----------
        name : str
            Mark for the atom. Note: if an atom with the same mark already
            exists in ``magnetic_atoms`` then it will be rewritten.
        x : int or float
            x coordinate of the atom, in real space (Angstroms).
        y : int or float
            y coordinate of the atom, in real space (Angstroms).
        z : int or float
            z coordinate of the atom, in real space (Angstroms).
        """
        self.magnetic_atoms[name] = np.array([x, y, z])

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
        atoms_to_clear = []
        if name in self.magnetic_atoms:
            del self.magnetic_atoms[name]
        if name in self.bonds:
            del self.bonds[name]
        for atom1 in self.bonds:
            if name in self.bonds[atom1]:
                del self.bonds[atom1][name]
            if len(self.bonds[atom1]) == 0:
                atoms_to_clear.append(atom1)
        for atom in atoms_to_clear:
            del self.bonds[atom]

    def get_atom_coordinates(self, atom, R=(0, 0, 0)):
        r"""
        Getter for the pair of atom coordinates.

        Parameters
        ----------
        atom : str
            Name of atom1 in R unit cell.
        R : tuple of ints, default (0, 0, 0)
            Radius vector of the unit cell for atom2.

        Returns
        -------
        x : float
            x coordinate of atom in the cell R.
        y : float
            y coordinate of atom in the cell R.
        z : float
            z coordinate of atom in the cell R.
        """
        R_vector = np.sum((R * self.cell.T).T, axis=0)
        x = R_vector[0] + self.magnetic_atoms[atom][0]
        y = R_vector[1] + self.magnetic_atoms[atom][1]
        z = R_vector[2] + self.magnetic_atoms[atom][2]
        return x, y, z

    def get_bond_coordinates(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for the middle point of the bond.

        Parameters
        ----------
        atom1 : str
            Name of atom1 in (0, 0, 0) unit cell.
        atom2 : str
            Name of atom2 in R unit cell.
        R : tuple of ints, default (0, 0, 0)
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

        x1, y1, z1 = self.get_atom_coordinates(atom1)
        x2, y2, z2 = self.get_atom_coordinates(atom2, R)

        return (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2

    @property
    def space_dimensions(self):
        r"""
        Model minimum and maximum coordiantes.

        Return the maximum and minimum coordinates for all 3 axes.
        Takes into account only the atoms which have at least one bond.

        Returns
        -------
        x_min : float
            Minimum x coordinate.
        y_min : float
            Minimum y coordinate.
        z_min : float
            Minimum z coordinate.
        x_max : float
            Maximum x coordinate.
        y_max : float
            Maximum y coordinate.
        z_max : float
            Maximum z coordinate.
        """
        x_max, y_max, z_max = -10E10, -10E10, -10E10
        x_min, y_min, z_min = 10E10, 10E10, 10E10
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:
                    x1, y1, z1 = self.get_atom_coordinates(atom1)
                    x2, y2, z2 = self.get_atom_coordinates(atom2, R)
                    x_max = max(x1, x2, x_max)
                    y_max = max(y1, y2, y_max)
                    z_max = max(z1, z2, z_max)
                    x_min = min(x1, x2, x_min)
                    y_min = min(y1, y2, y_min)
                    z_min = min(z1, z2, z_min)
        return x_min, y_min, z_min, x_max, y_max, z_max

    def remove_double_bonds(self):
        r"""
        Remove double bonds.

        If for atom pair atom1, atom2 exist bond 1-2 and 2-1 they 
        will be merged.
        All values will be calcuated as (1-2 + 2-1) / 2.

        .. note::
            This method is modifying the instance at which it is called.
        """

        pairs = set()

        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:
                    R_mirror = (-R[0], -R[1], -R[2])
                    if (atom2 in self.bonds
                        and atom1 in self.bonds[atom2]
                            and R_mirror in self.bonds[atom2][atom1]):
                        if ((atom2, atom1, R_mirror),
                                (atom1, atom2, R)) not in pairs:
                            pairs.add(((atom1, atom2, R),
                                       (atom2, atom1, R_mirror)))
        pairs = list(pairs)
        for bond1, bond2 in pairs:
            atom1, atom2, R = bond1
            bond1 = self.bonds[atom1][atom2][R]
            atom2, atom1, R_mirror = bond2
            bond2 = self.bonds[atom2][atom1][R_mirror]

            flip = False
            index = 0
            if R[0] == 0:
                if R[1] == 0:
                    if R[2] == 0:
                        x1, y1, z1 = self.get_atom_coordinates(atom1)
                        x2, y2, z2 = self.get_atom_coordinates(atom2)
                        x = x2 - x1
                        y = y2 - y1
                        z = z2 - z1
                        if x < 0:
                            flip = True
                        elif x == 0 and y < 0:
                            flip = True
                        elif x == 0 and y == 0 and z < 0:
                            flip = True
                    else:
                        index = 2
                else:
                    index = 1
            if R[index] < R_mirror[index] or flip:
                R, R_mirror = R_mirror, R
                atom1, atom2 = atom2, atom1
            bond = (bond1 + bond2) / 2
            self.remove_bond(atom2, atom1, R_mirror)
            self.bonds[atom1][atom2][R] = bond

    def removed_double_bonds(self):
        r"""
        Remove double bonds and create a new model.

        If for atom pair atom1, atom2 exist bond 1-2 and 2-1 they 
        will be merged.
        All values will be calcuated as (1-2 + 2-1) / 2.

        .. note::
            This method is not modifying the instance at which it is called.
            It will create a new instance with merged ``bonds`` and all the 
            other attributes will be copied (through deepcopy).

        Returns
        -------
        undoubled_model : ExchangeModel
            Exchange model after reoving double bonds.
        """

        unbounded_model = deepcopy(self)
        unbounded_model.remove_double_bonds()
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

        .. note::
            This method is modifying the instance at which it is called.

        Parameters
        ----------
        max_distance : float or int, optional
            Distance for sorting, the condition is <<less or equal>>.
        min_distance : float or int, optional
            Distance for sorting, the condition is <<more or equal>>.
        template : list or :py:class:`.ExchangeTemplate`
            List of pairs, which will remain in the model. ::

                [(atom1, atom2, R), ...]

        R_vector : tuple of ints or list of tuples of ints
            Tuple of 3 integers or list of tuples, specifying the R vectors,
            which will be kept after filtering.
        """
        if isinstance(template, ExchangeTemplate):
            template = template.get_list()
        if template is not None:
            template = set(template)
        if R_vector is not None:
            if type(R_vector) == tuple:
                R_vector = {R_vector}
            elif type(R_vector) == list:
                R_vector = set(R_vector)

        bonds_for_removal = set()
        for atom1 in self.bonds:
            for atom2 in self.bonds[atom1]:
                for R in self.bonds[atom1][atom2]:

                    dis = self.bonds[atom1][atom2][R].dis

                    if max_distance is not None and dis > max_distance:
                        bonds_for_removal.add((atom1, atom2, R))

                    if min_distance is not None and dis < min_distance:
                        bonds_for_removal.add((atom1, atom2, R))

                    if R_vector is not None and R not in R_vector:
                        bonds_for_removal.add((atom1, atom2, R))
                    if template is not None and (atom1, atom2, R) not in template:
                        bonds_for_removal.add((atom1, atom2, R))

        for atom1, atom2, R in bonds_for_removal:
            self.remove_bond(atom1, atom2, R)

    def filtered(self,
                 max_distance=None,
                 min_distance=None,
                 template=None,
                 R_vector=None):
        r"""
        Create filtered exchange model based on the given conditions.

        The result will be defined by logical conjugate of the conditions.
        Saying so the filtering will be performed for each given condition
        one by one.
        Note: this method is not modifying the instance at which it is called.
        It will create a new instance with sorted :py:attr:`bonds` and all the other
        attributes will be copied (through :py:func:`deepcopy`).

        .. note::
            This method is not modifying the instance at which it is called.
            It will create a new instance with merged ``bonds`` and all the 
            other attributes will be copied (through deepcopy).

        Parameters
        ----------
        max_distance : float or int, optional
            Distance for sorting, the condition is <<less or equal>>.
        min_distance : float or int, optional
            Distance for sorting, the condition is <<more or equal>>.
        template : list or :py:class:`.ExchangeTemplate`
            List of pairs, which will remain in the model. ::

                [(atom1, atom2, R), ...]

        R_vector : tuple of ints or list of tuples of ints
            Tuple of 3 integers or list of tuples, specifying the R vectors,
            which will be kept after filtering.

        Returns
        -------
        filtered_model : :py:class:`.ExchangeModel`
            Exchange model after filtering.
        """
        filtered_model = deepcopy(self)
        filtered_model.filter(max_distance=max_distance,
                              min_distance=min_distance,
                              template=template,
                              R_vector=R_vector)
        return filtered_model

    def summary_as_txt(self, template: ExchangeTemplate,
                       dmi_verbose=False, verbose=False):
        r"""
        Return exchange model based on the template file in .txt format.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template of the desired exchange model. 
            (see :ref:`rad-make-template`)
        dmi_verbose : bool, default False
            Whenever to write each individual DMI vector for 
            each average exchange. Ignored if ``verbose`` = True.
        verbose : bool, default False
            Whenever to write everything in a verbose manner.
            (see :ref:`tb2j-extractor_verbose-ref`)

        Returns
        -------
        summary : str
            Exchange information as a string.
        """

        summary = ""
        for name in template.names:
            summary += f"{name}\n"
            if verbose:
                for atom1, atom2, R in template.names[name]:

                    # Get values from the model
                    bond = self.bonds[atom1][atom2][R]
                    if not isinstance(bond, Bond):
                        raise TypeError
                    iso = bond.iso
                    aniso = bond.aniso
                    dmi = bond.dmi
                    abs_dmi = bond.dmi_module
                    rel_dmi = bond.dmi_vs_iso
                    matrix = bond.matrix

                    # Write the values
                    summary += (
                        f"  {atom1:3} {atom2:3} " +
                        f"({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n" +
                        f"    Isotropic: {iso:.4f}\n" +
                        f"    Anisotropic:\n" +
                        f"        {aniso[0][0]:7.4f}  " +
                        f"{aniso[0][1]:7.4f}  " +
                        f"{aniso[0][2]:7.4f}\n" +
                        f"        {aniso[1][0]:7.4f}  " +
                        f"{aniso[1][1]:7.4f}  " +
                        f"{aniso[1][2]:7.4f}\n" +
                        f"        {aniso[2][0]:7.4f}  " +
                        f"{aniso[2][1]:7.4f}  " +
                        f"{aniso[2][2]:7.4f}\n" +
                        f"    DMI: {dmi[0]:.4f} " +
                        f"{dmi[1]:.4f} " +
                        f"{dmi[2]:.4f}\n"
                        f"    |DMI|: {abs_dmi:.4f}\n" +
                        f"    |DMI/J| {rel_dmi:.4f}\n" +
                        f"    Matrix:\n" +
                        f"        {matrix[0][0]:7.4f}  " +
                        f"{matrix[0][1]:7.4f}  " +
                        f"{matrix[0][2]:7.4f}\n" +
                        f"        {matrix[1][0]:7.4f}  " +
                        f"{matrix[1][1]:7.4f}  " +
                        f"{matrix[1][2]:7.4f}\n" +
                        f"        {matrix[2][0]:7.4f}  " +
                        f"{matrix[2][1]:7.4f}  " +
                        f"{matrix[2][2]:7.4f}\n\n")
            else:
                # Compute mean values
                iso = 0
                aniso = np.zeros((3, 3), dtype=float)
                dmi = np.zeros(3, dtype=float)
                abs_dmi = 0
                rel_dmi = 0
                for atom1, atom2, R in template.names[name]:
                    bond = self.bonds[atom1][atom2][R]
                    if not isinstance(bond, Bond):
                        raise TypeError
                    iso += bond.iso
                    aniso += bond.aniso
                    dmi += bond.dmi
                    abs_dmi += bond.dmi_module
                    rel_dmi += bond.dmi_vs_iso
                iso /= len(template.names[name])
                aniso /= len(template.names[name])
                dmi /= len(template.names[name])
                abs_dmi /= len(template.names[name])
                rel_dmi /= len(template.names[name])

                # Write mean values
                summary += (
                    f"    Isotropic: {iso:.4f}\n" +
                    f"    Anisotropic:\n" +
                    f"        {aniso[0][0]:7.4f}  " +
                    f"{aniso[0][1]:7.4f}  " +
                    f"{aniso[0][2]:7.4f}\n" +
                    f"        {aniso[1][0]:7.4f}  " +
                    f"{aniso[1][1]:7.4f}  " +
                    f"{aniso[1][2]:7.4f}\n" +
                    f"        {aniso[2][0]:7.4f}  " +
                    f"{aniso[2][1]:7.4f}  " +
                    f"{aniso[2][2]:7.4f}\n" +
                    f"    |DMI|: {abs_dmi:.4f}\n" +
                    f"    |DMI/J| {abs(rel_dmi):.4f}\n")

                # Write additional info on DMI
                if dmi_verbose:
                    for atom1, atom2, R in template.names[name]:
                        bond = self.bonds[atom1][atom2][R]
                        if not isinstance(bond, Bond):
                            raise TypeError
                        dmi = bond.dmi
                        summary += (
                            f"    DMI: " +
                            f"{dmi[0]:7.4f} " +
                            f"{dmi[1]:7.4f} " +
                            f"{dmi[2]:7.4f} " +
                            f"({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n")
                    summary += "\n"
                # Write only mean value of DMI
                else:
                    summary += (
                        f"    DMI: " +
                        f"{dmi[0]:.4f} " +
                        f"{dmi[1]:.4f} " +
                        f"{dmi[2]:.4f}\n")
            summary += "\n"
        return summary

    def summary_as_py(self, template: ExchangeTemplate):
        r"""
        Return exchange model based on the template file in .py format.

        For the format see :ref:`tb2j-extractor_verbose-ref`.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template of the desired exchange model. 
            (see :ref:`rad-make-template`)

        Returns
        -------
        summary : str
            Exchange information as a python script.
        """

        output_python_iso = "iso = {\n"
        output_python_aniso = "aniso = {\n"
        output_python_dmi = "dmi = {\n"
        output_python_matrix = "matrix = {\n"
        for name in template.names:
            output_python_iso += f"    '{name}':\n" + "    {\n"
            output_python_aniso += f"    '{name}':\n" + "    {\n"
            output_python_dmi += f"    '{name}':\n" + "    {\n"
            output_python_matrix += f"    '{name}':\n" + "    {\n"

            for atom1, atom2, R in template.names[name]:
                bond = self.bonds[atom1][atom2][R]
                if not isinstance(bond, Bond):
                    raise TypeError
                iso = bond.iso
                aniso = bond.aniso
                dmi = bond.dmi
                matrix = bond.matrix
                output_python_iso += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): {iso},\n")
                output_python_aniso += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): np.array([" +
                    f"[{aniso[0][0]}, {aniso[0][1]}, {aniso[0][2]}], " +
                    f"[{aniso[1][0]}, {aniso[1][1]}, {aniso[1][2]}], " +
                    f"[{aniso[2][0]}, {aniso[2][1]}, {aniso[2][2]}]]),\n")
                output_python_dmi += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}):" +
                    f" np.array([{dmi[0]}, {dmi[1]}, {dmi[2]}]),\n")
                output_python_matrix += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): np.array([" +
                    f"[{matrix[0][0]}, {matrix[0][1]}, {matrix[0][2]}], " +
                    f"[{matrix[1][0]}, {matrix[1][1]}, {matrix[1][2]}], " +
                    f"[{matrix[2][0]}, {matrix[2][1]}, {matrix[2][2]}]]),\n")

            output_python_iso += "    },\n"
            output_python_aniso += "    },\n"
            output_python_dmi += "    },\n"
            output_python_matrix += "    },\n"
        output_python_iso += "}\n\n"
        output_python_aniso += "}\n\n"
        output_python_dmi += "}\n\n"
        output_python_matrix += "}\n\n"
        summary = ("import numpy as np\n" +
                   output_python_iso +
                   output_python_aniso +
                   output_python_dmi +
                   output_python_matrix)
        return summary
