r"""
Exchange model.

Write a tutorial with docstrings here.
"""

from copy import deepcopy
from math import cos, pi, sin, sqrt

import numpy as np

from rad_tools.exchange.bond import Bond
from rad_tools.exchange.template import ExchangeTemplate


class ExchangeModel:
    r"""
    Main entry point for the exchange model.

    Attributes
    ----------
    bonds : dict
        Dictionary of bonds.

        .. code-block:: python

            {(Atom_1, Atom_2, R): bond, ...}
    magnetic_atoms : dict
        Dictionary with keys : str - marks of atoms and value : 1 x 3 array
        - coordinate of the atom in Angstroms.
    nonmagnetic_atoms : dist
        Dictionary with keys : str - marks of atoms and value : 1 x 3 array
        - coordinate of the atom in Angstroms. Only non-magnetic atoms.
    cell : 3 x 3 array of floats
    a : 3 x 1 array
    b : 3 x 1 array
    c : 3 x 1 array
    len_a : float
    len_b : float
    len_c : float
    unit_cell_volume : float
    b1 : 3 x 1 array
    b2 : 3 x 1 array
    b3 : 3 x 1 array
    cell_list : list
    number_spins_in_unit_cell : int

    space_dimensions : tuple of floats
    """

    def __init__(self) -> None:
        self._cell = np.identity(3, dtype=float)
        self.magnetic_atoms = {}
        self.nonmagnetic_atoms = {}
        self.bonds = {}

    def __iter__(self):
        return ExchangeModelIterator(self)

    def __contains__(self, item):
        return item in self.bonds

    def __getitem__(self, key) -> Bond:
        return self.bonds[key]

    @property
    def cell(self):
        r"""
        Matrix of lattice vectors in Angstrom.

        .. code-block:: python

            [[a_x  a_y  a_z],
             [b_x  b_y  b_z],
             [c_x  x_y  c_z]]

        Default cell is orthonormal:

        .. code-block:: python

            np.identity(3, dtype=float)
        """

        return self._cell

    @cell.setter
    def cell(self, new_cell):
        new_cell = np.array(new_cell)
        if new_cell.shape != (3, 3):
            raise ValueError("Cell matrix shape have to be equal (3, 3)")
        if np.linalg.det(new_cell) == 0:
            raise ValueError("Lattice vectors are coplanar.")
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

        return sqrt(np.sum(self.a**2))

    @property
    def len_b(self):
        r"""
        Length of lattice vector :math:`\vec{b}`.
        """

        return sqrt(np.sum(self.b**2))

    @property
    def len_c(self):
        r"""
        Length of lattice vector :math:`\vec{c}`.
        """

        return sqrt(np.sum(self.c**2))

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        .. math::

            V = \delta_{\vec{A}}\cdot(\delta_{\vec{B}}\times\delta_{\vec{C}})
        """

        return np.dot(self.a, np.cross(self.b, self.c))

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{b}\times\vec{c}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.b, self.c)

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{c}\times\vec{a}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.c, self.a)

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}\times\vec{b}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return 2 * pi / self.unit_cell_volume * np.cross(self.a, self.b)

    @property
    def cell_list(self):
        r"""
        List of cells from the model.

        Return the list of R for all cell that are present in the model.
        Not ordered.

        Returns
        -------
        cells : n x 3 array
            Array of n unit cells.
        """

        cells = set()
        for atom1, atom2, R in self:
            cells.add(R)
        return np.array(list(cells))

    @property
    def number_spins_in_unit_cell(self):
        r"""
        Number of spins (or magnetic atoms) in the unit cell.
        """

        return len(self.magnetic_atoms)

    @property
    def space_dimensions(self):
        r"""
        Model minimum and maximum coordinates in real space.

        Takes into account only an atom which has at least
        one bond associated with it.

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

        x_max, y_max, z_max = None, None, None
        x_min, y_min, z_min = None, None, None
        for atom1, atom2, R in self:
            x1, y1, z1 = self.get_atom_coordinates(atom1)
            x2, y2, z2 = self.get_atom_coordinates(atom2, R)
            if x_max is None:
                x_min, y_min, z_min = x1, y1, z1
                x_max, y_max, z_max = x1, y1, z1
            x_max = max(x1, x2, x_max)
            y_max = max(y1, y2, y_max)
            z_max = max(z1, z2, z_max)
            x_min = min(x1, x2, x_min)
            y_min = min(y1, y2, y_min)
            z_min = min(z1, z2, z_min)
        return x_min, y_min, z_min, x_max, y_max, z_max

    def round(self, decimals=4):
        r"""Round exchange parameters.

        Parameters
        ----------
        decimals : int, default 4
            Number of decimals after the comma.
        """

        for atom1, atom2, R in self:
            self[(atom1, atom2, R)].round(decimals)

    def add_bond(self, bond: Bond, atom1, atom2, R):
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
            Vector of the unit cell for atom2.
            In the relative coordinates in real space.
        """

        self.bonds[(atom1, atom2, R)] = bond

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

        try:
            del self.bonds[(atom1, atom2, R)]
        except KeyError:
            pass

    def add_atom(self, name, a, b, c):
        r"""
        Add magnetic atom to the model.

        Parameters
        ----------
        name : str
            Mark for the atom. Note: if an atom with the same mark already
            exists in ``magnetic_atoms`` then it will be rewritten.
        a : int or float
            Relative coordinate of the atom along a lattice vector.
        b : int or float
            Relative coordinate of the atom along b lattice vector.
        c : int or float
            Relative coordinate of the atom along c lattice vector.
        """

        self.magnetic_atoms[name] = np.array([a, b, c])

    def remove_atom(self, name):
        r"""
        Remove magnetic atom from the model.

        Note: this method will remove atom from ``magnetic_atoms`` and all the
        bonds, which starts or ends in this atom.

        Parameters
        ----------
        name : str
            Mark for the atom.
        """

        bonds_for_removal = []
        if name in self.magnetic_atoms:
            del self.magnetic_atoms[name]
        for atom1, atom2, R in self:
            if atom1 == name or atom2 == name:
                bonds_for_removal.append((atom1, atom2, R))

        for bond in bonds_for_removal:
            self.remove_bond(*bond)

    def get_atom_coordinates(self, atom, R=(0, 0, 0)):
        r"""
        Getter for the atom coordinates.

        Parameters
        ----------
        atom : str
            Name of atom1 in R unit cell.
        R : 1 x 3 array, default (0, 0, 0)
            Radius vector of the unit cell for atom2.

        Returns
        -------
        coordinates : 1 x 3 array
            Coordinates of atom in the cell R in real space.
        """

        R = np.array(R)
        try:
            return np.matmul(R + self.magnetic_atoms[atom], self.cell)
        except KeyError:
            return np.matmul(R + self.nonmagnetic_atoms[atom], self.cell)

    def get_bond_centre_coordinates(self, atom1, atom2, R=(0, 0, 0)):
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
        coordinates : 1 x 3 array
            Coordinates of the middle of a bond in real space.
        """

        atom1 = self.get_atom_coordinates(atom1)
        atom2 = self.get_atom_coordinates(atom2, R)

        return (atom1 + atom2) / 2

    def get_bond_vector(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for vector between the atom1 and atom2.

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
        v : 1 x 3 array
            Vector from atom1 to atom2 of the bond.
        """

        atom1 = self.get_atom_coordinates(atom1)
        atom2 = self.get_atom_coordinates(atom2, R)

        return atom2 - atom1

    def get_distance(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for distance between the atom1 and atom2.

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
        distance : floats
            Distance between atom1 and atom2.
        """

        return sqrt(np.sum(self.get_bond_vector(atom1, atom2, R) ** 2))

    def filter(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        r"""
        Filter the exchange entries based on the given conditions.

        The result will be defined by logical conjugate of the conditions.
        Saying so the filtering will be performed for each given condition
        one by one.

        .. note::
            This method modifies the instance at which it is called.

        Parameters
        ----------
        max_distance : float or int, optional
            Distance for sorting, the condition is <<less or equal>>.
        min_distance : float or int, optional
            Distance for sorting, the condition is <<more or equal>>.
        template : list or :py:class:`.ExchangeTemplate`.
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
            else:
                raise TypeError("Type is not supported, supported: list.")

        bonds_for_removal = set()
        for atom1, atom2, R in self:
            dis = self.get_distance(atom1, atom2, R)

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

    def remove(self, template):
        r"""
        Remove bonds from the exchange model, based on the template.


        .. note::
            This method modifies the instance at which it is called.

        Parameters
        ----------
        template : list or :py:class:`.ExchangeTemplate`.
            List of pairs, which will remain in the model. ::

                [(atom1, atom2, R), ...]
        """

        if isinstance(template, ExchangeTemplate):
            template = template.get_list()
        template = set(template)

        for atom1, atom2, R in template:
            self.remove_bond(atom1, atom2, R)

    def filtered(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        r"""
        Create filtered exchange model based on the given conditions.

        The result will be defined by logical conjugate of the conditions.
        Saying so the filtering will be performed for each given condition
        one by one.
        Note: this method is not modifying the instance at which it is called.
        It will create a new instance with sorted :py:attr:`bonds` and all
        the other attributes will be copied (through :py:func:`deepcopy`).

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
        filtered_model.filter(
            max_distance=max_distance,
            min_distance=min_distance,
            template=template,
            R_vector=R_vector,
        )
        return filtered_model

    def force_symmetry(self, template):
        r"""
        Force the model to have the symmetries of the template.

        Takes mean values of the parameters.
        For DMI interaction directions are kept the same,
        but the length of the DMI vector is scaled.

        This method modifies an instance on which it was called.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template.
        """

        if not isinstance(template, ExchangeTemplate):
            raise TypeError

        self.filter(template=template)

        for name in template.names:
            bonds = template.names[name]
            symm_matrix = np.zeros((3, 3), dtype=float)
            dmi_module = 0
            for bond in bonds:
                symm_matrix = symm_matrix + self[bond].symm_matrix
                dmi_module = dmi_module + self[bond].dmi_module

            symm_matrix = symm_matrix / len(bonds)
            dmi_module = dmi_module / len(bonds)

            for bond in bonds:
                if self[bond].dmi_module != 0:
                    asymm_factor = dmi_module / self[bond].dmi_module
                else:
                    asymm_factor = 0
                self[bond].matrix = symm_matrix + self[bond].asymm_matrix * asymm_factor

    def forced_symmetry(self, template):
        r"""
        Force the model to have the symmetries of the template.

        Takes mean values of the parameters.
        Respect the direction of the DMI vectors.

        This method return a new instance.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template.

        Returns
        -------
        new_model : :py:class:`.ExchangeModel`
            Exchange model with forced symmetry.
        """

        new_model = deepcopy(self)
        new_model.force_symmetry(template=template)
        return new_model

    def summary_as_txt(
        self,
        template: ExchangeTemplate,
        decimals=4,
        force_symmetry=False,
        isotropic=False,
        anisotropic=False,
        out_matrix=False,
        out_dmi=False,
    ):
        r"""
        Return exchange model based on the template file in .txt format.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template of the desired exchange model.
            (see :ref:`rad-make-template`)
        accuracy : int, default 4
            Accuracy for the exchange values
        force_symmetry: bool, default False
            Whenever to force the symmetry of the template exchange model.
            If ``False`` then each individual bond is written, otherwise
            exchange parameters of the template model are written.
        isotropic : bool, default False
            Whenever to output isotropic exchange.
        anisotropic : bool, default False
            Whenever to output anisotropic exchange.
        out_matrix : bool, default False
            Whenever to output whole matrix exchange.
        out_dmi : bool, default False
            Whenever to output DMI exchange.

        Returns
        -------
        summary : str
            Exchange information as a string.
        """

        if force_symmetry:
            self.force_symmetry(template=template)
        else:
            self.filter(template=template)
        self.round(decimals=decimals)
        summary = ""
        for name in template.names:
            scalar_written = False
            for atom1, atom2, R in template.names[name]:
                # Get values from the model
                bond = self[(atom1, atom2, R)]
                iso = bond.iso
                aniso = bond.aniso
                dmi = bond.dmi
                abs_dmi = bond.dmi_module
                rel_dmi = abs_dmi / iso
                matrix = bond.matrix

                if not force_symmetry:
                    # Write the values
                    summary += (
                        f"{atom1:3} {atom2:3} "
                        + f"({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n"
                    )
                    if isotropic:
                        summary += f"    Isotropic: {iso:.{decimals}f}\n"
                    if anisotropic:
                        summary += (
                            f"    Anisotropic:\n"
                            + f"        {aniso[0][0]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[0][1]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[0][2]:{decimals+3}.{decimals}f}\n"
                            + f"        {aniso[1][0]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[1][1]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[1][2]:{decimals+3}.{decimals}f}\n"
                            + f"        {aniso[2][0]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[2][1]:{decimals+3}.{decimals}f}  "
                            + f"{aniso[2][2]:{decimals+3}.{decimals}f}\n"
                        )
                    if out_matrix:
                        summary += (
                            f"    Matrix:\n"
                            + f"        {matrix[0][0]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[0][1]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[0][2]:{decimals+3}.{decimals}f}\n"
                            + f"        {matrix[1][0]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[1][1]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[1][2]:{decimals+3}.{decimals}f}\n"
                            + f"        {matrix[2][0]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[2][1]:{decimals+3}.{decimals}f}  "
                            + f"{matrix[2][2]:{decimals+3}.{decimals}f}\n"
                        )
                    if out_dmi:
                        summary += (
                            f"    |DMI|: {abs_dmi:.{decimals}f}\n"
                            + f"    |DMI/J|: {rel_dmi:.{decimals}f}\n"
                            + f"    DMI: {dmi[0]:.{decimals}f} "
                            + f"{dmi[1]:.{decimals}f} "
                            + f"{dmi[2]:.{decimals}f}\n"
                        )
                    summary += "\n"
                else:
                    if not scalar_written:
                        scalar_written = True
                        summary += f"{name}\n"
                        if isotropic:
                            summary += f"    Isotropic: {iso:.{decimals}f}\n"
                        if anisotropic:
                            summary += (
                                f"    Anisotropic:\n"
                                + f"        {aniso[0][0]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[0][1]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[0][2]:{decimals+3}.{decimals}f}\n"
                                + f"        {aniso[1][0]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[1][1]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[1][2]:{decimals+3}.{decimals}f}\n"
                                + f"        {aniso[2][0]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[2][1]:{decimals+3}.{decimals}f}  "
                                + f"{aniso[2][2]:{decimals+3}.{decimals}f}\n"
                            )
                        if out_matrix:
                            summary += (
                                f"    Matrix:\n"
                                + f"        {matrix[0][0]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[0][1]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[0][2]:{decimals+3}.{decimals}f}\n"
                                + f"        {matrix[1][0]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[1][1]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[1][2]:{decimals+3}.{decimals}f}\n"
                                + f"        {matrix[2][0]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[2][1]:{decimals+3}.{decimals}f}  "
                                + f"{matrix[2][2]:{decimals+3}.{decimals}f}\n"
                            )
                        if out_dmi:
                            summary += (
                                f"    |DMI|: {abs_dmi:.{decimals}f}\n"
                                + f"    |DMI/J|: {rel_dmi:.{decimals}f}\n"
                            )
                    if out_dmi:
                        summary += (
                            f"    DMI: {dmi[0]:{decimals+3}.{decimals}f} "
                            + f"{dmi[1]:{decimals+3}.{decimals}f} "
                            + f"{dmi[2]:{decimals+3}.{decimals}f} "
                            + f"({atom1:3} {atom2:3} "
                            + f"{R[0]:2.0f} {R[1]:2.0f} {R[2]:2.0f})\n"
                        )
            if force_symmetry:
                summary += "\n"

        return summary

    def ferromagnetic_energy(self, theta=0, phi=0):
        r"""
        Compute energy of the model assuming ferromagnetic state.

        With Hamiltonian of the form:

        .. math::

            \hat{H} = - \sum_{i,j} \vec{S}_i J_{ij} \vec{S}_j

        where spin vectors are normalized to 1 and :math:`J_{ij}` is the
        exchange matrix.

        Parameters
        ----------
        theta : float, default 0
            Angle between z axis an direction of the magnetization.
            :math:`0 < \theta < 180`
        phi : float, default 0
            angle between x axis an projection of direction of the
            magnetization on xy plane.
            :math:`0 < \phi < 360`

        Returns
        -------
        energy : float
            Energy of ferromagnetic model with magnetization direction defined
            by ``theta`` and ``phi``. In the units of J values.
        """

        theta = theta / 180 * pi
        phi = phi / 180 * pi
        energy = np.zeros((3, 3), dtype=float)
        for bond in self:
            energy -= self[bond].matrix
        spin_vector = np.array(
            [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)]
        )
        energy = np.matmul(np.matmul(spin_vector, energy), spin_vector)
        return energy


class ExchangeModelIterator:
    def __init__(self, exchange_model: ExchangeModel) -> None:
        self._bonds = list(exchange_model.bonds)
        self._index = 0

    def __next__(self) -> Bond:
        if self._index < len(self._bonds):
            result = self._bonds[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self
