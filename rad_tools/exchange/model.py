r"""
Exchange model.

Write a tutorial with docstring here.
"""

from copy import deepcopy
from math import cos, pi, sin, sqrt

import numpy as np

from typing import Tuple


from rad_tools.exchange.bond import Bond
from rad_tools.exchange.template import ExchangeTemplate
from rad_tools.crystal.crystal import Crystal
from rad_tools.crystal.atom import Atom


class ExchangeModel:
    r"""
    Main entry point for the exchange model.

    By default the following notation of the Hamiltonian is assumed:

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}){i,j} \hat{\boldsymbol{S}}_j

    where double counting is present (:math:`ij` and :math:`ji` is in the sum).
    Spin vectors are **not** normalized.

    However, it can be changed via :py:meth:`.change_notation` method and
    checked via :py:attr:`.notation`.

    Parameters
    ----------
    crystal : :py:class:`.Crystal`, optional
        Crystal on which the model is build.
        By default it is orthonormal lattice
        ("CUB with :math:`a = 1`) with no atoms.

    Attributes
    ----------
    crystal : :py:class:`.Crystal`
    bonds : dict
        Dictionary of bonds.

        .. code-block:: python

            {(Atom_1, Atom_2, R): bond, ...}
    magnetic_atoms : list
        List of magnetic atoms.
    cell_list : list
    number_spins_in_unit_cell : int
    space_dimensions : tuple of floats
    """

    _notations = ["double-counting", "minus-sign", "spin-normalised"]

    def __init__(self, crystal: Crystal = None) -> None:
        self._crystal = None
        if crystal is None:
            crystal = Crystal()
        self.crystal = crystal

        self.magnetic_atoms = []
        self.bonds = {}
        self._double_counting = True
        self._spin_normalized = False
        self._factor_one_half = False
        self._factor_two = False
        self._positive_ferromagnetic = True

    def __iter__(self):
        return ExchangeModelIterator(self)

    def __contains__(self, item):
        return item in self.bonds

    def __getitem__(self, key) -> Bond:
        return self.bonds[key]

    @property
    def crystal(self) -> Crystal:
        r"""
        Crystal of the model.

        Crystal, which define the structure.
        See :py:class:`.Crystal`.
        """
        return self._crystal

    @crystal.setter
    def crystal(self, new_crystal: Crystal):
        if not isinstance(new_crystal, Crystal):
            raise TypeError(
                "New crystal is not a Crystal. " + f"Received {type(new_crystal)}."
            )
        self._crystal = new_crystal

    @property
    def notation(self):
        r"""
        Return a string with a simple comment about the Hamiltonian notation.
        """
        result = "H = "
        if self._positive_ferromagnetic:
            result += "-"
        if self._factor_one_half:
            result += "1/2 "
        if self._factor_two:
            result += "2 "
        result += "sum("
        if self._double_counting:
            result += "i,j)"
        else:
            result += "i>=j)"
        result += "S_i J_ij S_j\n"
        if self._double_counting:
            result += "Double counting is present"
        else:
            result += "No double counting"
        if self._spin_normalized:
            result += "Spin vectors are normalized to 1."
        return result

    def set_notation(
        self,
        double_counting=True,
        spin_normalized=False,
        factor_one_half=False,
        factor_two=False,
        positive_ferromagnetic=True,
    ):
        r"""
        Set the notation of the model.

        It is not changing the J values,
        but rather telling how to interpret the J values of the model.

        To **change** the notation of the model use :py:meth:`.change_notation`.
        """

    def change_notation(
        self,
        double_counting=True,
        spin_normalized=False,
        factor_one_half=False,
        factor_two=False,
        positive_ferromagnetic=True,
    ):
        r"""
        Change the notation of the Hamiltonian.

        This method changes J values with respect to the ``new_notation``
        and the notation of the model (use :py:attr:`.notation` to check it).
        In order to tell how to interpret J values (i.e. to set notation)
        use :py:meth:`.set_notation`.

        Parameters
        ----------
        double_counting : boll, default True
        spin_normalized : bool, default False
            Pay attention to the values of atom's spins.
        factor_one_half : bool, default False
        factor_two : bool, default False
        positive_ferromagnetic : bool, True
        """

        multiplier = 1

        if self._double_counting and not double_counting:
            multiplier *= 2
        elif not self._double_counting and double_counting:
            multiplier *= 0.5

        if self._factor_one_half and not factor_one_half:
            multiplier *= 0.5
        elif not self._factor_one_half and factor_one_half:
            multiplier *= 2

        if self._factor_two and not factor_two:
            multiplier *= 2
        elif not self._factor_two and factor_two:
            multiplier *= 0.5

        if self._positive_ferromagnetic ^ positive_ferromagnetic:
            multiplier *= -1

        for atom1, atom2, R in self:
            bond = self.bonds[(atom1, atom2, R)]
            if self._spin_normalized and not spin_normalized:
                factor = multiplier / atom1.spin / atom2.spin
            elif not self._spin_normalized and spin_normalized:
                factor = multiplier * atom1.spin * atom2.spin
            else:
                factor = multiplier
            self.bonds[(atom1, atom2, R)] = factor * bond

    @property
    def cell_list(self):
        r"""
        List of cells from the model.

        Return the list of R for all cell that are present in the model.
        Not ordered.

        Returns
        -------
        cells : (n, 3) :numpy:`ndarray`
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
            x1, y1, z1 = self.crystal.get_atom_coordinates(atom1)
            x2, y2, z2 = self.crystal.get_atom_coordinates(atom2, R)
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

    def add_bond(self, bond: Bond, atom1: Atom, atom2: Atom, R):
        r"""
        Add one bond to the model.

        Parameters
        ----------
        bond : :py:class:`.Bond`
            An instance of :py:class:`Bond` class with the information about
            exchange parameters.
        atom1 : :py:class:`Atom`
            Atom object in (0, 0, 0) unit cell.
        atom2 : :py:class:`Atom`
            Atom object in R unit cell.
        R : tuple of ints
            Vector of the unit cell for atom2.
            In the relative coordinates (i,j,k).
        """

        if atom1 not in self.magnetic_atoms:
            self.add_atom(atom1)
        if atom2 not in self.magnetic_atoms:
            self.add_atom(atom2)

        self.bonds[(atom1, atom2, R)] = bond

    def remove_bond(self, atom1: Atom, atom2: Atom, R):
        r"""
        Remove one bond from the model.

        Parameters
        ----------
        atom1 : py:class:`.Atom`
            Atom object in (0, 0, 0) unit cell.
        atom2 : py:class:`.Atom`
            Atom object in R unit cell.
        R : tuple of ints
            Radius vector of the unit cell for atom2 (i,j,k).
        """

        try:
            del self.bonds[(atom1, atom2, R)]
        except KeyError:
            pass

    def add_atom(self, atom: Atom):
        r"""
        Add atom to the model.

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        """

        self.magnetic_atoms.append(atom)
        self.crystal.add_atom(atom)

    def remove_atom(self, atom):
        r"""
        Remove magnetic atom from the model.

        Note: this method will remove atom itself and all the
        bonds, which starts or ends in this atom, if atom is magnetic.

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        """

        bonds_for_removal = []
        if atom in self.magnetic_atoms:
            self.magnetic_atoms.remove(atom)
        for atom1, atom2, R in self:
            if atom1 == atom or atom2 == atom:
                bonds_for_removal.append((atom1, atom2, R))

        for bond in bonds_for_removal:
            self.remove_bond(*bond)

        self.crystal.remove_atom(atom)

    def filter(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        r"""
        Filter the exchange entries based on the given conditions.

        The result will be defined by logical conjugate of the conditions.
        Saying so the filtering will be performed for each given condition
        one by one.

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

        See Also
        --------
        filtered : Returns new object.

        Notes
        -----
        This method modifies the instance at which it is called.
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
            dis = self.crystal.get_distance(atom1, atom2, R)

            if max_distance is not None and dis > max_distance:
                bonds_for_removal.add((atom1, atom2, R))

            if min_distance is not None and dis < min_distance:
                bonds_for_removal.add((atom1, atom2, R))

            if R_vector is not None and R not in R_vector:
                bonds_for_removal.add((atom1, atom2, R))
            # Here literals, not objects are compared, because in general
            # template only has information about literals and R.
            if (
                template is not None
                and (atom1.literal, atom2.literal, R) not in template
            ):
                bonds_for_removal.add((atom1, atom2, R))

        for atom1, atom2, R in bonds_for_removal:
            self.remove_bond(atom1, atom2, R)

    def remove(self, template):
        r"""
        Remove bonds from the exchange model, based on the template.


        Parameters
        ----------
        template : list or :py:class:`.ExchangeTemplate`.
            List of pairs, which will remain in the model. ::

                [(atom1, atom2, R), ...]

        Notes
        -----
        This method modifies the instance at which it is called.
        """

        if isinstance(template, ExchangeTemplate):
            template = template.get_list()
        template = set(template)

        # Here literals, not objects are compared, because in general
        # template only has information about literals and R.
        for atom1, atom2, R in template:
            self.remove_bond(
                self.crystal.get_atom(atom1), self.crystal.get_atom(atom2), R
            )

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

        See Also
        --------
        filter : Modifies current object.

        Notes
        -----
        This method is not modifying the instance at which it is called.
        It creates a new instance with merged :py:attr:`.bonds` and all the
        other attributes will be copied (through deepcopy).
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

        See Also
        --------
        forced_symmetry : Returns new object.
        """

        if not isinstance(template, ExchangeTemplate):
            raise TypeError

        self.filter(template=template)

        for name in template.names:
            bonds = template.names[name]
            symm_matrix = np.zeros((3, 3), dtype=float)
            dmi_module = 0
            for atom1, atom2, R in bonds:
                # Here literals, not objects are obtained, because in general
                # template only has information about literals and R.
                bond = (self.crystal.get_atom(atom1), self.crystal.get_atom(atom2), R)
                symm_matrix = symm_matrix + self[bond].symm_matrix
                dmi_module = dmi_module + self[bond].dmi_module

            symm_matrix = symm_matrix / len(bonds)
            dmi_module = dmi_module / len(bonds)

            for atom1, atom2, R in bonds:
                bond = (self.crystal.get_atom(atom1), self.crystal.get_atom(atom2), R)
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

        See Also
        --------
        force_symmetry: Modifies current object.
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
            Whether to force the symmetry of the template exchange model.
            If ``False`` then each individual bond is written, otherwise
            exchange parameters of the template model are written.
        isotropic : bool, default False
            Whether to output isotropic exchange.
        anisotropic : bool, default False
            Whether to output anisotropic exchange.
        out_matrix : bool, default False
            Whether to output whole matrix exchange.
        out_dmi : bool, default False
            Whether to output DMI exchange.

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
                atom1 = self.crystal.get_atom(atom1)
                atom2 = self.crystal.get_atom(atom2)
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
        energy = spin_vector @ energy @ spin_vector
        return energy

    # OLD METHODS AND ATTRIBUTES, KEPT FOR BACKWARD COMPATIBILITY

    @property
    def cell(self):
        r"""
        Matrix of lattice vectors.

        See :py:attr:`.Crystal.cell`.
        """

        return self.crystal.lattice.cell

    @cell.setter
    def cell(self, new_cell):
        self.crystal.cell = new_cell

    @property
    def a(self):
        r"""
        First lattice vector :math:`\vec{a} = (a_x, a_y, a_z)`.
        """

        return self.crystal.lattice.a1

    @property
    def b(self):
        r"""
        Second lattice vector :math:`\vec{b} = (b_x, b_y, b_z)`.
        """

        return self.crystal.lattice.a2

    @property
    def c(self):
        r"""
        Third lattice vector :math:`\vec{b} = (b_x, b_y, b_z)`.
        """

        return self.crystal.lattice.a3

    @property
    def len_a(self):
        r"""
        Length of lattice vector :math:`\vec{a}`.
        """

        return self.crystal.lattice.a

    @property
    def len_b(self):
        r"""
        Length of lattice vector :math:`\vec{b}`.
        """

        return self.crystal.lattice.b

    @property
    def len_c(self):
        r"""
        Length of lattice vector :math:`\vec{c}`.
        """

        return self.crystal.lattice.c

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        See :py:attr:`.Lattice.unit_cell_volume`
        """

        return self.crystal.lattice.unit_cell_volume

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector

        See :py:attr:`.Lattice.b1`
        """

        return self.crystal.lattice.b1

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector

        See :py:attr:`.Lattice.b2`
        """

        return self.crystal.lattice.b2

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector

        See :py:attr:`.Lattice.b3`
        """

        return self.crystal.lattice.b3

    def get_atom_coordinates(self, atom, R=(0, 0, 0)):
        r"""
        Getter for the atom coordinates.

        See :py:meth:`.Crystal.get_atom_coordinates`.
        """

        return self.crystal.get_atom_coordinates(atom, R)

    def get_bond_vector(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for distance between the atom1 and atom2.

        See :py:meth:`.Crystal.get_vector`
        """

        return self.crystal.get_vector(atom1, atom2, R)

    def get_distance(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for distance between the atom1 and atom2.

        See :py:meth:`.Crystal.get_distance`
        """

        return self.crystal.get_distance(atom1, atom2, R)


class ExchangeModelIterator:
    def __init__(self, exchange_model: ExchangeModel) -> None:
        self._bonds = list(exchange_model.bonds)
        self._index = 0

    def __next__(self) -> Tuple[Atom, Atom, tuple]:
        if self._index < len(self._bonds):
            result = self._bonds[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self
