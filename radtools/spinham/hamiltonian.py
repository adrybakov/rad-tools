r"""
Spin Hamiltonian module.

Write a tutorial with docstring here.
"""

__all__ = ["SpinHamiltonian", "ExchangeHamiltonian", "dump_spinham_txt"]

from copy import deepcopy
from typing import Iterable, Tuple

import numpy as np

from radtools.constants import TORADIANS
from radtools.crystal.atom import Atom
from radtools.crystal.crystal import Crystal
from radtools.decorate.array import print_2d_array
from radtools.exceptions import NotationError
from radtools.spinham.constants import PREDEFINED_NOTATIONS
from radtools.spinham.parameter import ExchangeParameter
from radtools.spinham.template import ExchangeTemplate
from radtools.decorate.array import print_2d_array
from radtools.spinham.constants import TXT_FLAGS
from radtools.decorate.stats import logo


class SpinHamiltonian(Crystal):
    r"""
    Spin Hamiltonian.

    By default the notation of the spin Hamiltonian is not defined
    and could be different in different context.
    However, in could always be checked via :py:attr:`.notation`.
    In user-specific cases it is the responsibility of the user to
    set the interpretation of the Hamiltonian`s notation.

    For the predefined notations see :py:meth:`.notation`.

    Child of the :py:class:`.Crystal` class.

    .. versionchanged:: 0.8.0 Renamed from SpinHamiltonian to SpinHamiltonian.

    Notes
    -----
    You can still use the old name, but it is deprecated and will be removed in the future.

    Parameters
    ----------
    crystal : :py:class:`.Crystal`, optional
        Crystal on which the spin Hamiltonian is build.
        By default it is cubic (:math:`a=1`) lattice with no atoms.
    notation : str or tuple of two bool and one float, optional
        One of the predefined notations or tuple with custom notation.
        See :py:attr:`.notation` for details.
    **kwargs
        Keyword arguments for the :py:class:`.Crystal` constructor.

    """

    def __init__(self, crystal: Crystal = None, notation=None, **kwargs) -> None:
        if crystal is not None:
            kwargs["atoms"] = crystal.atoms
            kwargs["cell"] = crystal.cell

        super().__init__(**kwargs)

        self._bonds = {}

        # Notation settings
        self._double_counting = None
        self._spin_normalized = None
        self._factor = None
        if notation is not None:
            self.notation = notation

    def __len__(self):
        return self._bonds.__len__()

    # Notation attributes
    @property
    def notation(self):
        r"""
        Tuple of the notation.

        It can be set with a

        * string
            One of the predefined notations: "standard", "TB2J", "SpinW"

        * iterable with 3 elements
            First two elements are converted to ``bool``,
            third element is interpreted as a float.
            Order: (double counting, spin normalized, factor).

        Predefined notations:

        * Standard
            (True, False, -1)

            .. math::
                H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are **not** normalized.
        * |TB2J|_
            (True, True, -1)

            .. math::
                H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are normalized to 1.
        * SpinW
            (True, False, 1)

            .. math::
                H = \sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are **not** normalized.


        Setting of this attribute either change the notation
        (if corresponding attributes are defined)
        or set the interpretation (if corresponding attribute are not defined).

        Returns
        -------
        notation : (3,) tuple of two bool and one float
            Values of notation properties:

            .. code-block:: python

                (double_counting, spin_normalized, factor)

        Examples
        --------
        Setting one of the predefined notations:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.SpinHamiltonian()
            >>> model.notation = "standard"
            >>> print(model.notation_string)
            H = - \sum_{i,j} S_i J_{ij} S_j
            Double counting is present.
            Spin vectors are not normalized.
            >>> model.notation
            (True, False, -1.0)
            >>> model.notation = "TB2J"
            >>> print(model.notation_string)
            H = - \sum_{i,j} e_i J_{ij} e_j
            Double counting is present.
            Spin vectors are normalized to 1.
            >>> model.notation
            (True, True, -1.0)
            >>> model.notation = "SpinW"
            >>> print(model.notation_string)
            H = \sum_{i,j} S_i J_{ij} S_j
            Double counting is present.
            Spin vectors are not normalized.
            >>> model.notation
            (True, False, 1.0)

        Setting the notation:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.SpinHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(Cr, Cr, (1, 0, 0), iso=1)
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # For the first time interpretation is set,
            >>> # values of exchange are not changed
            >>> model.notation = "standard"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # Once the notation is set the values
            >>> # are changing if the notation is changed again.
            >>> model.notation = "TB2J"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            2.25
            >>> model.notation = "SpinW"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0
            >>> model.notation = "standard"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0

        Setting individual properties:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.SpinHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(Cr, Cr, (1, 0, 0), iso=1)
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # For the first time interpretation is set,
            >>> # values of exchange are not changed
            >>> model.factor = -1
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # Once the property is set the values
            >>> # are changing if the property is changed again.
            >>> model.factor = 1
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0

        Changing individual properties:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.SpinHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(Cr, Cr, (1, 0, 0), iso=1)
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> model.notation = "standard"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> model.double_counting, model.spin_normalized, model.factor
            (True, False, -1.0)
            >>> model.factor = 1
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0
            >>> model.factor
            1.0
            >>> model.factor = 1/2
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -2.0
            >>> model.factor
            0.5
            >>> model.factor = 2
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -0.5
            >>> model.spin_normalized = True
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.125
            >>> model.double_counting = False
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -2.25


        See Also
        --------
        double_counting
        spin_normalized
        factor
        set_interpretation
        """

        return (self.double_counting, self.spin_normalized, self.factor)

    @notation.setter
    def notation(self, new_notation):
        # Set the notation from predefined notations
        if isinstance(new_notation, str):
            new_notation = new_notation.lower()
            if new_notation not in PREDEFINED_NOTATIONS:
                raise ValueError(
                    f"Predefine notations are: "
                    + f"{list(PREDEFINED_NOTATIONS)}, got: {new_notation}"
                )
            new_notation = PREDEFINED_NOTATIONS[new_notation]
        # Set the notation from five values, converted to bool
        elif isinstance(new_notation, Iterable) and len(new_notation) == 3:
            new_notation = (
                bool(new_notation[0]),
                bool(new_notation[1]),
                float(new_notation[2]),
            )
        else:
            raise ValueError(
                "New notation has to be either a string "
                + "or an iterable with three elements, "
                + f"got: {new_notation}"
            )
        (self.double_counting, self.spin_normalized, self.factor) = new_notation

    @property
    def notation_string(self):
        r"""
        Human-readable representation of the notation.

        Returns
        -------
        notation_string : str
            Human-readable representation of the notation.
        """

        text_result = "H = "

        if self.factor != 1:
            if self.factor == -1:
                text_result += "- "
            elif self.factor == int(self.factor):
                text_result += f"{self.factor:.0f} "
            elif self.factor == 0.5:
                text_result += "1/2 "
            elif self.factor == -0.5:
                text_result += "-1/2 "
            else:
                text_result += f"{self.factor} "

        text_result += "\sum_{"
        if self.double_counting:
            text_result += "i,j}"
        else:
            text_result += "i>=j}"

        if self.spin_normalized:
            text_result += " e_i J_{ij} e_j\n"
        else:
            text_result += " S_i J_{ij} S_j\n"

        if self.double_counting:
            text_result += "Double counting is present.\n"
        else:
            text_result += "No double counting.\n"

        if self.spin_normalized:
            text_result += "Spin vectors are normalized to 1."
        else:
            text_result += "Spin vectors are not normalized."

        return text_result

    def set_interpretation(
        self, double_counting: bool = None, spin_normalized: bool = None, factor=None
    ):
        r"""
        Set the interpretation of the Hamiltonian`s notation.

        It is not changing the J values,
        but rather indicate how to interpret the J values of the Hamiltonian.

        Parameters
        ----------
        double_counting : bool, optional
        spin_normalized : bool, optional
        factor : float or int, optional

        See Also
        --------
        double_counting
        spin_normalized
        factor
        notation
        """

        if double_counting is not None:
            self._double_counting = bool(double_counting)
        if spin_normalized is not None:
            self._spin_normalized = bool(spin_normalized)
        if factor is not None:
            self._factor = factor

    @property
    def double_counting(self) -> bool:
        r"""
        Whether double counting is present in the Hamiltonian.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        double_counting : bool
            ``True`` if double counting is present, ``False`` otherwise.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._double_counting is None:
            raise NotationError("double_counting")
        return self._double_counting

    def _ensure_double_counting(self):
        r"""
        Ensure double counting.

        Check and fix if needed, that for each (atom1, atom2, R)
        (atom2, atom1, -R) is present in the Hamiltonian.
        """

        bonds = list(self)
        for atom1, atom2, (i, j, k), J in bonds:
            if (atom2, atom1, (-i, -j, -k)) not in self:
                self.add_bond(atom2, atom1, (-i, -j, -k), J=J.T)

    def _ensure_no_double_counting(self):
        r"""
        Ensure that there is no double counting in the model.

        If the bond is (atom1, atom2, (i,j,k)), then the bond is kept in the model
        if one of the conditions is true:

        * i > 0
        * i = 0 and j > 0
        * i = 0, j = 0 and k > 0
        * i = 0, j = 0, k = 0 and atom1.index <= atom2.index
        """

        bonds = list(self)

        for atom1, atom2, (i, j, k), J in bonds:
            if (
                i > 0
                or (i == 0 and j > 0)
                or (i == 0 and j == 0 and k > 0)
                or (i == 0 and j == 0 and k == 0 and atom1.index <= atom2.index)
            ):
                if (atom2, atom1, (-i, -j, -k)) in self:
                    self.remove_bond(atom2, atom1, (-i, -j, -k))
            else:
                if (atom2, atom1, (-i, -j, -k)) not in self:
                    self.add_bond(atom2, atom1, (-i, -j, -k), J=J.T)
                    self.remove_bond(atom1, atom2, (i, j, k))

    @double_counting.setter
    def double_counting(self, new_value: bool):
        # Need to be at the beginning, otherwise _ensure methods do not work
        old_value = self._double_counting
        self._double_counting = bool(new_value)

        if old_value is not None:
            factor = 1
            if old_value and not new_value:
                factor = 2
            elif not old_value and new_value:
                factor = 0.5
            if factor != 1:
                if new_value:
                    self._ensure_double_counting()
                else:
                    self._ensure_no_double_counting()
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * factor
        else:
            if new_value:
                self._ensure_double_counting()
            else:
                self._ensure_no_double_counting()

    @property
    def spin_normalized(self) -> bool:
        r"""
        Whether spin is normalized.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        spin_normalized : bool
            ``True`` if spin is normalized, ``False`` otherwise.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._spin_normalized is None:
            raise NotationError("spin_normalized")
        return self._spin_normalized

    @spin_normalized.setter
    def spin_normalized(self, new_value: bool):
        if self._spin_normalized is not None:
            if self._spin_normalized and not new_value:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J / atom1.spin / atom2.spin
            elif not self._spin_normalized and new_value:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * atom1.spin * atom2.spin
        self._spin_normalized = bool(new_value)

    @property
    def factor(self) -> float:
        r"""
        Whether any factor is present in the Hamiltonian.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        factor : float
            value of the factor before the sum in the Hamiltonian.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._factor is None:
            raise NotationError("factor")
        return self._factor

    @factor.setter
    def factor(self, new_factor: bool):
        if self._factor is not None:
            factor = self._factor / new_factor

            if factor != 1:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * factor

        self._factor = float(new_factor)

    def __iter__(self):
        return SpinHamiltonianIterator(self)

    def __contains__(self, key):
        atom1, atom2, R = key
        # If atom is a string, get the atom object
        if isinstance(atom1, str):
            atom1 = self.get_atom(atom1)
        if isinstance(atom2, str):
            atom2 = self.get_atom(atom2)

        key = (atom1, atom2, R)
        return key in self._bonds

    def __getitem__(self, key) -> ExchangeParameter:
        atom1, atom2, R = key
        # If atom is a string, get the atom object
        if isinstance(atom1, str):
            atom1 = self.get_atom(atom1)
        if isinstance(atom2, str):
            atom2 = self.get_atom(atom2)

        key = (atom1, atom2, R)
        return self._bonds[key]

    def __getattr__(self, name):
        # Fix copy/deepcopy RecursionError
        if name in ["__setstate__"]:
            raise AttributeError(name)
        raise AttributeError(name)

    @property
    def crystal(self) -> Crystal:
        r"""
        Crystal of the Hamiltonian.

        Return an independent instance of a crystal.
        You can use it to play with the model`s crystal independently,
        but it will not affect the model itself.

        Returns
        -------
        crystal : :py:class:`.Crystal`
            Crystal of the Hamiltonian.

        See Also
        --------
        Crystal
        """
        return Crystal(self.lattice, self.atoms)

    @property
    def cell_list(self):
        r"""
        List of cells from the Hamiltonian.

        Return the list of R for all cell that are present in the Hamiltonian.
        Not ordered.

        Returns
        -------
        cells : (n, 3) :numpy:`ndarray`
            Array of n unit cells.
        """

        cells = set()
        for atom1, atom2, R, J in self:
            cells.add(R)
        return np.array(list(cells))

    @property
    def magnetic_atoms(self):
        r"""
        Magnetic atoms of the model.

        Atoms with at least bond starting or finishing in it.

        Atoms are ordered with respect to the :py:attr:`.Atom.index`.

        Returns
        -------
        magnetic_atoms : list of :py:class:`.Atom`
            List of magnetic atoms.
        """
        result = set()
        for atom1, atom2, R, J in self:
            result.add(atom1)
            result.add(atom2)

        return sorted(list(result), key=lambda x: x.index)

    @property
    def number_spins_in_unit_cell(self):
        r"""
        Number of spins (or magnetic atoms) in the unit cell.

        Returns
        -------
        number_spins_in_unit_cell : int
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
        for atom1, atom2, R in self._bonds:
            x1, y1, z1 = self.get_atom_coordinates(atom1, relative=False)
            x2, y2, z2 = self.get_atom_coordinates(atom2, R, relative=False)
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

    def __setitem__(self, key, value):
        self.add_bond(*key, value)

    def add_bond(
        self, atom1: Atom, atom2: Atom, R, J: ExchangeParameter = None, **kwargs
    ):
        r"""
        Add one bond to the Hamiltonian.

        More straightforward syntax is advised:

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.SpinHamiltonian(rad.Crystal())
            >>> model[Cr, Cr, (1,0,0)] = J
            >>> (Cr, Cr, (1,0,0)) in model
            True

        It is equivalent to

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.SpinHamiltonian(rad.Crystal())
            >>> model.add_bond(Cr, Cr, (1,0,0), J=J)
            >>> (Cr, Cr, (1,0,0)) in model
            True

        Parameters
        ----------
        atom1 : :py:class:`Atom` or str
            Atom object in (0, 0, 0) unit cell.
            str works only if atom is already in the Hamiltonian.
        atom2 : :py:class:`Atom` or str
            Atom object in R unit cell.
            str works only if atom is already in the Hamiltonian.
        R : tuple of ints
            Vector of the unit cell for atom2.
            In the relative coordinates (i,j,k).
        J : :py:class:`.ExchangeParameter`, optional
            An instance of :py:class:`ExchangeParameter`.
        ** kwargs
            Keyword arguments for the constructor of :py:class:`ExchangeParameter`.
            Ignored if J is given.
        """

        if isinstance(atom1, str):
            atom1 = self.get_atom(atom1)
        elif atom1 not in self.atoms:
            self.add_atom(atom1)

        if isinstance(atom2, str):
            atom2 = self.get_atom(atom2)
        elif atom2 not in self.atoms:
            self.add_atom(atom2)

        if J is None:
            J = ExchangeParameter(**kwargs)

        self._bonds[(atom1, atom2, R)] = J

        # Check for double counting
        double_counting = False
        try:
            double_counting = self.double_counting
        except NotationError:
            pass
        i, j, k = R
        if double_counting and (atom2, atom1, (-i, -j, -k)) not in self:
            self._bonds[(atom2, atom1, (-i, -j, -k))] = J.T

    def __delitem__(self, key):
        self.remove_bond(*key)

    def remove_bond(self, atom1: Atom, atom2: Atom, R):
        r"""
        Remove one bond from the Hamiltonian.

        More straightforward syntax is advised:

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.SpinHamiltonian(rad.Crystal())
            >>> model[Cr, Cr, (1,0,0)] = J
            >>> (Cr, Cr, (1,0,0)) in model
            True
            >>> del model[Cr, Cr, (1,0,0)]
            >>> (Cr, Cr, (1,0,0)) in model
            False


        It is the same as

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.SpinHamiltonian(rad.Crystal())
            >>> model[Cr, Cr, (1,0,0)] = J
            >>> (Cr, Cr, (1,0,0)) in model
            True
            >>> model.remove_bond(Cr, Cr, (1,0,0))
            >>> (Cr, Cr, (1,0,0)) in model
            False

        Parameters
        ----------
        atom1 : py:class:`.Atom`
            Atom object in (0, 0, 0) unit cell.
        atom2 : py:class:`.Atom`
            Atom object in R unit cell.
        R : tuple of ints
            Radius vector of the unit cell for atom2 (i,j,k).
        """

        # If atom is a string, get the atom object
        if isinstance(atom1, str):
            atom1 = self.get_atom(atom1)
        if isinstance(atom2, str):
            atom2 = self.get_atom(atom2)

        try:
            del self._bonds[(atom1, atom2, R)]
        except KeyError:
            raise KeyError(
                f"Bond ({atom2.fullname}, {atom2.fullname}, {R}) is not present in the model."
            )

        # Check for double counting
        double_counting = False
        try:
            double_counting = self.double_counting
        except NotationError:
            pass
        i, j, k = R
        if double_counting and (atom2, atom1, (-i, -j, -k)) in self:
            del self._bonds[(atom2, atom1, (-i, -j, -k))]

    def remove_atom(self, atom):
        r"""
        Remove magnetic atom from the Hamiltonian.

        Note: this method will remove atom itself and all the
        bonds, which starts or ends in this atom, if atom is magnetic.

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        """

        # If atom given as a string, get the atom object
        if isinstance(atom, str):
            atom = self.get_atom(atom)

        bonds_for_removal = []
        for atom1, atom2, R, J in self:
            if atom1 == atom or atom2 == atom:
                bonds_for_removal.append((atom1, atom2, R))

        for bond in bonds_for_removal:
            del self[bond]

        super().remove_atom(atom)

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
            List of pairs, which will remain in the Hamiltonian. ::

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
        for atom1, atom2, R, J in self:
            i, j, k = R
            dis = self.get_distance(atom1, atom2, R)

            if max_distance is not None and dis > max_distance:
                bonds_for_removal.add((atom1, atom2, R))

            if min_distance is not None and dis < min_distance:
                bonds_for_removal.add((atom1, atom2, R))

            # This behaviour should depend on the notation of the Hamiltonian
            condition = R_vector is not None and R not in R_vector
            try:
                if self.double_counting:
                    condition = condition and (-i, -j, -k) not in R_vector
            except NotationError:
                pass
            if condition:
                bonds_for_removal.add((atom1, atom2, R))

            # Here names, not objects are compared, because in general
            # template only has information about names (or fullnames) and R.
            if (
                template is not None
                and (atom1.name, atom2.name, R) not in template
                and (atom1.fullname, atom2.fullname, R) not in template
            ):
                bonds_for_removal.add((atom1, atom2, R))

        for atom1, atom2, R in bonds_for_removal:
            try:
                del self[atom1, atom2, R]
            except KeyError:
                pass

    def filtered(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        r"""
        Create filtered spin Hamiltonian based on the given conditions.

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
            List of pairs, which will remain in the Hamiltonian. ::

                [(atom1, atom2, R), ...]

        R_vector : tuple of ints or list of tuples of ints
            Tuple of 3 integers or list of tuples, specifying the R vectors,
            which will be kept after filtering.

        Returns
        -------
        filtered_model : :py:class:`.SpinHamiltonian`
            Spin Hamiltonian after filtering.

        See Also
        --------
        filter : Modifies current object.

        Notes
        -----
        This method is not modifying the instance at which it is called.
        It creates a new instance (through deepcopy).
        """

        filtered_model = deepcopy(self)
        filtered_model.filter(
            max_distance=max_distance,
            min_distance=min_distance,
            template=template,
            R_vector=R_vector,
        )
        return filtered_model

    # TODO rewrite with new template logic (J1 through getattr)
    def form_model(self, template):
        r"""
        Force the Hamiltonian to have the symmetries of the template.

        Takes mean values of the parameters.
        For DMI interaction directions are kept the same,
        but the length of the DMI vector is scaled.

        This method modifies an instance on which it was called.

        .. versionchanged:: 0.8.0 Renamed from ``force_symmetry``

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template.

        See Also
        --------
        formed_model : Returns new object.
        """

        if not isinstance(template, ExchangeTemplate):
            raise TypeError

        self.filter(template=template)

        for name in template.names:
            bonds = template.names[name]
            symm_matrix = np.zeros((3, 3), dtype=float)
            dmi_module = 0
            for atom1, atom2, R in bonds:
                # Here names, not objects are obtained, because in general
                # template only has information about names and R.
                J = self[self.get_atom(atom1), self.get_atom(atom2), R]
                symm_matrix = symm_matrix + J.symm_matrix
                dmi_module = dmi_module + J.dmi_module

            symm_matrix = symm_matrix / len(bonds)
            dmi_module = dmi_module / len(bonds)

            for atom1, atom2, R in bonds:
                J = self[self.get_atom(atom1), self.get_atom(atom2), R]
                if J.dmi_module != 0:
                    asymm_factor = dmi_module / J.dmi_module
                else:
                    asymm_factor = 0
                J.matrix = symm_matrix + J.asymm_matrix * asymm_factor

    def formed_model(self, template):
        r"""
        Form the model from the Hamiltoian based on the template.

        Takes mean values of the parameters.
        Respect the direction of the DMI vectors.

        This method return a new instance.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template.

        Returns
        -------
        new_model : :py:class:`.SpinHamiltonian`
            Spin Hamiltonian with forced symmetry.

        See Also
        --------
        form_model: Modifies current object.
        """

        new_model = deepcopy(self)
        new_model.form_model(template=template)
        return new_model

    def dump_txt(self, *args, **kwargs):
        """
        Save the Hamiltonian in a Human-readable format.

        Parameters
        ----------
        *args
            Positional arguments for :py:func:`.dump_spinham_txt`.
        ** kwargs
            Keyword arguments for :py:func:`.dump_spinham_txt`.
        """

        dump_spinham_txt(self, *args, **kwargs)

    def dump_pickle(self, filename):
        """
        Save the Hamiltonian in a binary format.

        Parameters
        ----------
        filename : str
            Name of the file for the Hamiltonian to be saved in.
            ".pickle" is added automatically.
        """
        import pickle

        with open(f"{filename}.pickle", "wb") as file:
            pickle.dump(self, file)

    def ferromagnetic_energy(self, theta=0, phi=0):
        r"""
        Compute energy of the Hamiltonian assuming ferromagnetic state.

        Notes
        -----
        Notation of the Hamiltonian has to be known.
        See :py:meth:`.set_interpretation` and :py:attr:`.notation`.

        Parameters
        ----------
        theta : float or (N,) |array_like|_, default 0
            Angle between z axis and direction of the magnetization.
            :math:`0^{\circ} \le \theta \le 180^{\circ}`.
        phi : float or (N,) |array_like|_, default 0
            angle between x axis and projection of
            magnetization`s direction on xy plane.
            :math:`0^{\circ} \le \phi < 360^{\circ}`.

        Returns
        -------
        energy : float or (N,) :numpy:`ndarray`
            Energy of ferromagnetic Hamiltonian with magnetization direction defined
            by ``theta`` and ``phi``. In the units of J values.
        """

        # Compute spin direction
        theta = np.array(theta) * TORADIANS
        phi = np.array(phi) * TORADIANS
        if theta.shape == ():
            theta = np.array([theta])
        if phi.shape == ():
            phi = np.array([phi])
        spin_direction = np.array(
            [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)]
        )

        energy = np.zeros((3, 3), dtype=float)
        for atom1, atom2, R, J in self:
            if self.spin_normalized:
                energy += self.factor * J.matrix
            else:
                energy += self.factor * J.matrix * atom1.spin * atom2.spin

        energy = np.einsum("ni,ij,jn->n", spin_direction.T, energy, spin_direction)
        if len(energy) == 1:
            energy = energy[0]
        return energy

    def input_for_magnons(self, nodmi=False, noaniso=False, custom_mask=None):
        r"""
        Input from the spin Hamiltonian.

        This function prepare the list of exchange parameters to
        be used as an input for magnon dispersion calculation.

        Parameters
        ----------
        nodmi : bool, default=False
            If True, then DMI is not included in the dispersion.
        noaniso : bool, default=False
            If True, then anisotropy is not included in the dispersion.
        custom_mask : func
            Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
            as an input and returns (3,3) numpy:`ndarray` as an output.

        Returns
        -------
        Jij : list
        i : list
        j : list
        dij : list
        """

        Jij = []
        i = []
        j = []
        dij = []
        magnetic_atoms = self.magnetic_atoms
        atom_index = dict([(atom, i) for i, atom in enumerate(magnetic_atoms)])
        for atom1, atom2, R, J in self:
            if custom_mask is not None:
                result = custom_mask(J.matrix)
            else:
                result = J.matrix
                if nodmi:
                    result -= J.dmi_matrix
                if noaniso:
                    result -= J.aniso
            Jij.append(result)
            i.append(atom_index[atom1])
            j.append(atom_index[atom2])
            dij.append(self.get_vector(atom1, atom1, R))

        return Jij, i, j, dij


class ExchangeHamiltonian(SpinHamiltonian):
    pass


class SpinHamiltonianIterator:
    def __init__(self, exchange_model: SpinHamiltonian) -> None:
        self._bonds = list(
            map(
                lambda x: (x[0], x[1], x[2], exchange_model._bonds[x]),
                exchange_model._bonds,
            )
        )
        self._index = 0

    def __next__(self) -> Tuple[Atom, Atom, tuple, ExchangeParameter]:
        if self._index < len(self._bonds):
            result = self._bonds[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self


def dump_spinham_txt(
    spinham: SpinHamiltonian,
    filename=None,
    anisotropic=True,
    matrix=True,
    dmi=True,
    template=None,
    decimals=4,
    additional_stats=None,
):
    """
    Save the Hamiltonian in a Human-readable format.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name of the file for the Hamiltonian to be saved in.
        If not given, the Hamiltonian will be printed in the console.
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    matrix : bool, default True
        Whether to output whole matrix exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    template : :py:class:`.ExchangeTemplate`, optional
        If provided, then not the SpinHamiltonian will be written, but the model
        based on the template.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    additional_stats : str, optional
        Additional info, which will be printed right after the logo, before the
        main separator.
    """

    main_separator = "=" * 80 + "\n"
    separator = "-" * 80 + "\n"
    spinham_txt = []
    fmt = f"{decimals+4}.{decimals}f"

    spinham_txt.append(main_separator)
    spinham_txt.append(logo(date_time=True, line_length=80) + "\n")
    if additional_stats is not None:
        spinham_txt.append(additional_stats)
    spinham_txt.append(main_separator)
    spinham_txt.append(TXT_FLAGS["cell"] + "\n")
    spinham_txt.append(
        print_2d_array(
            spinham.cell,
            borders=False,
            fmt="^.8f",
            print_result=False,
            header_row=["x", "y", "z"],
        )
        + "\n"
    )
    spinham_txt.append(main_separator)
    spinham_txt.append(TXT_FLAGS["atoms"] + "\n")
    header_column = []
    header_row = [
        f"Index Name",
        "a1 (rel)",
        "a2 (rel)",
        "a3 (rel)",
        "x",
        "y",
        "z",
    ]
    atom_data = np.zeros((len(spinham.atoms), 6))
    for a_i, atom in enumerate(spinham.atoms):
        header_column.append(f"{atom.index:<5} {atom.name:<4}")
        atom_data[a_i, :3] = spinham.get_atom_coordinates(atom, relative=True)
        atom_data[a_i, 3:] = spinham.get_atom_coordinates(atom, relative=False)
    spinham_txt.append(
        print_2d_array(
            atom_data,
            borders=False,
            fmt="^.8f",
            print_result=False,
            header_row=header_row,
            header_column=header_column,
        )
        + "\n"
    )
    spinham_txt.append(main_separator)

    write_spins = True
    write_magmoms = True
    for atom in spinham.magnetic_atoms:
        try:
            atom.magmom
        except ValueError:
            write_magmoms = False
        try:
            atom.spin
            atom.spin_vector
        except ValueError:
            write_spins = False
    if write_magmoms:
        spinham_txt.append(TXT_FLAGS["magmoms"] + "\n")
        header_column = []
        header_row = [f"Index Name", "m_x", "m_y", "m_z"]
        magmom_data = np.zeros((len(spinham.magnetic_atoms), 3))
        for a_i, atom in enumerate(spinham.magnetic_atoms):
            header_column.append(f"{atom.index:<5} {atom.name:<4}")
            magmom_data[a_i] = atom.magmom
        spinham_txt.append(
            print_2d_array(
                magmom_data,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=header_row,
                header_column=header_column,
            )
            + "\n"
        )
        spinham_txt.append(main_separator)
    if write_spins:
        spinham_txt.append(TXT_FLAGS["spins"] + "\n")
        header_column = []
        header_row = [f"Index Name", "S", "S_x", "S_y", "S_z"]
        spin_data = np.zeros((len(spinham.magnetic_atoms), 4))
        for a_i, atom in enumerate(spinham.magnetic_atoms):
            header_column.append(f"{atom.index:<5} {atom.name:<4}")
            spin_data[a_i, 0] = atom.spin
            spin_data[a_i, 1:] = atom.spin_vector
        spinham_txt.append(
            print_2d_array(
                spin_data,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=header_row,
                header_column=header_column,
            )
            + "\n"
        )
        spinham_txt.append(main_separator)

    spinham_txt.append(TXT_FLAGS["notation"] + "\n")
    spinham_txt.append(f"{spinham.notation}\n")
    spinham_txt.append(spinham.notation_string + "\n")
    spinham_txt.append(main_separator)

    if template is None:
        spinham_txt.append(TXT_FLAGS["spinham"] + "\n")
        spinham_txt.append(
            f"{'Atom1':6} {'Atom2':6} (  i,   j,   k) {'J_iso':^{decimals+4}} {'Distance':^8}\n"
        )
        bonds_data = []
        for atom1, atom2, (i, j, k), J in spinham:
            data_entry = []
            distance = spinham.get_distance(atom1, atom2, (i, j, k))
            atom1 = f"{atom1.name}({atom1.index})"
            atom2 = f"{atom2.name}({atom2.index})"
            data_entry.append(separator)
            data_entry.append(
                f"{atom1:6} {atom2:6} "
                + f"({i:>3}, {j:>3}, {k:>3}) "
                + f"{J.iso:^{decimals+4}.{decimals}f} "
                + f"{distance:^8.4f}\n"
            )
            if matrix:
                data_entry.append(TXT_FLAGS["matrix"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.matrix, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if anisotropic:
                data_entry.append(TXT_FLAGS["aniso"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.aniso, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if dmi:
                data_entry.append(TXT_FLAGS["dmi_module"] + "\n")
                data_entry.append(f"  {J.dmi_module:{decimals+4}.{decimals}f}\n")
                data_entry.append(TXT_FLAGS["dmi_relative"] + "\n")
                data_entry.append(f"  {J.rel_dmi:{decimals+4}.{decimals}f}\n")
                data_entry.append(TXT_FLAGS["dmi"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.dmi, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            bonds_data.append(["".join(data_entry), distance])
        bonds_data = sorted(bonds_data, key=lambda x: x[1])
        bonds_data = [x[0] for x in bonds_data]
        spinham_txt.append("".join(bonds_data))
    else:
        spinham = spinham.formed_model(template)
        spinham_txt.append(TXT_FLAGS["spinmodel"] + "\n")
        for parameter in template.names:
            bonds = template.names[parameter]
            J = spinham[bonds[0]]
            spinham_txt.append(separator)
            spinham_txt.append(f"{parameter} contains {len(bonds)} bonds\n")
            spinham_txt.append(TXT_FLAGS["iso"])
            spinham_txt.append(f"  {J.iso:{decimals+4}.{decimals}f}\n")
            if matrix:
                spinham_txt.append(TXT_FLAGS["matrix"] + "\n")
                spinham_txt.append(
                    print_2d_array(
                        J.matrix, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if anisotropic:
                spinham_txt.append(TXT_FLAGS["aniso"] + "\n")
                spinham_txt.append(
                    print_2d_array(
                        J.aniso, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if dmi:
                spinham_txt.append(TXT_FLAGS["dmi_module"] + "\n")
                spinham_txt.append(f"  {J.dmi_module:{decimals+4}.{decimals}f}\n")
                spinham_txt.append(TXT_FLAGS["dmi_relative"] + "\n")
                spinham_txt.append(f"  {J.rel_dmi:{decimals+4}.{decimals}f}\n")
                spinham_txt.append(TXT_FLAGS["dmis"] + "\n")
                for atom1, atom2, (i, j, k) in bonds:
                    J = spinham[atom1, atom2, (i, j, k)]
                    atom1 = spinham.get_atom(atom1)
                    atom2 = spinham.get_atom(atom2)
                    atom1 = f"{atom1.name}({atom1.index})"
                    atom2 = f"{atom2.name}({atom2.index})"
                    spinham_txt.append(
                        print_2d_array(
                            J.dmi, fmt=fmt, borders=False, print_result=False, shift=2
                        )
                    )
                    spinham_txt.append(
                        f"   ({atom1:6} {atom2:6} {i:>3}, {j:>3}, {k:>3})\n"
                    )
            else:
                spinham_txt.append(TXT_FLAGS["bonds"] + "\n")
                spinham_txt.append(f"   {'Atom1':6} {'Atom2':6}   i,   j,   k\n")
                for atom1, atom2, (i, j, k) in bonds:
                    atom1 = spinham.get_atom(atom1)
                    atom2 = spinham.get_atom(atom2)
                    atom1 = f"{atom1.name}({atom1.index})"
                    atom2 = f"{atom2.name}({atom2.index})"
                    spinham_txt.append(
                        f"   {atom1:6} {atom2:6} {i:>3}, {j:>3}, {k:>3}\n"
                    )
    spinham_txt.append(main_separator)
    spinham_txt = "".join(spinham_txt)
    if filename is not None:
        file = open(filename, "w")
        file.write(spinham_txt)
        file.close()
    else:
        print(spinham_txt)
