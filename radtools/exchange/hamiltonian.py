r"""
Exchange module.

Write a tutorial with docstring here.
"""

from copy import deepcopy
from typing import Iterable, Tuple

import numpy as np

from radtools.crystal.atom import Atom
from radtools.crystal.crystal import Crystal
from radtools.exchange.parameter import ExchangeParameter
from radtools.exchange.template import ExchangeTemplate
from radtools.routines import toradians, print_2d_array


class NotationError(ValueError):
    r"""
    Raised when the notation (or individual property) is not defined.

    Gives a summary of the notation (or individual property) and how to set it.

    Parameters
    ----------
    name : str
        Name of the corresponding attribute.

    """

    def __init__(self, name):
        self.message = (
            f"\n\nNotation`s interpretation is not set for the property {name}.\n"
            + f"Set the notation first:\n"
            + f"    ExchangeHamiltonian.{name} = True  "
            + f"or  ExchangeHamiltonian.{name} = False\n\n"
            + f"Note: When the attribute is set for the first time it sets the interpretation, "
            + "afterwards it change the notation.\n\n"
            + f"If you want to set the interpretation again, use \n"
            + f"    ExchangeHamiltonian.set_interpretation({name} = True)"
            + "\nor\n"
            + f"    ExchangeHamiltonian.set_interpretation({name} = False)\n"
        )

    def __str__(self):
        return self.message


class ExchangeHamiltonian:
    r"""
    Exchange Hamiltonian.

    By default the notation of the exchange Hamiltonian is not defined
    and could be different in different context.
    However, in could always be checked via :py:attr:`.notation`.
    In user-specific cases it is the responsibility of the user to
    set the interpretation of the Hamiltonian`s notation.

    For the predefined notations see :py:meth:`.notation`:

    Parameters
    ----------
    crystal : :py:class:`.Crystal`, optional
        Crystal on which the exchange Hamiltonian is build.
        By default it is orthonormal lattice
        (:py:class:`.CUB`, :math:`a = 1`) with no atoms.
    notation : str or tuple of bool, optional
        One of the predefined notations or list of 5 bool.
        See :py:attr:`.notation` for details.

    Attributes
    ----------
    crystal : :py:class:`.Crystal`
        Crystal on which the ExchangeHamiltonian is build.
    """

    def __init__(self, crystal: Crystal = None, notation=None) -> None:
        self._predefined_notations = {
            "standard": (True, False, False, False, True),
            "tb2j": (True, True, False, False, True),
            "spinw": (True, False, False, False, False),
        }
        self._crystal = None
        if crystal is None:
            crystal = Crystal()
        self.crystal = crystal

        self._bonds = {}

        # Notation settings
        self._double_counting = None
        self._spin_normalized = None
        self._factor_one_half = None
        self._factor_two = None
        self._minus_sign = None
        if notation is not None:
            self.notation = notation

    def __len__(self):
        return self._bonds.__len__()

    # Notation attributes
    @property
    def notation(self):
        r"""
        Return a string with a simple comment about the Hamiltonian notation.

        It can be set with a

        * string
            One of the predefined notations: "standard", "TB2J", "SpinW"
        * iterable with 5 elements
            Each element is converted to bool and the parameters are set with values.
            Order: (double counting, spin normalized, factor 1/2, factor 2, minus sign).

        Predefined notations:

        * Standard
            (True, False, False, False, True)

            .. math::
                H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are **not** normalized.
        * |TB2J|_
            (True, True, False, False, True)

            .. math::
                H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are normalized to 1.
        * SpinW
            (True, False, False, False, False)

            .. math::
                H = \sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

            where double counting is present (:math:`ij` and :math:`ji` is in the sum).
            Spin vectors are **not** normalized.


        Setting of this attribute either change the notation
        (if corresponding attributes are defined)
        or set the interpretation (if corresponding attribute are not defined).

        Returns
        -------
        notation : (5,) tuple of bool
            ``True``/``False`` values of notation properties:

            .. code-block:: python

                (double_counting, spin_normalized, factor_one_half, factor_two, minus_sign)

        Examples
        --------
        Setting one of the predefined notations:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.ExchangeHamiltonian()
            >>> model.notation = "standard"
            >>> model.notation
            H = -sum_{i,j} S_i J_ij S_j
            Double counting is present.
            Spin vectors are not normalized.
            (True, False, False, False, True)
            >>> model.notation = "TB2J"
            >>> model.notation
            H = -sum_{i,j} S_i J_ij S_j
            Double counting is present.
            Spin vectors are normalized to 1.
            (True, True, False, False, True)
            >>> model.notation = "SpinW"
            >>> model.notation
            H = sum_{i,j} S_i J_ij S_j
            Double counting is present.
            Spin vectors are not normalized.
            (True, False, False, False, False)

        Setting the notation:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.ExchangeHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
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
            >>> model = rad.ExchangeHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # For the first time interpretation is set,
            >>> # values of exchange are not changed
            >>> model.minus_sign = True
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> # Once the property is set the values
            >>> # are changing if the property is changed again.
            >>> model.minus_sign = False
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0

        Changing individual properties:

        .. doctest::

            >>> import radtools as rad
            >>> model = rad.ExchangeHamiltonian()
            >>> Cr = rad.Atom("Cr", spin=1.5)
            >>> model.add_bond(rad.ExchangeParameter(iso=1), Cr, Cr, (1, 0, 0))
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> model.notation = "standard"
            >>> model[Cr, Cr, (1, 0, 0)].iso
            1.0
            >>> model.double_counting, model.spin_normalized, model.factor_one_half, model.factor_two, model.minus_sign
            (True, False, False, False, True)
            >>> model.minus_sign = False
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0
            >>> model.factor_one_half, model.factor_two
            (False, False)
            >>> model.factor_one_half = True
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -2.0
            >>> model.factor_one_half, model.factor_two
            (True, False)
            >>> model.factor_two = True
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -1.0
            >>> # Note that the values are switched to False,
            >>> # since factor one half and two are cancelling each other
            >>> model.factor_one_half, model.factor_two
            (False, False)
            >>> model.spin_normalized = True
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -2.25
            >>> model.double_counting = False
            >>> model[Cr, Cr, (1, 0, 0)].iso
            -4.5


        See Also
        --------
        double_counting
        spin_normalized
        factor_one_half
        factor_two
        minus_sign
        set_interpretation
        """

        text_result = "H = "
        if self.minus_sign:
            text_result += "-"

        if self.factor_one_half and not self._factor_two:
            text_result += "1/2 "
        if self._factor_two and not self.factor_one_half:
            text_result += "2 "

        text_result += "sum_{"
        if self.double_counting:
            text_result += "i,j} "
        else:
            text_result += "i>=j} "

        text_result += "S_i J_ij S_j\n"

        if self.double_counting:
            text_result += "Double counting is present.\n"
        else:
            text_result += "No double counting.\n"

        if self.spin_normalized:
            text_result += "Spin vectors are normalized to 1."
        else:
            text_result += "Spin vectors are not normalized."

        print(text_result)

        return (
            self.double_counting,
            self.spin_normalized,
            self.factor_one_half,
            self.factor_two,
            self.minus_sign,
        )

    @notation.setter
    def notation(self, new_notation):
        if isinstance(new_notation, str):
            new_notation = new_notation.lower()
            if new_notation not in self._predefined_notations:
                raise ValueError(
                    f"Predefine notations are: "
                    + f"{list(self._predefined_notations)}, got: {new_notation}"
                )
            new_notation = self._predefined_notations[new_notation]
        elif isinstance(new_notation, Iterable) and len(new_notation) == 5:
            new_notation = tuple(map(bool, new_notation))
        else:
            raise ValueError(
                "New notation has to be either a string "
                + "or an iterable with five elements, "
                + f"got: {new_notation}"
            )
        (
            self.double_counting,
            self.spin_normalized,
            self.factor_one_half,
            self.factor_two,
            self.minus_sign,
        ) = new_notation

    def set_interpretation(
        self,
        double_counting: bool = None,
        spin_normalized: bool = None,
        factor_one_half: bool = None,
        factor_two: bool = None,
        minus_sign: bool = None,
    ):
        r"""
        Set the interpretation of the Hamiltonian`s notation.

        It is not changing the J values,
        but rather indicate how to interpret the J values of the Hamiltonian.

        Parameters
        ----------
        double_counting : bool, optional
        spin_normalized : bool, optional
        factor_one_half : bool, optional
        factor_two : bool, optional
        minus_sign : bool, optional

        See Also
        --------
        double_counting
        spin_normalized
        factor_one_half
        factor_two
        minus_sign
        """

        if double_counting is not None:
            self._double_counting = bool(double_counting)
        if spin_normalized is not None:
            self._spin_normalized = bool(spin_normalized)
        if factor_one_half is not None:
            self._factor_one_half = bool(factor_one_half)
        if factor_two is not None:
            self._factor_two = bool(factor_two)
        if minus_sign is not None:
            self._minus_sign = bool(minus_sign)

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
                self.add_bond(J.T, atom2, atom1, (-i, -j, -k))

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
                    self.add_bond(J.T, atom2, atom1, (-i, -j, -k))
                    self.remove_bond(atom1, atom2, (i, j, k))

    @double_counting.setter
    def double_counting(self, new_value: bool):
        if self._double_counting is not None:
            factor = 1
            if self._double_counting and not new_value:
                factor = 2
            elif not self._double_counting and new_value:
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
        self._double_counting = bool(new_value)

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
            if self._spin_normalized ^ new_value:
                for atom1, atom2, R, J in self:
                    if self._spin_normalized and not new_value:
                        self[atom1, atom2, R] = J / atom1.spin / atom2.spin
                    elif not self._spin_normalized and new_value:
                        self[atom1, atom2, R] = J * atom1.spin * atom2.spin
        self._spin_normalized = bool(new_value)

    @property
    def factor_one_half(self) -> bool:
        r"""
        Whether factor 1/2 is present in the Hamiltonian.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        factor_one_half : bool
            ``True`` if factor 1/2 is present, ``False`` otherwise.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._factor_one_half is None:
            raise NotationError("factor_one_half")
        return self._factor_one_half

    @factor_one_half.setter
    def factor_one_half(self, new_value: bool):
        if self._factor_one_half is not None:
            factor = 1
            if self._factor_one_half and not new_value:
                factor = 0.5
            elif not self._factor_one_half and new_value:
                factor = 2
            if factor != 1:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * factor
        self._factor_one_half = bool(new_value)
        if self._factor_one_half and self.factor_two:
            self._factor_one_half = False
            self._factor_two = False

    @property
    def factor_two(self) -> bool:
        r"""
        Whether factor 2 is present in the Hamiltonian.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        factor_two : bool
            ``True`` if factor 2 is present, ``False`` otherwise.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._factor_two is None:
            raise NotationError("factor_two")
        return self._factor_two

    @factor_two.setter
    def factor_two(self, new_value: bool):
        if self._factor_two is not None:
            factor = 1
            if self._factor_two and not new_value:
                factor = 2
            elif not self._factor_two and new_value:
                factor = 0.5
            if factor != 1:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * factor
        self._factor_two = bool(new_value)
        if self._factor_one_half and self.factor_two:
            self._factor_one_half = False
            self._factor_two = False

    @property
    def minus_sign(self) -> bool:
        r"""
        Whether the minus sign is present in the Hamiltonian.

        Access this attribute to check the current notation of the Hamiltonian,
        set this attribute to change the notation of the Hamiltonian.

        Returns
        -------
        minus_sign : bool
            ``True`` if minus sign is present, ``False`` otherwise.

        Raises
        ------
        :py:exc:`.NotationError`
            If the interpretation of the Hamiltonian`s notation is not set.

        See Also
        --------
        set_interpretation : to change the interpretation of the Hamiltonian`s notation
            without the modification of the parameters.
        """

        if self._minus_sign is None:
            raise NotationError("minus_sign")
        return self._minus_sign

    @minus_sign.setter
    def minus_sign(self, new_value: bool):
        if self._minus_sign is not None:
            factor = 1
            if self._minus_sign ^ new_value:
                factor = -1

            if factor != 1:
                for atom1, atom2, R, J in self:
                    self[atom1, atom2, R] = J * factor
        self._minus_sign = bool(new_value)

    def __iter__(self):
        return ExchangeHamiltonianIterator(self)

    def __contains__(self, key):
        atom1, atom2, R = key
        if isinstance(atom1, str):
            atom1 = self.crystal.get_atom(atom1)
        if isinstance(atom2, str):
            atom2 = self.crystal.get_atom(atom2)
        key = (atom1, atom2, R)
        return key in self._bonds

    def __getitem__(self, key) -> ExchangeParameter:
        atom1, atom2, R = key
        if isinstance(atom1, str):
            atom1 = self.crystal.get_atom(atom1)
        if isinstance(atom2, str):
            atom2 = self.crystal.get_atom(atom2)
        key = (atom1, atom2, R)
        return self._bonds[key]

    def __getattr__(self, name):
        # Fix copy/deepcopy RecursionError
        if name in ["__setstate__"]:
            raise AttributeError(name)
        return getattr(self.crystal, name)

    @property
    def crystal(self) -> Crystal:
        r"""
        Crystal of the Hamiltonian.

        Crystal, which define the structure.
        See :py:class:`.Crystal`.

        Returns
        -------
        crystal : :py:class:`.Crystal`
            Crystal of the Hamiltonian.
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
        for atom1, atom2, R, J in self:
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

    def __setitem__(self, key, value):
        self.add_bond(value, *key)

    def add_bond(self, J: ExchangeParameter, atom1: Atom, atom2: Atom, R):
        r"""
        Add one bond to the Hamiltonian.

        More straightforward syntax is advised:

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.ExchangeHamiltonian(rad.Crystal())
            >>> model[Cr, Cr, (1,0,0)] = J
            >>> (Cr, Cr, (1,0,0)) in model
            True

        It is the same as

        .. doctest::

            >>> import radtools as rad
            >>> Cr = rad.Atom("Cr")
            >>> J = rad.ExchangeParameter(iso=1)
            >>> model = rad.ExchangeHamiltonian(rad.Crystal())
            >>> model.add_bond(J, Cr, Cr, (1,0,0))
            >>> (Cr, Cr, (1,0,0)) in model
            True

        Parameters
        ----------
        J : :py:class:`.ExchangeParameter`
            An instance of :py:class:`ExchangeParameter`.
        atom1 : :py:class:`Atom` or str
            Atom object in (0, 0, 0) unit cell.
            str works only if atom is already in the Hamiltonian.
        atom2 : :py:class:`Atom` or str
            Atom object in R unit cell.
            str works only if atom is already in the Hamiltonian.
        R : tuple of ints
            Vector of the unit cell for atom2.
            In the relative coordinates (i,j,k).
        """

        if isinstance(atom1, str):
            atom1 = self.crystal.get_atom(atom1)
        elif atom1 not in self.crystal:
            self.add_atom(atom1)

        if isinstance(atom2, str):
            atom2 = self.crystal.get_atom(atom2)
        elif atom2 not in self.crystal:
            self.add_atom(atom2)

        self._bonds[(atom1, atom2, R)] = J

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
            >>> model = rad.ExchangeHamiltonian(rad.Crystal())
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
            >>> model = rad.ExchangeHamiltonian(rad.Crystal())
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

        if isinstance(atom1, str):
            atom1 = self.crystal.get_atom(atom1)

        if isinstance(atom2, str):
            atom2 = self.crystal.get_atom(atom2)

        try:
            del self._bonds[(atom1, atom2, R)]
        except KeyError:
            raise KeyError(
                f"Bond ({atom2.fullname}, {atom2.fullname}, {R}) is not present in the model."
            )

    def add_atom(self, atom: Atom):
        r"""
        Add atom to the Hamiltonian.

        see :py:meth:`.Crystal.add_atom`

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        """

        self.crystal.add_atom(atom)

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

        if isinstance(atom, str):
            atom = self.crystal.get_atom(atom)
        bonds_for_removal = []
        for atom1, atom2, R, J in self:
            if atom1 == atom or atom2 == atom:
                bonds_for_removal.append((atom1, atom2, R))

        for bond in bonds_for_removal:
            del self[bond]

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
            dis = self.crystal.get_distance(atom1, atom2, R)

            if max_distance is not None and dis > max_distance:
                bonds_for_removal.add((atom1, atom2, R))

            if min_distance is not None and dis < min_distance:
                bonds_for_removal.add((atom1, atom2, R))

            if R_vector is not None and R not in R_vector:
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
            del self[atom1, atom2, R]

    def filtered(
        self, max_distance=None, min_distance=None, template=None, R_vector=None
    ):
        r"""
        Create filtered exchange Hamiltonian based on the given conditions.

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
        filtered_model : :py:class:`.ExchangeHamiltonian`
            Exchange Hamiltonian after filtering.

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
    def force_symmetry(self, template):
        r"""
        Force the Hamiltonian to have the symmetries of the template.

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
                # Here names, not objects are obtained, because in general
                # template only has information about names and R.
                J = self[self.crystal.get_atom(atom1), self.crystal.get_atom(atom2), R]
                symm_matrix = symm_matrix + J.symm_matrix
                dmi_module = dmi_module + J.dmi_module

            symm_matrix = symm_matrix / len(bonds)
            dmi_module = dmi_module / len(bonds)

            for atom1, atom2, R in bonds:
                J = self[self.crystal.get_atom(atom1), self.crystal.get_atom(atom2), R]
                if J.dmi_module != 0:
                    asymm_factor = dmi_module / J.dmi_module
                else:
                    asymm_factor = 0
                J.matrix = symm_matrix + J.asymm_matrix * asymm_factor

    def forced_symmetry(self, template):
        r"""
        Force the Hamiltonian to have the symmetries of the template.

        Takes mean values of the parameters.
        Respect the direction of the DMI vectors.

        This method return a new instance.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template.

        Returns
        -------
        new_model : :py:class:`.ExchangeHamiltonian`
            Exchange Hamiltonian with forced symmetry.

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
        matrix=False,
        dmi=False,
    ):
        r"""
        Return exchange Hamiltonian based on the template file in .txt format.

        Parameters
        ----------
        template : :py:class:`.ExchangeTemplate`
            Template of the desired exchange Hamiltonian.
            (see :ref:`rad-make-template`)
        accuracy : int, default 4
            Accuracy for the exchange values
        force_symmetry: bool, default False
            Whether to force the symmetry of the template on the exchange Hamiltonian.
            If ``False`` then each individual bond is written, otherwise
            exchange parameters of the template are written.
        isotropic : bool, default False
            Whether to output isotropic exchange.
        anisotropic : bool, default False
            Whether to output anisotropic exchange.
        matrix : bool, default False
            Whether to output whole matrix exchange.
        dmi : bool, default False
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
        summary = ""
        for name in template.names:
            scalar_written = False
            for atom1, atom2, R in template.names[name]:
                # Get values from the model
                atom1 = self.crystal.get_atom(atom1)
                atom2 = self.crystal.get_atom(atom2)
                J = self[atom1, atom2, R]

                if not force_symmetry:
                    # Write the values
                    summary += (
                        f"{atom1:3} {atom2:3} "
                        + f"({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n"
                    )
                    if isotropic:
                        summary += f"    Isotropic: {J.iso:.{decimals}f}\n"
                    if anisotropic:
                        summary += (
                            f"    Anisotropic:\n"
                            + print_2d_array(
                                J.aniso,
                                fmt=f"{decimals+3}.{decimals}f",
                                borders=False,
                                shift=7,
                                print_result=False,
                            )
                            + "\n"
                        )
                    if matrix:
                        summary += (
                            f"    Matrix:\n"
                            + print_2d_array(
                                J.matrix,
                                fmt=f"{decimals+3}.{decimals}f",
                                borders=False,
                                shift=7,
                                print_result=False,
                            )
                            + "\n"
                        )
                    if dmi:
                        summary += (
                            f"    |DMI|: {J.dmi_module:.{decimals}f}\n"
                            + f"    |DMI/J|: {J.rel_dmi:.{decimals}f}\n"
                            + f"    DMI: {J.dmi[0]:.{decimals}f} "
                            + f"{J.dmi[1]:.{decimals}f} "
                            + f"{J.dmi[2]:.{decimals}f}\n"
                        )
                    summary += "\n"
                else:
                    if not scalar_written:
                        scalar_written = True
                        summary += f"{name}\n"
                        if isotropic:
                            summary += f"    Isotropic: {J.iso:.{decimals}f}\n"
                        if anisotropic:
                            summary += (
                                f"    Anisotropic:\n"
                                + print_2d_array(
                                    J.aniso,
                                    fmt=f"{decimals+3}.{decimals}f",
                                    borders=False,
                                    shift=7,
                                    print_result=False,
                                )
                                + "\n"
                            )
                        if matrix:
                            summary += (
                                f"    Matrix:\n"
                                + print_2d_array(
                                    J.matrix,
                                    fmt=f"{decimals+3}.{decimals}f",
                                    borders=False,
                                    shift=7,
                                    print_result=False,
                                )
                                + "\n"
                            )
                        if dmi:
                            summary += (
                                f"    |DMI|: {J.dmi_module:.{decimals}f}\n"
                                + f"    |DMI/J|: {J.rel_dmi:.{decimals}f}\n"
                            )
                    if dmi:
                        summary += (
                            f"    DMI: {J.dmi[0]:{decimals+3}.{decimals}f} "
                            + f"{J.dmi[1]:{decimals+3}.{decimals}f} "
                            + f"{J.dmi[2]:{decimals+3}.{decimals}f} "
                            + f"({atom1:3} {atom2:3} "
                            + f"{R[0]:2.0f} {R[1]:2.0f} {R[2]:2.0f})\n"
                        )
            if force_symmetry:
                summary += "\n"

        return summary

    def ferromagnetic_energy(self, theta=0, phi=0):
        r"""
        Compute energy of the Hamiltonian assuming ferromagnetic state.

        Notes
        -----
        Notation of the Hamiltonian nas to be known.
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

        theta = np.array(theta) * toradians
        phi = np.array(phi) * toradians
        if theta.shape == ():
            theta = np.array([theta])
        if phi.shape == ():
            phi = np.array([phi])
        spin_direction = np.array(
            [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)]
        )

        energy = np.zeros((3, 3), dtype=float)
        factor = 1
        if self.minus_sign:
            factor *= -1
        if self.factor_one_half:
            factor /= 2
        if self.factor_two:
            factor *= 2
        for atom1, atom2, R, J in self:
            if self.spin_normalized:
                energy += factor * J.matrix
            else:
                energy += factor * J.matrix * atom1.spin * atom2.spin

        energy = np.einsum("ni,ij,jn->n", spin_direction.T, energy, spin_direction)
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


class ExchangeHamiltonianIterator:
    def __init__(self, exchange_model: ExchangeHamiltonian) -> None:
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


if __name__ == "__main__":
    model = ExchangeHamiltonian()
    Cr = Atom("Cr", (0, 0, 0), spin=3 / 2)
    model[Cr, Cr, (1, 0, 0)] = ExchangeParameter(iso=1)
    assert model[Cr, Cr, (1, 0, 0)].iso == 1

    model.notation = "standard"
    assert model[Cr, Cr, (1, 0, 0)].iso == 1
    model.notation = "standard"
    assert model[Cr, Cr, (1, 0, 0)].iso == 1

    assert model.double_counting
    model.double_counting = False
    assert model[Cr, Cr, (1, 0, 0)].iso == 2
    model.double_counting = True
    assert model[Cr, Cr, (1, 0, 0)].iso == 1
