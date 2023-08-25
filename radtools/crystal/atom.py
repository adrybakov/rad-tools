r"""Atom class"""

from typing import Iterable

import numpy as np

from radtools.crystal.constants import ATOM_TYPES


class Atom:
    r"""
    Atom class.

    Notes
    -----
    "==" (and "!=") operation compare two atoms based on their names and indexes.
    If index of one atom is not define, then comparison raises ``ValueError``.
    For the check of the atom type use :py:attr:`Atom.type`.
    In most cases :py:attr:`Atom.name` = :py:attr:`Atom.type`.

    Parameters
    ----------
    name : str, default X
        Name of the atom.
    position : (3,) |array_like|_, default [0,0,0]
        Position of the atom in absolute coordinates.
    spin : float or (3,) |array_like|_, optional
        Spin or spin vector of the atom.
    magmom : (3,) |array_like|_, optional
        Magnetic moment of the atom.
    charge : float, optional
        Charge of the atom.
    index : int, optional
        Custom index of an atom, used differently in different scenarios.
        Combination of :py:attr:`.name` and :py:attr:`.index`
        is meant to be unique, when an atom belongs to some group
        (i.e. to :py:class:`.Crystal` or :py:class:`.SpinHamiltonian`).
    """

    def __init__(
        self,
        name="X",
        position=None,
        spin=None,
        magmom=None,
        charge=None,
        index=None,
    ) -> None:
        # Set name
        self._name = "X"
        self.name = name

        # Set index
        self._index = None
        if index is not None:
            self.index = index

        # Set position
        self._position = np.array([0.0, 0.0, 0.0])
        if position is not None:
            self.position = np.array(position)

        # Set magmom
        self._magmom = None
        if magmom is not None:
            self.magmom = magmom

        # Set charge
        self._charge = None
        if charge is not None:
            self.charge = charge

        # Set spin
        self._spin = None
        self._spin_direction = None
        self.spin_direction = [0, 0, 1]
        if isinstance(spin, Iterable):
            self.spin_vector = spin
        elif spin is not None:
            self.spin = spin

        # Set type placeholder
        self._type = None

    def __str__(self):
        return self.name

    def __format__(self, format_spec):
        return format(str(self), format_spec)

    # ==
    def __eq__(self, other) -> bool:
        if not isinstance(other, Atom):
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for ==: '{other.__class__.__name__}' and 'Atom'"
            )
        return self.name == other.name and self.index == other.index

    def __hash__(self):
        return hash(str(self.name) + str(self.index))

    # !=
    def __neq__(self, other):
        return not self == other

    @property
    def position(self):
        r"""
        Position of the atom.

        Returns
        -------
        position : (3,) :numpy:`ndarray`
            Position of the atom in absolute coordinates.
        """
        return self._position

    @position.setter
    def position(self, new_position):
        try:
            new_position = np.array(new_position, dtype=float)
        except:
            raise ValueError(
                f"New position is not array-like, new_position = {new_position}"
            )
        if new_position.shape != (3,):
            raise ValueError(
                f"New position has to be a 3 x 1 vector, shape: {new_position.shape}"
            )
        self._position = new_position

    @property
    def name(self):
        r"""
        Name of the atom.

        Returns
        -------
        name : str
            Name of the atom.
        """
        return self._name

    @name.setter
    def name(self, new_name):
        if new_name.startswith("__") or new_name.endswith("__"):
            raise ValueError(
                f"Name of the atom ({new_name}) is not valid. It cannot start/end with '__'."
            )
        self._name = new_name
        # Reset type
        self._type = None

    @property
    def type(self):
        r"""
        Type of an atom (i.e. Cr, Ni, ...).

        Returns
        -------
        type : str
            Type of the atom.
        """
        if self._type is None:
            self._type = "X"
            for i in ATOM_TYPES:
                if i.lower() in self._name.lower():
                    self._type = i
                    if len(i) == 2:
                        break
        return self._type

    @property
    def index(self):
        r"""
        Index of an atom, meant to be unique for some group of atoms.

        Returns
        -------
        index : int
            Index of the atom.

        Raises
        ------
        ValueError
            If index is not defined for the atom.
        """

        if self._index is None:
            raise ValueError(f"Index is not defined for the atom {self}.")
        return self._index

    @index.setter
    def index(self, new_index):
        self._index = new_index

    @property
    def spin(self):
        r"""
        Spin value of the atom.

        Independent of :py:attr:`.Atom.spin_direction`.

        Returns
        -------
        spin : float
            Spin value of the atom.

        Raises
        ------
        ValueError
            If spin is not defined for the atom.
        """

        if self._spin is None:
            raise ValueError(f"Spin value is not defined for the atom {self.fullname}.")
        return self._spin

    @spin.setter
    def spin(self, new_spin):
        self._spin = float(new_spin)

    @property
    def spin_direction(self):
        r"""
        Classical spin direction of the atom.

        .. math::

            \vec{n} = (n_x, n_y, n_z), \vert \vec{n}\vert = 1

        Returns
        -------
        spin_direction : (3,) :numpy:`ndarray`
            Classical spin direction of the atom.

        Raises
        ------
        ValueError
            If spin direction is not defined for the atom.
        """

        return self._spin_direction

    @spin_direction.setter
    def spin_direction(self, new_spin_direction):
        try:
            new_spin_direction = np.array(new_spin_direction, dtype=float)
            new_spin_direction /= np.linalg.norm(new_spin_direction)
        except BufferError:
            raise ValueError(
                f"New spin direction is not array-like, new_spin_direction = {new_spin_direction}"
            )
        if new_spin_direction.shape != (3,):
            raise ValueError(
                f"New spin direction has to be a 3 x 1 vector, shape: {new_spin_direction.shape}"
            )
        self._spin_direction = new_spin_direction

    @property
    def spin_vector(self):
        r"""
        Classical spin vector of the atom.

        .. math::

            \vec{S} = (S_x, S_y, S_z), \vert \vec{S}\vert = S

        Returns
        -------
        spin_vector : (3,) :numpy:`ndarray`
            Classical spin vector of the atom.

        Raises
        ------
        ValueError
            If :py:attr:`spin` or :py:meth:`spin_direction` is not defined for the atom.
        """

        return self.spin_direction * self.spin

    @spin_vector.setter
    def spin_vector(self, new_spin_vector):
        try:
            new_spin_vector = np.array(new_spin_vector, dtype=float)
        except:
            raise ValueError(
                f"New spin vector is not array-like, new_spin_direction = {new_spin_vector}"
            )
        if new_spin_vector.shape != (3,):
            raise ValueError(
                f"New spin vector has to be a 3 x 1 vector, shape: {new_spin_vector.shape}"
            )
        self._spin_direction = new_spin_vector / np.linalg.norm(new_spin_vector)
        self._spin = np.linalg.norm(new_spin_vector)

    @property
    def magmom(self):
        r"""
        Magnetic moment of the atom.

        Implementation is fully independent of the atom`s spin.

        .. code-block:: python

            magmom = [m_x, m_y, m_z]

        units - :math:`\mu_B`

        Returns
        -------
        magmom : (3,) :numpy:`ndarray`
            Magnetic moment of the atom.
        """

        if self._magmom is None:
            raise ValueError(
                f"Magnetic moment is not defined for the atom {self.fullname}."
            )
        return self._magmom

    @magmom.setter
    def magmom(self, new_magmom):
        try:
            new_magmom = np.array(new_magmom, dtype=float)
        except:
            raise ValueError(
                f"New magnetic moment value is not array-like, new_magmom = {new_magmom}"
            )
        if new_magmom.shape != (3,):
            raise ValueError(
                f"New magnetic moment has to be a 3 x 1 vector, shape: {new_magmom.shape}"
            )
        self._magmom = new_magmom

    @property
    def charge(self):
        r"""
        Charge of the atom.

        Returns
        -------
        charge : float
            Charge of the atom.
        """

        if self._charge is None:
            raise ValueError(f"Charge is not defined for the atom {self.fullname}.")
        return self._charge

    @charge.setter
    def charge(self, new_charge):
        self._charge = float(new_charge)

    @property
    def fullname(self):
        r"""
        Fullname (name__index) of an atom.

        Double "_" is used intentionally, so the user can use "_" for
        the name of the atom.

        If index is not defined, then only name is returned.

        Returns
        -------
        fullname : str
            Fullname of the atom.

        Raises
        ------
        ValueError
            If index is not defined for the atom.
        """
        try:
            return f"{self.name}__{self.index}"
        except ValueError:
            return self.name
