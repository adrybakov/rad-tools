r"""Atom class"""

import numpy as np

from radtools.crystal.atom_types import atom_types


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
    spin : float, optional
        Spin of the atom.
    spin_vector : (3,) |array_like|_, optional
        Classical spin vector of the atom.
        By default oriented alond z axis.
        Only the direction matters, for the value use ``spin``.
    magmom : (3,) |array_like|_, optional
        Magnetic moment of the atom.
    charge : float, optional
        Charge of the atom.
    index : int, optional
        Custom index of an atom, used differently in different scenarios.
        Combination of :py:attr:`.name` and :py:attr:`.index`
        is meant to be unique, when an atom belongs to some group
        (i.e. to :py:class:`.Crystal` or :py:class:`.ExchangeHamiltonian`).

    Attributes
    ----------
    name : str
    type : str
    index : int
    position : (3,) :numpy:`ndarray`
        Position of the atom in absolute coordinates.
    spin : (3,) :numpy:`ndarray`
    magmom : (3,) :numpy:`ndarray`
    charge : float
    """

    def __init__(
        self,
        name="X",
        position=None,
        spin=None,
        spin_vector=None,
        magmom=None,
        charge=None,
        index=None,
    ) -> None:
        self._name = "X"
        self._index = None
        self._type = None
        if position is None:
            position = (0, 0, 0)
        self.position = np.array(position)
        self._spin = None
        self._spin_vector = None
        self._magmom = None
        self._charge = None

        self.name = name
        if spin is not None:
            self.spin = spin
        if spin_vector is not None:
            self.spin_vector = spin_vector
        if magmom is not None:
            self.magmom = magmom
        if charge is not None:
            self.charge = charge
        if index is not None:
            self.index = index

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
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name
        self._type = None

    @property
    def type(self):
        r"""
        Type of an atom (i.e. Cr, Ni, ...).
        """
        if self._type is None:
            self._type = "X"
            for i in atom_types:
                if i in self._name:
                    self._type = i
                    break
        return self._type

    @property
    def index(self):
        r"""
        Index of an atom, meant to be unique for some group of atoms.
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
        Classical spin vector of the atom.

        .. code-block:: python

            S_vec = [S_x, S_y, S_z]
            S = np.linalg.norm(np.array(S))
        """

        if self._spin is None:
            raise ValueError(f"Spin is not defined for the atom {self.fullname}.")
        return self._spin

    @spin.setter
    def spin(self, new_spin):
        self._spin = new_spin

    @property
    def spin_vector(self):
        r"""
        Classical spin vector of the atom.

        If :py:attr:`spin` is set:

        .. math::

            \vec{S} = (S_x, S_y, S_z)

        If :py:attr:`spin` is not set:

        .. math::

            \vec{n} = (n_x, n_y, n_z): \vec{n} \parallel \vec{S}, \vert \vec{n}\vert = 1
        """

        if self._spin_vector is None:
            raise ValueError(
                f"Spin vector is not defined for the atom {self.fullname}."
            )
        if self._spin is not None:
            return self._spin * self._spin_vector
        return self._spin_vector

    @spin_vector.setter
    def spin_vector(self, new_spin_vector):
        try:
            new_spin_vector = np.array(new_spin_vector)
            new_spin_vector /= np.linalg.norm(new_spin_vector)
        except:
            raise ValueError(
                f"New spin vector is not array-like, new_spin = {new_spin_vector}"
            )
        if new_spin_vector.shape != (3,):
            raise ValueError(
                f"New spin vector has to be a 3 x 1 vector, shape: {new_spin_vector.shape}"
            )
        self._spin_vector = new_spin_vector

    @property
    def magmom(self):
        r"""
        Magnetic moment of the atom.

        .. code-block:: python

            magmom = [m_x, m_y, m_z]
        """

        if self._magmom is None:
            raise ValueError(
                f"Magnetic moment is not defined for the atom {self.fullname}."
            )
        return self._magmom

    @magmom.setter
    def magmom(self, new_magmom):
        try:
            new_magmom = np.array(new_magmom)
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
        """

        if self._charge is None:
            raise ValueError(f"Charge is not defined for the atom {self.fullname}.")
        return self._charge

    @charge.setter
    def charge(self, new_charge):
        self._charge = new_charge

    @property
    def fullname(self):
        r"""Return fullname (name + index) of an atom."""
        try:
            return f"{self.name}_{self.index}"
        except ValueError:
            return self.name
