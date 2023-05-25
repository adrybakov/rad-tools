r"""Atom class"""

import numpy as np


class Atom:
    r"""
    Atom.

    Parameters
    ----------
    literal : str, default X
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
        Meant to be unique, when an atom belongs to some group
        (i.e. to :py:class:`.Crystal` or :py:class:`.ExchangeModel`).

    Attributes
    ----------
    literal : str
    position : (3,) :numpy:`ndarray`
        Position of the atom in absolute coordinates coordinates.
    spin : (3,) :numpy:`ndarray`
    magmom : (3,) :numpy:`ndarray`
    charge : float
    """

    def __init__(
        self,
        literal="X",
        position=None,
        spin=None,
        spin_vector=None,
        magmom=None,
        charge=None,
        index=None,
    ) -> None:
        self.literal = literal
        if position is None:
            position = (0, 0, 0)
        self.position = np.array(position)
        self._spin = None
        self._spin_vector = None
        self._magmom = None
        self._charge = None
        self._index = None

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
        return self.literal

    def __format__(self, format_spec):
        return format(str(self), format_spec)

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
    def index(self):
        r"""
        Index of the atom.
        """

        if self._index is None:
            raise ValueError(f"Index is not defined for the atom {self}.")
        return self._index

    @index.setter
    def index(self, new_index):
        self._index = new_index

    @property
    def fullname(self):
        r"""Return fullname (literal + index) of an atom."""
        try:
            return f"{self.literal} #{self.index}"
        except ValueError:
            return self.literal
