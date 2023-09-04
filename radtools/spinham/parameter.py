r"""
Exchange parameter class.
"""
__all__ = ["ExchangeParameter"]
# SELF is introduced in python 3.11, it is too fresh for my taste
from typing import TypeVar

Self = TypeVar("Self", bound="ExchangeParameter")

import numpy as np

from radtools.decorate.array import print_2d_array


class ExchangeParameter:
    r"""
    Exchange parameter (:math:`\boldsymbol{J}`) class.

    Basically it is a wrapper around :numpy:`ndarray`
    with predefined exchange-specific attributes.
    Any function, which work on :numpy:`ndarray` should work on
    :py:class:`.ExchangeParameter` as well. It will act on the
    :py:attr:`.matrix` attribute.

    If ``matrix`` is specified then ``iso``, ``aniso`` and ``dmi`` are
    ignored and derived from ``matrix``. If ``matrix`` is not specified
    then it is derived from ``iso``, ``aniso`` and ``dmi``.
    If nothing is provided exchange matrix is set to zero matrix.

    Parameters
    ----------
    iso : int or float, optional
        Value of isotropic exchange parameter.
    aniso : (3, 3) |array_like|_, optional
        3 x 3 matrix of symmetric anisotropic exchange.
    dmi : (3,) |array_like|_, optional
        Dzyaroshinsky-Moria interaction vector :math:`(D_x, D_y, D_z)`.
    matrix : (3, 3) |array_like|_, optional
        Exchange matrix.
    """

    def __init__(self, matrix=None, iso=None, aniso=None, dmi=None) -> None:
        self._matrix = np.zeros((3, 3), dtype=float)
        self._iso = 0.0
        self._aniso = np.zeros((3, 3), dtype=float)
        self._dmi = np.zeros(3, dtype=float)

        if matrix is None:
            if iso is not None:
                self.iso = iso
            if dmi is not None:
                self.dmi = dmi
            if aniso is not None:
                self.aniso = aniso
        else:
            self.matrix = matrix

    def __format__(self, fmt):
        return self.__str__(fmt=fmt)

    def __str__(self, fmt=".4f"):
        return print_2d_array(self.matrix, fmt=fmt, print_result=False, borders=False)

    def __repr__(self):
        return f"ExchangeParameter({self.matrix.__repr__()})"

    @property
    def __array_interface__(self):
        return self.matrix.__array_interface__

    @property
    def T(self):
        r"""
        Transposes a matrix of the exchange parameter.

        It returns new instance of the :py:class`.ExchangeParameter`
        with transposed matrix.

        Returns
        -------
        ExchangeParameter : :py:class`.ExchangeParameter`
            New instance of the :py:class`.ExchangeParameter`
            with transposed matrix.
        """
        return ExchangeParameter(matrix=self.matrix.T)

    @property
    def matrix(self) -> np.ndarray:
        r"""
        Full exchange matrix.

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_yx, J_yy, J_yz],
             [J_zx, J_zy, J_zz]]

        Returns
        -------
        matrix : (3, 3) :numpy:`ndarray`
            Full exchange matrix.

        See Also
        --------
        symm_matrix
            Symmetric part.
        asymm_matrix
            Asymmetric part.
        """

        return self._matrix

    @matrix.setter
    def matrix(self, new_matrix):
        try:
            new_matrix = np.array(new_matrix, dtype=float)
        except:
            raise ValueError(f"New matrix is not array_like: \n{new_matrix}")
        if new_matrix.shape != (3, 3):
            raise ValueError("Matrix shape has to be equal to (3, 3)")
        self._matrix = new_matrix

    @property
    def symm_matrix(self) -> np.ndarray:
        r"""
        Symmetric part of exchange :py:attr:`.matrix`.

        .. math::

            J_{symm} = \dfrac{J + J^T}{2}

        Returns
        -------
        symm_matrix : (3, 3) :numpy:`ndarray`
            Symmetric part of exchange matrix.

        See Also
        --------
        asymm_matrix
            Asymmetric part.
        matrix
            Full matrix.
        """

        return (self.matrix + self.matrix.T) / 2

    @property
    def asymm_matrix(self) -> np.ndarray:
        """
        Asymmetric part of exchange :py:attr:`.matrix`.

        .. math::

            J_{asymm} = \dfrac{J - J^T}{2}

        Returns
        -------
        asymm_matrix : (3, 3) :numpy:`ndarray`
            Asymmetric part of exchange matrix.

        See Also
        --------
        symm_matrix
            Symmetric part.
        matrix
            Full matrix.
        """

        return (self.matrix - self.matrix.T) / 2

    @property
    def iso(self) -> float:
        r"""
        Value of isotropic exchange parameter.

        Derived from the exchange matrix (:math:`\mathbf{J}`) as

        .. math::
            J_{iso} = \dfrac{1}{3}Tr(\mathbf{J})

        Returns
        -------
        iso : float
            Value of isotropic exchange parameter.

        See Also
        --------
        iso_matrix
            Isotropic exchange in a matrix form.
        """

        return np.trace(self.symm_matrix) / 3

    @iso.setter
    def iso(self, new_iso):
        new_iso = new_iso - self.iso
        self.matrix = self.matrix + (new_iso) * np.identity(3, dtype=float)

    @property
    def iso_matrix(self) -> np.ndarray:
        r"""
        Isotropic part of the exchange :py:attr:`.matrix`.

        Matrix form:

        .. code-block:: python

            [[J, 0, 0],
             [0, J, 0],
             [0, 0, J]]

        Returns
        -------
        iso_matrix : (3, 3) :numpy:`ndarray`
            Isotropic part of the exchange matrix.

        See Also
        --------
        iso
            Isotropic exchange parameter.
        """

        return self.iso * np.identity(3, dtype=float)

    @property
    def aniso(self) -> np.ndarray:
        r"""
        3 x 3 matrix of symmetric anisotropic exchange.

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_xy, J_yy, J_yz],
             [J_xz, J_yz, J_zz]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as

        .. math::
            \mathbf{J}_{aniso} = \mathbf{J}_{symm} - \dfrac{1}{3}Tr(\mathbf{J})
            \cdot \mathbf{I}

        where :math:`\mathbf{I}` is a :math:`3\times3` identity matrix.

        Returns
        -------
        aniso : (3, 3) :numpy:`ndarray`
            Matrix of symmetric anisotropic exchange.

        See Also
        --------
        aniso_diagonal
            Diagonal part of the symmetric anisotropic exchange.
        aniso_diagonal_matrix
            Diagonal part of the symmetric anisotropic exchange in a matrix form.
        """

        return self.symm_matrix - self.iso * np.identity(3, dtype=float)

    @aniso.setter
    def aniso(self, new_aniso):
        try:
            new_aniso = np.array(new_aniso, dtype=float)
        except:
            raise ValueError(f"New aniso is not array_like: \n{new_aniso}")
        if new_aniso.shape != (3, 3):
            raise ValueError("Aniso matrix shape have " + "to be equal to (3, 3)")
        # correction is correct (see matrix update)
        new_aniso = new_aniso - self.aniso
        self.matrix = self.matrix + new_aniso

    @property
    def aniso_diagonal(self) -> np.ndarray:
        r"""
        Diagonal part of the symmetric anisotropic exchange.

        .. code-block:: python

            [J_xx, J_yy, J_zz]

        Returns
        -------
        aniso_diagonal : (3,) :numpy:`ndarray`
            Diagonal part of the symmetric anisotropic exchange.

        See Also
        --------
        aniso
            3 x 3 matrix of symmetric anisotropic exchange.
        aniso_diagonal_matrix
            Diagonal part of the symmetric anisotropic exchange in a matrix form.
        """

        return np.diag(self.aniso)

    @property
    def aniso_diagonal_matrix(self) -> np.ndarray:
        r"""
        Diagonal part of the symmetric anisotropic exchange.

        .. code-block:: python

            [[J_xx, 0, 0],
             [0, J_yy, 0],
             [0, 0, J_zz]]

        Returns
        -------
        aniso_diagonal_matrix : (3, 3) :numpy:`ndarray`
            Diagonal part of the symmetric anisotropic exchange in matrix form.

        See Also
        --------
        aniso
            3 x 3 matrix of symmetric anisotropic exchange.
        aniso_diagonal
            Diagonal part of the symmetric anisotropic exchange.
        """

        return np.diag(np.diag(self.aniso))

    @property
    def dmi(self) -> np.ndarray:
        r"""
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz).

        .. code-block:: python

            [D_x, D_y, D_z]

        Returns
        -------
        dmi : (3,) :numpy:`ndarray`
            Dzyaroshinsky-Moria interaction vector.

        See Also
        --------
        dmi_matrix
            Dzyaroshinsky-Moria interaction in a matrix form.
        dmi_module
            Module of Dzyaroshinsky-Moria vector.
        rel_dmi
            Relative strength of Dzyaroshinsky-Moria interaction.
        """

        return np.array(
            [self.asymm_matrix[1][2], self.asymm_matrix[2][0], self.asymm_matrix[0][1]],
            dtype=float,
        )

    @dmi.setter
    def dmi(self, new_dmi):
        try:
            new_dmi = np.array(new_dmi, dtype=float)
        except:
            raise ValueError(f"New DMI is not array_like: \n{new_dmi}")
        if new_dmi.shape != (3,):
            raise ValueError(f"DMI have to be a 3 component vector. {new_dmi}")
        new_dmi = new_dmi - self.dmi
        dmi_matrix = np.array(
            [
                [0, new_dmi[2], -new_dmi[1]],
                [-new_dmi[2], 0, new_dmi[0]],
                [new_dmi[1], -new_dmi[0], 0],
            ],
            dtype=float,
        )
        self.matrix = self.matrix + dmi_matrix

    @property
    def dmi_matrix(self) -> np.ndarray:
        r"""
        Asymmetric part of the exchange :py:attr:`.matrix`.

        .. code-block:: python

            [[0, D_z, -D_y],
             [-D_z, 0, D_x],
             [D_y, -D_x, 0]]

        Returns
        -------
        dmi_matrix : (3, 3) :numpy:`ndarray`
            Asymmetric part of the exchange matrix.

        See Also
        --------
        dmi
            Dzyaroshinsky-Moria interaction vector.
        dmi_module
            Module of Dzyaroshinsky-Moria vector.
        rel_dmi
            Relative strength of Dzyaroshinsky-Moria interaction.
        """

        return self.asymm_matrix

    @property
    def dmi_module(self) -> float:
        r"""
        Length of the DMI vector in the units of exchange interaction.

        Returns
        -------
        dmi_module : float
            Length of the DMI vector in the units of exchange interaction.

        See Also
        --------
        dmi
            Dzyaroshinsky-Moria interaction vector.
        dmi_matrix
            Dzyaroshinsky-Moria interaction in a matrix form.
        rel_dmi
            Relative strength of Dzyaroshinsky-Moria interaction.
        """

        return np.linalg.norm(self.dmi)

    @property
    def rel_dmi(self) -> float:
        r"""
        Relative value of DMI.

        .. math::
            \frac{\vert DMI\vert}{\vert J\vert}

        Returns
        -------
        rel_dmi : float
            Relative value of DMI.

        See Also
        --------
        dmi
            Dzyaroshinsky-Moria interaction vector.
        dmi_matrix
            Dzyaroshinsky-Moria interaction in a matrix form.
        dmi_module
            Module of Dzyaroshinsky-Moria vector.
        """

        return self.dmi_module / abs(self.iso)

    @property
    def xx(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{xx}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.xx
            1.0
            >>> J.matrix[0][0]
            1.0

        Returns
        -------
        xx : float
            Value of exchange parameter :math:`J_{xx}`.
        """
        return self.matrix[0][0]

    @xx.setter
    def xx(self, new_xx):
        self._matrix[0][0] = float(new_xx)

    @property
    def xy(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{xy}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.xy
            2.0
            >>> J.matrix[0][1]
            2.0

        Returns
        -------
        xy : float
            Value of exchange parameter :math:`J_{xy}`.
        """
        return self.matrix[0][1]

    @xy.setter
    def xy(self, new_xy):
        self._matrix[0][1] = float(new_xy)

    @property
    def xz(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{xz}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.xz
            3.0
            >>> J.matrix[0][2]
            3.0

        Returns
        -------
        xz : float
            Value of exchange parameter :math:`J_{xz}`.
        """
        return self.matrix[0][2]

    @xz.setter
    def xz(self, new_xz):
        self._matrix[0][2] = float(new_xz)

    @property
    def yx(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{yx}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.yx
            4.0
            >>> J.matrix[1][0]
            4.0

        Returns
        -------
        yx : float
            Value of exchange parameter :math:`J_{yx}`.
        """
        return self.matrix[1][0]

    @yx.setter
    def yx(self, new_yx):
        self._matrix[1][0] = float(new_yx)

    @property
    def yy(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{yy}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.yy
            5.0
            >>> J.matrix[1][1]
            5.0

        Returns
        -------
        yy : float
            Value of exchange parameter :math:`J_{yy}`.
        """
        return self.matrix[1][1]

    @yy.setter
    def yy(self, new_yy):
        self._matrix[1][1] = float(new_yy)

    @property
    def yz(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{yz}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.yz
            6.0
            >>> J.matrix[1][2]
            6.0

        Returns
        -------
        yz : float
            Value of exchange parameter :math:`J_{yz}`.
        """
        return self.matrix[1][2]

    @yz.setter
    def yz(self, new_yz):
        self._matrix[1][2] = float(new_yz)

    @property
    def zx(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{zx}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.zx
            7.0
            >>> J.matrix[2][0]
            7.0

        Returns
        -------
        zx : float
            Value of exchange parameter :math:`J_{zx}`.
        """
        return self.matrix[2][0]

    @zx.setter
    def zx(self, new_zx):
        self._matrix[2][0] = float(new_zx)

    @property
    def zy(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{zy}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.zy
            8.0
            >>> J.matrix[2][1]
            8.0

        Returns
        -------
        zy : float
            Value of exchange parameter :math:`J_{zy}`.
        """
        return self.matrix[2][1]

    @zy.setter
    def zy(self, new_zy):
        self._matrix[2][1] = float(new_zy)

    @property
    def zz(self) -> float:
        r"""
        Value of exchange parameter :math:`J_{zz}`.

        .. versionadded:: 0.8.1

        .. doctest::

            >>> from radtools.spinham import ExchangeParameter
            >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> J.zz
            9.0
            >>> J.matrix[2][2]
            9.0

        Returns
        -------
        zz : float
            Value of exchange parameter :math:`J_{zz}`.
        """
        return self.matrix[2][2]

    @zz.setter
    def zz(self, new_zz):
        self._matrix[2][2] = float(new_zz)

    # Definition of arithmetic operations5t

    # + (add)
    def __add__(self, other) -> Self:
        if isinstance(other, ExchangeParameter):
            matrix = self.matrix + other.matrix
        else:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for +: 'ExchangeParameter' and '{other.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # - (sub)
    def __sub__(self, other) -> Self:
        if isinstance(other, ExchangeParameter):
            matrix = self.matrix - other.matrix
        else:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for -: 'ExchangeParameter' and '{other.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # *
    def __mul__(self, number) -> Self:
        try:
            matrix = self.matrix * number
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for *: 'ExchangeParameter' and '{number.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # @
    def __matmul__(self, other) -> np.ndarray:
        try:
            result = self.matrix @ other
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for @: 'ExchangeParameter' and '{other.__class__.__name__}'"
            )
        return result

    # /
    def __truediv__(self, number) -> Self:
        try:
            matrix = self.matrix / number
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for /: 'ExchangeParameter' and '{number.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # //
    def __floordiv__(self, number) -> Self:
        try:
            matrix = self.matrix // number
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for //: 'ExchangeParameter' and '{number.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # %
    def __mod__(self, number) -> Self:
        try:
            matrix = self.matrix % number
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for %: 'ExchangeParameter' and '{number.__class__.__name__}'"
            )
        return ExchangeParameter(matrix=matrix)

    # *
    def __rmul__(self, number) -> Self:
        try:
            matrix = number * self.matrix
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for *: '{number.__class__.__name__}' and 'ExchangeParameter'"
            )
        return ExchangeParameter(matrix=matrix)

    # @
    def __rmatmul__(self, other) -> np.ndarray:
        try:
            result = other @ self.matrix
        except:
            raise TypeError(
                f"TypeError: unsupported operand type(s) "
                + f"for @: '{other.__class__.__name__}' and 'ExchangeParameter'"
            )
        return result

    # - (neg)
    def __neg__(self) -> Self:
        return ExchangeParameter(matrix=-self.matrix)

    # + (pos)
    def __pos__(self) -> Self:
        return ExchangeParameter(matrix=self.matrix)

    # abs()
    def __abs__(self) -> Self:
        return ExchangeParameter(matrix=np.abs(self.matrix))

    # ==
    def __eq__(self, other):
        if isinstance(other, ExchangeParameter):
            return (self.matrix == other.matrix).all()
        try:
            other = np.array(other, dtype=float)
        except:
            raise ValueError(f"Other is not array_like: \n{other}")
        if other.shape != (3, 3):
            raise ValueError("Other has to be equal to (3, 3)")
        return (self.matrix == other).all()

    # !=
    def __neq__(self, other):
        return not self == other


if __name__ == "__main__":
    J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    print(J.xy)
