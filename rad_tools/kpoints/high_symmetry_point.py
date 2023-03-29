from math import cos, sin, tan

import numpy as np


class HighSymmetryPoints:
    r"""
    Generator of the high symmetry k-points.

    Parameters
    ----------
    kpoints : dict
        Dictionary of the high symmetry k-points and coordinate
        (fractions of the reciprocal vectors).

        .. code-block:: python

            kpoints = {"name" : [xb1, xb2, xb3]}
    """

    def __init__(self, kpoints=None, path=None) -> None:
        self._kpoints = kpoints
        self._path = path

        self._PLOT_LITERALS = {
            "Gamma": "$\\Gamma$",
            "M": "$M$",
            "R": "$R$",
            "X": "$X$",
            "K": "$K$",
            "L": "$L$",
            "U": "$U$",
            "W": "$W$",
            "H": "$H$",
            "P": "$P$",
            "N": "$N$",
            "A": "$A$",
            "Z": "$Z$",
            "Z1": "$Z_1$",
            "Sigma": "$\\Sigma$",
            "Sigma1": "$\\Sigma_1$",
            "Y": "$Y$",
            "Y1": "$Y_1$",
            "S": "$S$",
            "T": "$T$",
            "A1": "$A_1$",
            "X1": "$X_1$",
            "C": "$C$",
            "C1": "$C_1$",
            "D": "$D$",
            "D1": "$D_1$",
            "H1": "$H_1$",
            "L1": "$L_1$",
            "L2": "$L_2$",
            "B": "$B$",
            "B1": "$B_1$",
            "F": "$F$",
            "P1": "$P_1$",
            "P2": "$P_2$",
            "Q": "$Q$",
            "Q1": "$Q_1$",
            "E": "$E$",
            "H2": "$H_2$",
            "M1": "$M_1$",
            "M2": "$M_2$",
            "N1": "$N_1$",
            "F1": "$F_1$",
            "F2": "$F_2$",
            "F3": "$F_3$",
            "I": "$I$",
            "I1": "$I_1$",
            "X2": "$X_2$",
            "Y2": "$Y_2$",
            "Y3": "$Y_3$",
        }

    @property
    def kpoints(self):
        r"""
        Dictionary of K-points.

        .. code-block:: python

            kpoints = {Point1 : [k1, k2, k3], ...}

        where ``k1``, ``k2``, ``k3`` are
        the relative coordinates of the ``Point1``.
        """

        if self._kpoints is None:
            return {}
        else:
            return self._kpoints

    def add_kpoint(self, name, coordinates, plot_name=None):
        r"""
        Add one kpoint to the set of High symmetry kpoints.

        Parameters
        ----------
        name : str
            Name of the kpoints, which is used by the code as a key.
        coordinates : 1x3 array
            Relative coordinates of the kpoint.
        plot_name : str
            name of the kpoints to be display on the plot.
            Equal to the ``name`` by default.
        """

        coordinates = np.array(coordinates)
        if coordinates.shape != (3,):
            raise ValueError
        if self._kpoints is None:
            self._kpoints = {}
        self._kpoints[name] = coordinates

        if plot_name is None and name not in self._PLOT_LITERALS:
            self._PLOT_LITERALS[name] = name
        elif plot_name is not None:
            self._PLOT_LITERALS[name] = plot_name

    @property
    def path(self):
        r"""
        Path of K-points.

        Always is a list of lists.

        .. code-block:: python

            path = [subpath1, subpath2, ...]
            subpath1 = [hs_kpoint1, hs_kpoint2, ...]
        """

        if self._path is None:
            return []
        else:
            return self._path

    @path.setter
    def path(self, new_path):
        if new_path is None:
            self._path = []

        if isinstance(new_path, str):
            if "|-" in new_path or "-|" in new_path:
                raise ValueError("Check the format of the new path")
            subpaths = new_path.split("|")
            if subpaths[-1] == "":
                subpaths = subpaths[0:-1]
            if subpaths[0] == "":
                subpaths = subpaths[1:]

            for i in range(0, len(subpaths)):
                subpaths[i] = subpaths[i].split("-")
                for j in range(0, len(subpaths[i])):
                    if subpaths[i][j] not in self.kpoints:
                        raise ValueError(
                            f"Provided  high symmetry "
                            + f"k-point {subpaths[i][j]} is not known."
                        )
            self._path = subpaths

        if isinstance(new_path, list):
            if len(new_path) == 0:
                self._path = []
            elif isinstance(new_path[0], list):
                for i in range(0, len(new_path)):
                    for j in range(0, len(new_path[i])):
                        if new_path[i][j] not in self.kpoints:
                            raise ValueError(
                                f"Provided  high symmetry "
                                + f"k-point {new_path[i][j]} is not known."
                            )
                if len(new_path) == 1 and len(new_path[0]) == 0:
                    self._path = []
                else:
                    self._path = new_path
            elif isinstance(new_path[0], str):
                for i in range(0, len(new_path)):
                    for j in range(0, len(new_path[i])):
                        if new_path[i][j] not in self.kpoints:
                            raise ValueError(
                                f"Provided  high symmetry "
                                + f"k-point {new_path[i][j]} is not known."
                            )
                self._path = [new_path]

            else:
                raise ValueError(f"Path format is not supported: {new_path}")

    def cub(self):
        r"""Cubic (CUB, cP)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
        }

        path = [["Gamma", "X", "M", "Gamma", "R", "X"], ["M", "R"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def fcc(self):
        r"""Face-centred cubic (FCC, cF)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "K": np.array([3 / 8, 3 / 8, 3 / 4]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "U": np.array([5 / 8, 1 / 4, 5 / 8]),
            "W": np.array([1 / 2, 1 / 4, 3 / 4]),
            "X": np.array([1 / 2, 0, 1 / 2]),
        }

        path = [["Gamma", "X", "W", "K", "Gamma", "L", "U", "W", "L", "K"], ["U", "X"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def bcc(self):
        r"""Body-centered cubic (BCC, cl)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "H": np.array([1 / 2, -1 / 2, 1 / 2]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "N": np.array([0, 0, 1 / 2]),
        }

        path = [["Gamma", "H", "N", "Gamma", "P", "H"], ["P", "N"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def tet(self):
        r"""Tetragonal (TET, tP)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([0, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "X", "M", "Gamma", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def bct1(self, a, c):
        r"""Body-centred tetragonal (BCT, tI),

        .. math::
            c < a
        """

        eta = (1 + c**2 / a**2) / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "M": np.array([-1 / 2, 1 / 2, 1 / 2]),
            "N": np.array([0, 1 / 2, 0]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "X": np.array([0, 0, 1 / 2]),
            "Z": np.array([eta, eta, -eta]),
            "Z1": np.array([-eta, 1 - eta, eta]),
        }

        path = [["Gamma", "X", "M", "Gamma", "Z", "P", "N", "Z1", "M"], ["X", "P"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def bct2(self, a, c):
        r"""Body-centred tetragonal (BCT, tI),

        .. math::
            c > a
        """

        eta = (1 + c**2 / a**2) / 4
        zeta = a**2 / (2 * c**2)
        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "N": np.array([0, 1 / 2, 0]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "Sigma": np.array([-eta, eta, eta]),
            "Sigma1": np.array([eta, 1 - eta, -eta]),
            "X": np.array([0, 0, 1 / 2]),
            "Y": np.array([-zeta, zeta, 1 / 2]),
            "Y1": np.array([1 / 2, 1 / 2, -zeta]),
            "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
        }

        path = [
            ["gamma", "X", "Y", "Sigma", "Gamma", "Z", "Sigma1", "N", "P", "Y1", "Z"],
            ["X", "P"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orc(self):
        r"""Orthorhombic (ORC, oP)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "S": np.array([1 / 2, 1 / 2, 0]),
            "T": np.array([0, 1 / 2, 1 / 2]),
            "U": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "X", "S", "Y", "Gamma", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orcf1(self, a, b, c):
        r"""Face-centred orthorhombic (ORCF, oF),

        .. math::
            \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}
        """

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
            "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "T": np.array([1, 1 / 2, 1 / 2]),
            "X": np.array([0, eta, eta]),
            "X1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }

        path = [
            ["Gamma", "Y", "T", "Z", "Gamma", "X", "A1", "Y"],
            ["T", "X1"],
            ["X", "A", "Z"],
            ["L", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orcf2(self, a, b, c):
        r"""Face-centred orthorhombic (ORCF, oF),

        .. math::
            \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}
        """

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        delta = (1 + b**2 / a**2 - b**2 / c**2) / 4
        phi = (1 + c**2 / b**2 - c**2 / a**2) / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "C": np.array([1 / 2, 1 / 2 - eta, 1 - eta]),
            "C1": np.array([1 / 2, 1 / 2 + eta, 1 - eta]),
            "D": np.array([1 / 2 - delta, 1 / 2, 1 - delta]),
            "D1": np.array([1 / 2 + delta, 1 / 2, delta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "H": np.array([1 - phi, 1 / 2 - phi, 1 / 2]),
            "H1": np.array([phi, 1 / 2 + phi, 1 / 2]),
            "X": np.array([0, 1 / 2, 1 / 2]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }

        path = [
            ["Gamma", "Y", "C", "D", "X", "Gamma", "Z", "D1", "H", "C"],
            ["C1", "Z"],
            ["X", "H1"],
            ["H", "Y"],
            ["L", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orcf3(self, a, b, c):
        r"""Face-centred orthorhombic (ORCF, oF),

        .. math::
            \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}
        """

        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
            "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "T": np.array([1, 1 / 2, 1 / 2]),
            "X": np.array([0, eta, eta]),
            "X1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }

        path = [
            ["Gamma", "Y", "T", "Z", "Gamma", "X", "A1", "Y"],
            ["X", "A", "Z"],
            ["L", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orci(self, a, b, c):
        r"""Body-centred orthorhombic (ORCI, oI)"""

        eta = (1 + a**2 / c**2) / 4
        zeta = (1 + b**2 / c**2) / 4
        delta = (b**2 - a**2) / (4 * c**2)
        mu = (a**2 + b**2) / (4 * c**2)

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "L": np.array([-mu, mu, 1 / 2 - delta]),
            "L1": np.array([mu, -mu, 1 / 2 + delta]),
            "L2": np.array([1 / 2 - delta, 1 / 2 + delta, -mu]),
            "R": np.array([0, 1 / 2, 0]),
            "S": np.array([1 / 2, 0, 0]),
            "T": np.array([0, 0, 1 / 2]),
            "W": np.array([1 / 4, 1 / 4, 1 / 4]),
            "X": np.array([-zeta, zeta, zeta]),
            "X1": np.array([zeta, 1 - zeta, -zeta]),
            "Y": np.array([eta, -eta, eta]),
            "Y1": np.array([1 - eta, eta, -eta]),
            "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
        }

        path = [
            ["Gamma", "X", "L", "T", "W", "R", "X1", "Z", "Gamma", "Y", "S", "W"],
            ["L1", "Y"],
            ["Y1", "Z"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def orcc(self, a, b):
        r"""C-centred orthorhombic (ORCC, oS)"""

        zeta = (1 + a**2 / b**2) / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([zeta, zeta, 1 / 2]),
            "A1": np.array([-zeta, 1 - zeta, 1 / 2]),
            "R": np.array([0, 1 / 2, 1 / 2]),
            "S": np.array([0, 1 / 2, 0]),
            "T": np.array([-1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([zeta, zeta, 0]),
            "X1": np.array([-zeta, 1 - zeta, 0]),
            "Y": np.array([-1 / 2, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "X", "S", "R", "A", "Z", "Gamma", "Y", "X1", "A1", "T", "Y"],
            ["Z", "T"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def hex(self):
        r"""Hexagonal (HEX, hP)"""

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([0, 0, 1 / 2]),
            "H": np.array([1 / 3, 1 / 3, 1 / 2]),
            "K": np.array([1 / 3, 1 / 3, 0]),
            "L": np.array([1 / 2, 0, 1 / 2]),
            "M": np.array([1 / 2, 0, 0]),
        }

        path = [
            ["Gamma", "M", "K", "Gamma", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def rhl1(self, alpha):
        r"""Rhombohedral (RHL, hR),

        .. math::
            \alpha < 90^{\circ}
        """

        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3 / 4 - eta / 2

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "B": np.array([eta, 1 / 2, 1 - eta]),
            "B1": np.array([1 / 2, 1 - eta, eta - 1]),
            "F": np.array([1 / 2, 1 / 2, 0]),
            "L": np.array([1 / 2, 0, 0]),
            "L1": np.array([0, 0, -1 / 2]),
            "P": np.array([eta, nu, nu]),
            "P1": np.array([1 - nu, 1 - nu, 1 - eta]),
            "P2": np.array([nu, nu, eta - 1]),
            "Q": np.array([1 - nu, nu, 0]),
            "X": np.array([nu, 0, -nu]),
            "Z": np.array([1 / 2, 1 / 2, 1 / 2]),
        }

        path = [
            ["Gamma", "L", "B1"],
            ["B", "Z", "Gamma", "X"],
            ["Q", "F", "P1", "Z"],
            ["L", "P"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def rhl2(self, alpha):
        r"""Rhombohedral (RHL, hR),

        .. math::
            \alpha > 90^{\circ}
        """

        eta = 1 / (2 * tan(alpha / 2) ** 2)
        nu = 3 / 4 - eta / 2

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "F": np.array([1 / 2, -1 / 2, 0]),
            "L": np.array([1 / 2, 0, 0]),
            "P": np.array([1 - nu, -nu, 1 - nu]),
            "P1": np.array([nu, nu - 1, nu - 1]),
            "Q": np.array([eta, eta, eta]),
            "Q1": np.array([1 - eta, -eta, -eta]),
            "Z": np.array([1 / 2, -1 / 2, 1 / 2]),
        }

        path = [["Gamma", "P", "Z", "Q", "Gamma", "F", "P1", "Q1", "L", "Z"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mcl(self, b, c, alpha):
        r"""Monoclinic (MCL, mP)"""

        eta = (1 - b * cos(alpha) / c) / (2 * sin(alpha) ** 2)
        nu = 1 / 2 - eta * c * cos(alpha) / b

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2, 0]),
            "C": np.array([0, 1 / 2, 1 / 2]),
            "D": np.array([1 / 2, 0, 1 / 2]),
            "D1": np.array([1 / 2, 0, -1 / 2]),
            "E": np.array([1 / 2, 1 / 2, 1 / 2]),
            "H": np.array([0, eta, 1 - nu]),
            "H1": np.array([0, 1 - eta, nu]),
            "H2": np.array([0, eta, -nu]),
            "M": np.array([1 / 2, eta, 1 - nu]),
            "M1": np.array([1 / 2, 1 - eta, nu]),
            "M2": np.array([1 / 2, eta, -nu]),
            "X": np.array([0, 1 / 2, 0]),
            "Y": np.array([0, 0, 1 / 2]),
            "Y1": np.array([0, 0, -1 / 2]),
            "Z": np.array([1 / 2, 0, 0]),
        }

        path = [
            ["Gamma", "Y", "H", "C", "E", "M1", "A", "X", "H1"],
            ["M", "D", "Z"],
            ["Y", "D"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mclc1(self, a, b, c, alpha):
        r"""C-centred monoclinic (MCLC, mS),

        .. math::
            k_{\gamma} > 90^{\circ}
        """

        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 1 / 2 + 2 * eta * c * cos(alpha) / b
        psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha) ** 2)
        phi = psi + (3 / 4 - psi) * b * cos(alpha) / c

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F1": np.array([zeta, zeta, eta]),
            "F2": np.array([-zeta, -zeta, 1 - eta]),
            "F3": np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 1 / 2]),
            "I1": np.array([1 - phi, phi - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 - psi, psi - 1, 0]),
            "X1": np.array([psi, 1 - psi, 0]),
            "X2": np.array([psi - 1, -psi, 0]),
            "Y": np.array([1 / 2, 1 / 2, 0]),
            "Y1": np.array([-1 / 2, -1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "Y", "F", "L", "I"],
            ["I1", "Z", "F1"],
            ["Y", "X1"],
            ["X", "Gamma", "N"],
            ["M", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mclc2(self, a, b, c, alpha):
        r"""C-centred monoclinic (MCLC, mS),

        .. math::
            k_{\gamma} = 90^{\circ}
        """

        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 1 / 2 + 2 * eta * c * cos(alpha) / b
        psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha) ** 2)
        phi = psi + (3 / 4 - psi) * b * cos(alpha) / c

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F1": np.array([zeta, zeta, eta]),
            "F2": np.array([-zeta, -zeta, 1 - eta]),
            "F3": np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 1 / 2]),
            "I1": np.array([1 - phi, phi - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 - psi, psi - 1, 0]),
            "X1": np.array([psi, 1 - psi, 0]),
            "X2": np.array([psi - 1, -psi, 0]),
            "Y": np.array([1 / 2, 1 / 2, 0]),
            "Y1": np.array([-1 / 2, -1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [["Gamma", "Y", "F", "L", "I"], ["I1", "Z", "F1"], ["N", "Gamma", "M"]]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mclc3(self, a, b, c, alpha):
        r"""C-centred monoclinic (MCLC, mS),

        .. math::
            k_{\gamma} < 90^{\circ},
            \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1
        """

        mu = (1 + b**2 / a**2) / 4
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 1 / 4 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F1": np.array([phi, phi - 1, psi]),
            "F2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([1 / 2, -1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "Y", "F", "H", "Z", "I", "F1"],
            ["H1", "Y1", "X", "Gamma", "N"],
            ["M", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mclc4(self, a, b, c, alpha):
        r"""C-centred monoclinic (MCLC, mS),

        .. math::
            k_{\gamma} < 90^{\circ},
            \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1
        """

        mu = (1 + b**2 / a**2) / 4
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 1 / 4 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F1": np.array([phi, phi - 1, psi]),
            "F2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([1 / 2, -1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "Y", "F", "H", "Z", "I"],
            ["H1", "Y1", "X", "Gamma", "N"],
            ["M", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def mclc5(self, a, b, c, alpha):
        r"""C-centred monoclinic (MCLC, mS),

        .. math::
            k_{\gamma} < 90^{\circ},
            \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1
        """

        zeta = (b**2 / a**2 + (1 - b * cos(alpha) / c) / sin(alpha) ** 2) / 4
        mu = zeta / 2 + b**2 / (4 * a**2) - b * c * cos(alpha) / (2 * a**2)
        nu = 2 * mu - zeta
        rho = 1 - zeta * a**2 / b**2
        omega = (
            (4 * nu - 1 - b**2 * sin(alpha) ** 2 / a**2) * c / (2 * b * cos(alpha))
        )
        eta = 1 / 2 + 2 * zeta * c * cos(alpha) / b
        delta = zeta * c * cos(alpha) / b + omega / 2 - 1 / 4

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "F": np.array([nu, nu, omega]),
            "F1": np.array([1 - nu, 1 - nu, 1 - omega]),
            "F2": np.array([nu, nu - 1, omega]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([rho, 1 - rho, 1 / 2]),
            "I1": np.array([1 - rho, rho - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["Gamma", "Y", "F", "L", "I"],
            ["I1", "Z", "H", "F1"],
            ["H1", "Y1", "X", "Gamma", "N"],
            ["M", "Gamma"],
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def tri1a(self, a):
        r"""Triclinic (TRI, aP),

        .. math::
            k_{\alpha} > 90^{\circ},
            k_{\beta} > 90^{\circ},
            k_{\gamma} > 90^{\circ}
        .. math::
            k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})
        """

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "L": np.array([1 / 2, 1 / 2, 0]),
            "M": np.array([0, 1 / 2, 1 / 2]),
            "N": np.array([1 / 2, 0, 1 / 2]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["X", "Gamma", "Y"],
            ["L", "Gamma", "Z"],
            ["N", "Gamma", "M"],
            "R",
            "Gamma",
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def tri1b(self, a):
        r"""Triclinic (TRI, aP),

        .. math::
            k_{\alpha} < 90^{\circ},
            k_{\beta} < 90^{\circ},
            k_{\gamma} < 90^{\circ}
        .. math::
            k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})
        """

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "L": np.array([1 / 2, -1 / 2, 0]),
            "M": np.array([0, 0, 1 / 2]),
            "N": np.array([-1 / 2, -1 / 2, 1 / 2]),
            "R": np.array([0, -1 / 2, 1 / 2]),
            "X": np.array([0, -1 / 2, 0]),
            "Y": np.array([1 / 2, 0, 0]),
            "Z": np.array([-1 / 2, 0, 1 / 2]),
        }

        path = [
            ["X", "Gamma", "Y"],
            ["L", "Gamma", "Z"],
            ["N", "Gamma", "M"],
            "R",
            "Gamma",
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def tri2a(self, a):
        r"""Triclinic (TRI, aP),

        .. math::
            k_{\alpha} > 90^{\circ},
            k_{\beta} > 90^{\circ},
            k_{\gamma} = 90^{\circ}
        """

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "L": np.array([1 / 2, 1 / 2, 0]),
            "M": np.array([0, 1 / 2, 1 / 2]),
            "N": np.array([1 / 2, 0, 1 / 2]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        path = [
            ["X", "Gamma", "Y"],
            ["L", "Gamma", "Z"],
            ["N", "Gamma", "M"],
            "R",
            "Gamma",
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)

    def tri2b(self, a):
        r"""Triclinic (TRI, aP),

        .. math::
            k_{\alpha} < 90^{\circ},
            k_{\beta} < 90^{\circ},
            k_{\gamma} = 90^{\circ}
        """

        kpoints = {
            "Gamma": np.array([0, 0, 0]),
            "L": np.array([1 / 2, -1 / 2, 0]),
            "M": np.array([0, 0, 1 / 2]),
            "N": np.array([-1 / 2, -1 / 2, 1 / 2]),
            "R": np.array([0, -1 / 2, 1 / 2]),
            "X": np.array([0, -1 / 2, 0]),
            "Y": np.array([1 / 2, 0, 0]),
            "Z": np.array([-1 / 2, 0, 1 / 2]),
        }

        path = [
            ["X", "Gamma", "Y"],
            ["L", "Gamma", "Z"],
            ["N", "Gamma", "M"],
            "R",
            "Gamma",
        ]

        self._kpoints = kpoints
        self.path = path

        return HighSymmetryPoints(kpoints, path)
