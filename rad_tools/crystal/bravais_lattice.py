r"""14 Bravais lattice"""

from math import sqrt, sin, cos, tan, pi
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

__all__ = ["BravaisLattice", "CUB"]

TOLERANCE = 10e-6

RED = "#FF4D67"
GREEN = "#58EC2E"
ORANGE = "#F7CB3D"
BLUE = "#274DD1"
PURPLE = "#DC5CFF"


def deduct_zone(
    cell,
):
    r"""
    Construct Brillouin zone (or Wigner-Seitz cell), from the given set of basis vectors.

    cell : (3,3) array_like
        Matrix of the basis vectors, rows are interpreted as vectors,
        columns as cartesian coordinates:

        .. code-block:: python
            cell = [
                [b1_x, b1_y, b1_z],
                [b2_x, b2_y, b2_z],
                [b3_x, b3_y, b3_z]
                ]
    Returns
    -------
    planes: (N, 3) array
        Array of N planes, constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v`, which is perpendicular to the plane
        and gives the coordinate of the point from the plane.
        The vector is given in relative coordinates, with respect to the basis vectors
    edges : (K, 2, 3) array
        Array of K edges of the Brillouin zone (or Wigner-Seitz cell).
        Define by the two points.
    corners: (M, 3) array
        Array of M corners of the Brillouin zone (or Wigner-Seitz cell).
    """

    cell = np.array(cell)

    # Compute planes
    planes_candidates = []
    other_points = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if i**2 + j**2 + k**2 != 0:
                    rel_vector = np.array((i, j, k))
                    vector = np.matmul(rel_vector, cell)
                    distance_gamma = np.linalg.norm(vector) / 2
                    planes_candidates.append((rel_vector, vector, distance_gamma))
                    other_points.append(vector)

    planes_candidates.sort(key=lambda x: x[2])
    planes = []
    for plane in planes_candidates:
        takeit = True
        for other_point in other_points:
            if (
                not (other_point == plane[1]).all()
                and (plane[2] - np.linalg.norm(plane[1] / 2 - other_point)) > -TOLERANCE
            ):
                takeit = False
                break
        if takeit:
            planes.append(plane[0])

    # Compute corners
    corners_candidates = []
    for f, fplane in enumerate(planes):
        nf = np.matmul(fplane, cell) / 2
        for s in range(f + 1, len(planes)):
            splane = planes[s]
            ns = np.matmul(splane, cell) / 2
            for t in range(s + 1, len(planes)):
                tplane = planes[t]
                nt = np.matmul(tplane, cell) / 2
                A = np.array([nf, ns, nt])
                b = np.array(
                    [
                        np.linalg.norm(nf) ** 2,
                        np.linalg.norm(ns) ** 2,
                        np.linalg.norm(nt) ** 2,
                    ]
                )
                try:
                    x = np.linalg.solve(A, b)
                    corners_candidates.append(x)
                except:
                    pass

    corners = []

    # Check if corner is closer to the Gamma point then to any other point of lattice.
    for corner in corners_candidates:
        takeit = True
        for other_point in other_points:
            if (
                np.linalg.norm(corner) - np.linalg.norm(corner - other_point)
            ) > TOLERANCE:
                takeit = False
                break
        if takeit:
            corners.append(corner)

    # Filter equal entries
    new_corners = []
    for i in range(len(corners)):
        takeit = True
        for j in range(len(new_corners)):
            if (np.abs(corners[i] - new_corners[j]) < TOLERANCE).all():
                takeit = False
                break
        if takeit:
            new_corners.append(corners[i])
    corners = new_corners

    # Compute corners
    corners_planes = []
    for corner in corners:
        corners_planes.append([])
        for plane in planes:
            A, B, C = tuple(np.matmul(plane, cell) / 2)
            D = -np.linalg.norm(np.matmul(plane, cell) / 2) ** 2
            if abs(A * corner[0] + B * corner[1] + C * corner[2] + D) < TOLERANCE:
                corners_planes[-1].append(plane)
    edges = []
    for fc, fcorner in enumerate(corners):
        for sc in range(fc + 1, len(corners)):
            scorner = corners[sc]
            n = 0
            for fplane in corners_planes[fc]:
                for splane in corners_planes[sc]:
                    if np.linalg.norm(fplane - splane) < TOLERANCE:
                        n += 1
            if n == 2:
                edges.append([fcorner, scorner])
    return planes, edges, corners


class Lattice:
    PLOT_NAMES = {
        "G": "$\\Gamma$",
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
        "Y": "$Y$",
        "Y1": "$Y_1$",
        "S": "$S$",  # it is overwritten to sigma if needed.
        "S1": "$S_1$",  # it is overwritten to sigma if needed.
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

    def __init__(self, a1, a2, a3) -> None:
        self.cell = np.array([a1, a2, a3])
        self.primitive = False
        self.points = {}
        self.path = []

    @property
    def a1(self):
        return self.cell[0]

    @property
    def a2(self):
        return self.cell[1]

    @property
    def a3(self):
        return self.cell[2]

    @property
    def unit_cell_volume(self):
        r"""
        Volume of the unit cell.

        .. math::

            V = \delta_{\vec{A}}\cdot(\delta_{\vec{B}}\times\delta_{\vec{C}})
        """

        return np.dot(self.a1, np.cross(self.a2, self.a3))

    @property
    def reciprocal_cell(self):
        # TODO Rethink the issue with the primitive-non primitive
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True

        result = np.array(
            [
                2 * pi / self.unit_cell_volume * np.cross(self.a2, self.a3),
                2 * pi / self.unit_cell_volume * np.cross(self.a3, self.a1),
                2 * pi / self.unit_cell_volume * np.cross(self.a1, self.a2),
            ]
        )
        if reverse:
            self.make_conventional()
        return result

    @property
    def b1(self):
        r"""
        First reciprocal lattice vector.

        .. math::

            \vec{b}_1 = \frac{2\pi}{V}\vec{b}\times\vec{c}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[0]

    @property
    def b2(self):
        r"""
        Second reciprocal lattice vector.

        .. math::

            \vec{b}_2 = \frac{2\pi}{V}\vec{c}\times\vec{a}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[1]

    @property
    def b3(self):
        r"""
        Third reciprocal lattice vector.

        .. math::

            \vec{b}_3 = \frac{2\pi}{V}\vec{a}\times\vec{b}

        where :math:`V = \vec{a}\cdot\vec{b}\times\vec{c}`
        """

        return self.reciprocal_cell[2]

    @property
    def k_alpha(self):
        v1 = self.b2 / np.linalg.norm(self.b2)
        v2 = self.b3 / np.linalg.norm(self.b3)
        return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) / pi * 180

    @property
    def k_beta(self):
        v1 = self.b1 / np.linalg.norm(self.b1)
        v2 = self.b3 / np.linalg.norm(self.b3)
        return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) / pi * 180

    @property
    def k_gamma(self):
        v1 = self.b1 / np.linalg.norm(self.b1)
        v2 = self.b2 / np.linalg.norm(self.b2)
        return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) / pi * 180

    def make_primitive(self):
        r"""
        This method does nothing.
        Conventional and primitive lattice are the same.
        """
        pass

    def make_conventional(self):
        r"""
        This method does nothing.
        Conventional and primitive lattice are the same.
        """
        pass

    @property
    def variant(self):
        r"""There is no variations of the Lattice"""
        return self.__class__.__name__

    @property
    def pearson_symbol(self):
        raise NotImplementedError

    def prepare_figure(self, background=True, focal_length=0.2) -> None:
        self._fig = plt.figure(figsize=(6, 6))
        rcParams["axes.linewidth"] = 0
        rcParams["xtick.color"] = "#B3B3B3"
        self._ax = self._fig.add_subplot(projection="3d")
        self._ax.set_proj_type("persp", focal_length=focal_length)
        if background:
            self._ax.axes.linewidth = 0
            self._ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            self._ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            self._ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
            self._ax.set_xlabel("x", fontsize=15, alpha=0.5)
            self._ax.set_ylabel("y", fontsize=15, alpha=0.5)
            self._ax.set_zlabel("z", fontsize=15, alpha=0.5)
            self._ax.tick_params(axis="both", zorder=0, color="#B3B3B3")
        else:
            self._ax.axis("off")

    def plot(self, kind="conventional", **kwargs):
        if isinstance(kind, str):
            kinds = [kind]
        else:
            kinds = kind
        try:
            a = self._ax
            del a
        except AttributeError:
            raise ValueError("Use .prepare_figure() first.")
        try:
            for kind in kinds:
                getattr(self, f"plot_{kind}")(self._ax, **kwargs)
        except AttributeError:
            d = dir(self)
            for i in d:
                print(i, i == f"plot_{kind}")
            raise ValueError(f"Plot kind '{kind}' does not exist!")

    def show(self):
        self._ax.set_aspect("equal")
        plt.show()
        del self._fig
        del self._ax
        plt.close()

    def savefig(self, output_name="graph.png", elev=30, azim=-60, **kwargs):
        self._ax.set_aspect("equal")
        self._ax.view_init(elev=elev, azim=azim)
        self._fig.savefig(output_name, **kwargs)

    def legend(self, **kwargs):
        self._ax.legend(**kwargs)

    def plot_real_space(
        self, ax, vectors=True, colour=BLUE, label=None, vector_pad=1.1
    ):
        if label is not None:
            ax.scatter(0, 0, 0, color=colour, label=label)
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            ax.text(
                self.cell[0][0] * vector_pad[0],
                self.cell[0][1] * vector_pad[0],
                self.cell[0][2] * vector_pad[0],
                "$a_1$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.cell[1][0] * vector_pad[2],
                self.cell[1][1] * vector_pad[2],
                self.cell[1][2] * vector_pad[2],
                "$a_2$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.cell[2][0] * vector_pad[2],
                self.cell[2][1] * vector_pad[2],
                self.cell[2][2] * vector_pad[2],
                "$a_3$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            for i in self.cell:
                ax.quiver(
                    0,
                    0,
                    0,
                    *tuple(i),
                    arrow_length_ratio=0.2,
                    color=colour,
                    alpha=0.7,
                    linewidth=2,
                )
                ax.scatter(*tuple(i), s=0)

        def plot_line(line, shift):
            ax.plot(
                [shift[0], shift[0] + line[0]],
                [shift[1], shift[1] + line[1]],
                [shift[2], shift[2] + line[2]],
                color=colour,
            )

        for i in range(0, 3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            plot_line(self.cell[i], np.zeros(3))
            plot_line(self.cell[i], self.cell[j])
            plot_line(self.cell[i], self.cell[k])
            plot_line(self.cell[i], self.cell[j] + self.cell[k])

    def plot_reciprocal_space(
        self, ax, vectors=True, colour=RED, label=None, vector_pad=1.1
    ):
        if label is not None:
            ax.scatter(0, 0, 0, color=colour, label=label)
        planes, edges, corners = deduct_zone([self.b1, self.b2, self.b3])

        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            ax.text(
                self.b1[0] * vector_pad[0],
                self.b1[1] * vector_pad[0],
                self.b1[2] * vector_pad[0],
                "$b_1$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.b2[0] * vector_pad[1],
                self.b2[1] * vector_pad[1],
                self.b2[2] * vector_pad[1],
                "$b_2$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.b3[0] * vector_pad[2],
                self.b3[1] * vector_pad[2],
                self.b3[2] * vector_pad[2],
                "$b_3$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            for i in [self.b1, self.b2, self.b3]:
                ax.quiver(
                    0,
                    0,
                    0,
                    *tuple(i),
                    arrow_length_ratio=0.2,
                    color=colour,
                    alpha=0.7,
                    linewidth=2,
                )
                ax.scatter(*tuple(i), s=0)
        for a, b in edges:
            ax.plot(
                [a[0], b[0]],
                [a[1], b[1]],
                [a[2], b[2]],
                color=colour,
            )

    def plot_conventional(self, ax, **kwargs):
        reverse = False
        if self.primitive:
            self.make_conventional()
            reverse = True
        self.plot_real_space(ax, **kwargs)
        if reverse:
            self.make_primitive()

    def plot_primitive(self, ax, **kwargs):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        self.plot_real_space(ax, **kwargs)
        if reverse:
            self.make_conventional()

    def plot_brillouin(self, ax, **kwargs):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        self.plot_reciprocal_space(ax, **kwargs)
        if reverse:
            self.make_conventional()

    def plot_kpath(self, ax, colour="black", label=None):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        for point in self.points:
            ax.scatter(
                *tuple(self.points[point] @ self.reciprocal_cell),
                s=36,
                color=colour,
            )
            ax.text(
                *tuple(
                    self.points[point] @ self.reciprocal_cell
                    + 0.025 * self.b1
                    + +0.025 * self.b2
                    + 0.025 * self.b3
                ),
                self.PLOT_NAMES[point],
                fontsize=20,
                color=colour,
            )
        if label is not None:
            ax.scatter(
                *tuple(self.points[point] @ self.reciprocal_cell),
                color=colour,
                label=label,
            )

        for subpath in self.path:
            for i in range(len(subpath) - 1):
                ax.plot(
                    *tuple(
                        np.concatenate(
                            (
                                self.points[subpath[i]] @ self.reciprocal_cell,
                                self.points[subpath[i + 1]] @ self.reciprocal_cell,
                            )
                        )
                        .reshape(2, 3)
                        .T
                    ),
                    color=colour,
                    alpha=0.5,
                    linewidth=3,
                )
        if reverse:
            self.make_conventional()

    def plot_brillouin_kpath(self, ax, zone_colour=RED, path_colour="black", **kwargs):
        self.plot_brillouin(ax, colour=zone_colour, **kwargs)
        self.plot_kpath(ax, colour=path_colour, **kwargs)

    def plot_wigner_seitz(
        self, ax, vectors=True, colour="black", label=None, vector_pad=1.1
    ):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        planes, edges, corners = deduct_zone([self.a1, self.a2, self.a3])

        if label is not None:
            ax.scatter(0, 0, 0, color=colour, label=label)
        if vectors:
            if not isinstance(vector_pad, Iterable):
                vector_pad = [vector_pad, vector_pad, vector_pad]
            ax.text(
                self.a1[0] * vector_pad[0],
                self.a1[1] * vector_pad[0],
                self.a1[2] * vector_pad[0],
                "$a_1$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.a2[0] * vector_pad[1],
                self.a2[1] * vector_pad[1],
                self.a2[2] * vector_pad[1],
                "$a_2$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            ax.text(
                self.a3[0] * vector_pad[2],
                self.a3[1] * vector_pad[2],
                self.a3[2] * vector_pad[2],
                "$a_3$",
                fontsize=20,
                color=colour,
                ha="center",
                va="center",
            )
            for i in [self.a1, self.a2, self.a3]:
                ax.quiver(
                    0, 0, 0, *tuple(i), arrow_length_ratio=0.2, color=colour, alpha=0.5
                )
                ax.scatter(*tuple(i), s=0)
        for a, b in edges:
            ax.plot(
                [a[0], b[0]],
                [a[1], b[1]],
                [a[2], b[2]],
                color=colour,
            )

        if reverse:
            self.make_conventional()


# 1
class CUB(Lattice):
    r"""Cubic (CUB, cP)"""

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive = True
        self.points = {
            "G": np.array([0, 0, 0]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
        }
        self.path = [["G", "X", "M", "G", "R", "X"], ["M", "R"]]


# 2
class FCC(Lattice):
    r"""Face-centred cubic (FCC, cF)"""

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive = False
        self.points = {
            "G": np.array([0, 0, 0]),
            "K": np.array([3 / 8, 3 / 8, 3 / 4]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "U": np.array([5 / 8, 1 / 4, 5 / 8]),
            "W": np.array([1 / 2, 1 / 4, 3 / 4]),
            "X": np.array([1 / 2, 0, 1 / 2]),
        }
        self.path = [
            ["G", "X", "W", "K", "G", "L", "U", "W", "L", "K"],
            ["U", "X"],
        ]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.a, self.a])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array([[0, self.a, self.a], [self.a, 0, self.a], [self.a, self.a, 0]])
            / 2
        )
        self.primitive = True


# 3
class BCC(Lattice):
    r"""Body-centered cubic (BCC, cl)"""

    def __init__(self, a: float) -> None:
        self.a = a
        self.cell = np.diag([a, a, a])
        self.primitive = False
        self.points = {
            "G": np.array([0, 0, 0]),
            "H": np.array([1 / 2, -1 / 2, 1 / 2]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "N": np.array([0, 0, 1 / 2]),
        }
        self.path = [["G", "H", "N", "G", "P", "H"], ["P", "N"]]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.a, self.a])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array(
                [
                    [-self.a, self.a, self.a],
                    [self.a, -self.a, self.a],
                    [self.a, self.a, -self.a],
                ]
            )
            / 2
        )
        self.primitive = True


# 4
class TET(Lattice):
    r"""Tetragonal (TET, tP)"""

    def __init__(self, a: float, c: float) -> None:
        self.a = a
        self.c = c
        self.cell = np.diag([a, a, c])
        self.primitive = True
        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 1 / 2, 0]),
            "R": np.array([0, 1 / 2, 1 / 2]),
            "X": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }
        self.path = [
            ["G", "X", "M", "G", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"],
        ]


# 5
class BCT(Lattice):
    r"""Body-centred tetragonal (BCT, tI)"""

    def __init__(self, a: float, c: float) -> None:
        self.PLOT_NAMES["S"] = "$\\Sigma$"
        self.PLOT_NAMES["S1"] = "$\\Sigma_1$"
        if a == c:
            raise ValueError("Are you trying to create BCC Lattice (a == c)?")
        self.a = a
        self.c = c
        self.cell = np.diag([a, a, c])
        self.primitive = False
        if self.variant == "BCT1":
            eta = (1 + self.c**2 / self.a**2) / 4
            self.points = {
                "G": np.array([0, 0, 0]),
                "M": np.array([-1 / 2, 1 / 2, 1 / 2]),
                "N": np.array([0, 1 / 2, 0]),
                "P": np.array([1 / 4, 1 / 4, 1 / 4]),
                "X": np.array([0, 0, 1 / 2]),
                "Z": np.array([eta, eta, -eta]),
                "Z1": np.array([-eta, 1 - eta, eta]),
            }

            self.path = [
                ["G", "X", "M", "G", "Z", "P", "N", "Z1", "M"],
                ["X", "P"],
            ]

        elif self.variant == "BCT2":
            eta = (1 + self.a**2 / self.c**2) / 4
            zeta = self.a**2 / (2 * self.c**2)
            self.points = {
                "G": np.array([0, 0, 0]),
                "N": np.array([0, 1 / 2, 0]),
                "P": np.array([1 / 4, 1 / 4, 1 / 4]),
                "S": np.array([-eta, eta, eta]),
                "S1": np.array([eta, 1 - eta, -eta]),
                "X": np.array([0, 0, 1 / 2]),
                "Y": np.array([-zeta, zeta, 1 / 2]),
                "Y1": np.array([1 / 2, 1 / 2, -zeta]),
                "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
            }

            self.path = [
                [
                    "G",
                    "X",
                    "Y",
                    "S",
                    "G",
                    "Z",
                    "S1",
                    "N",
                    "P",
                    "Y1",
                    "Z",
                ],
                ["X", "P"],
            ]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.a, self.c])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array(
                [
                    [-self.a, self.a, self.c],
                    [self.a, -self.a, self.c],
                    [self.a, self.a, -self.c],
                ]
            )
            / 2
        )
        self.primitive = True

    @property
    def variant(self):
        r"""
        Two variants of the Lattice.

        :math:`\text{BCT}_1: c < a` and :math:`\text{BCT}_2: c > a`
        """
        if self.a > self.c:
            return "BCT1"
        elif self.a < self.c:
            return "BCT2"


# 6
class ORC(Lattice):
    r"""
    Orthorhombic (ORC, oP)

    :math:`a < b < c`
    """

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct CUB Lattice (a = b = c)?")
        if a == b or b == c:
            raise ValueError(
                "Are you trying to construct TET Lattice (a = b != c or a != b == c)?"
            )

        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive = True
        self.points = {
            "G": np.array([0, 0, 0]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "S": np.array([1 / 2, 1 / 2, 0]),
            "T": np.array([0, 1 / 2, 1 / 2]),
            "U": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

        self.path = [
            ["G", "X", "S", "Y", "G", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"],
        ]


# 7
class ORCF(Lattice):
    r"""
    Face-centred orthorhombic (ORCF, oF)

    :math:`a < b < c`
    """

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct FCC Lattice (a = b = c)?")
        if a == b:
            raise ValueError("FIXME, dont know which lattice it will be.")
        if b == c:
            raise ValueError("FIXME, dont know which lattice it will be.")
        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive = False
        if self.variant == "ORCF1":
            eta = (1 + self.a**2 / self.b**2 + self.a**2 / self.c**2) / 4
            zeta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4
            self.points = {
                "G": np.array([0, 0, 0]),
                "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
                "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "T": np.array([1, 1 / 2, 1 / 2]),
                "X": np.array([0, eta, eta]),
                "X1": np.array([1, 1 - eta, 1 - eta]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self.path = [
                ["G", "Y", "T", "Z", "G", "X", "A1", "Y"],
                ["T", "X1"],
                ["X", "A", "Z"],
                ["L", "G"],
            ]
        elif self.variant == "ORCF2":
            eta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4
            delta = (1 + self.b**2 / self.a**2 - self.b**2 / self.c**2) / 4
            phi = (1 + self.c**2 / self.b**2 - self.c**2 / self.a**2) / 4

            self.points = {
                "G": np.array([0, 0, 0]),
                "C": np.array([1 / 2, 1 / 2 - eta, 1 - eta]),
                "C1": np.array([1 / 2, 1 / 2 + eta, eta]),
                "D": np.array([1 / 2 - delta, 1 / 2, 1 - delta]),
                "D1": np.array([1 / 2 + delta, 1 / 2, delta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "H": np.array([1 - phi, 1 / 2 - phi, 1 / 2]),
                "H1": np.array([phi, 1 / 2 + phi, 1 / 2]),
                "X": np.array([0, 1 / 2, 1 / 2]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self.path = [
                ["G", "Y", "C", "D", "X", "G", "Z", "D1", "H", "C"],
                ["C1", "Z"],
                ["X", "H1"],
                ["H", "Y"],
                ["L", "G"],
            ]
        elif self.variant == "ORCF3":
            eta = (1 + self.a**2 / self.b**2 + self.a**2 / self.c**2) / 4
            zeta = (1 + self.a**2 / self.b**2 - self.a**2 / self.c**2) / 4

            self.points = {
                "G": np.array([0, 0, 0]),
                "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
                "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
                "L": np.array([1 / 2, 1 / 2, 1 / 2]),
                "T": np.array([1, 1 / 2, 1 / 2]),
                "X": np.array([0, eta, eta]),
                "Y": np.array([1 / 2, 0, 1 / 2]),
                "Z": np.array([1 / 2, 1 / 2, 0]),
            }

            self.path = [
                ["G", "Y", "T", "Z", "G", "X", "A1", "Y"],
                ["X", "A", "Z"],
                ["L", "G"],
            ]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.b, self.c])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array(
                [
                    [0, self.b, self.c],
                    [self.a, 0, self.c],
                    [self.a, self.b, 0],
                ]
            )
            / 2
        )
        self.primitive = True

    @property
    def variant(self):
        r"""
        Three variants of the Lattice.

        :math:`\text{ORCF}_1: \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_2: \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        :math:`\text{ORCF}_3: \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
        """
        if 1 / self.a**2 > 1 / self.b**2 + 1 / self.c**2:
            return "ORCF1"
        elif 1 / self.a**2 < 1 / self.b**2 + 1 / self.c**2:
            return "ORCF2"
        else:
            return "ORCF3"


# 8
class ORCI(Lattice):
    r"""
    Body-centred orthorhombic (ORCI, oI)

    :math:`a < b < c`
    """

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct BCC Lattice (a = b = c)?")
        if a == b:
            raise ValueError("Are you trying to construct BCT2 Lattice (a = b < c)?")
        if b == c:
            raise ValueError("Are you trying to construct BCT1 Lattice (a < b = c)?")
        self.a = a
        self.b = b
        self.c = c
        self.cell = np.diag([a, b, c])
        self.primitive = False
        zeta = (1 + self.a**2 / self.c**2) / 4
        eta = (1 + self.b**2 / self.c**2) / 4
        delta = (self.b**2 - a**2) / (4 * self.c**2)
        mu = (self.a**2 + self.b**2) / (4 * self.c**2)

        self.points = {
            "G": np.array([0, 0, 0]),
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

        self.path = [
            ["G", "X", "L", "T", "W", "R", "X1", "Z", "G", "Y", "S", "W"],
            ["L1", "Y"],
            ["Y1", "Z"],
        ]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.b, self.c])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array(
                [
                    [-self.a, self.b, self.c],
                    [self.a, -self.b, self.c],
                    [self.a, self.b, -self.c],
                ]
            )
            / 2
        )
        self.primitive = True


# 9
class ORCC(Lattice):
    r"""
    C-centred orthorhombic (ORCC, oS)

    :math:`a < b`
    """

    def __init__(self, a: float, b: float, c: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if a == b == c:
            raise ValueError("Are you trying to construct TET Lattice (a = b = c)?")
        self.a = a
        self.b = c
        self.c = b
        self.cell = np.diag([a, b, c])
        self.primitive = False
        zeta = (1 + self.a**2 / self.b**2) / 4

        self.points = {
            "G": np.array([0, 0, 0]),
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

        self.path = [
            ["G", "X", "S", "R", "A", "Z", "G", "Y", "X1", "A1", "T", "Y"],
            ["Z", "T"],
        ]

    def make_conventional(self):
        self.cell = np.diag([self.a, self.b, self.c])
        self.primitive = False

    def make_primitive(self):
        self.cell = (
            np.array(
                [
                    [self.a, -self.b, 0],
                    [self.a, self.b, 0],
                    [0, 0, 2 * self.c],
                ]
            )
            / 2
        )
        self.primitive = True


# 10
class HEX(Lattice):
    r"""Hexagonal (HEX, hP)"""

    def __init__(self, a: float, c: float) -> None:
        self.a = a
        self.c = c
        self.cell = np.array(
            [[a / 2, -a * sqrt(3) / 2, 0], [a / 2, a * sqrt(3) / 2, 0], [0, 0, c]]
        )
        self.primitive = True

        self.points = {
            "G": np.array([0, 0, 0]),
            "A": np.array([0, 0, 1 / 2]),
            "H": np.array([1 / 3, 1 / 3, 1 / 2]),
            "K": np.array([1 / 3, 1 / 3, 0]),
            "L": np.array([1 / 2, 0, 1 / 2]),
            "M": np.array([1 / 2, 0, 0]),
        }

        self.path = [
            ["G", "M", "K", "G", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]


# 11
class RHL(Lattice):
    r"""Rhombohedral (RHL, hR)"""

    def __init__(self, a: float, alpha: float) -> None:
        if alpha == 90:
            raise ValueError("Are you trying to construct CUB Lattice (alpha == 90)?")
        if alpha >= 120:
            raise ValueError("alpha has to be < 120 degrees.")
        self.a = a
        self.alpha = alpha
        self.cell = np.array(
            [
                [a * cos(alpha / 180 * pi / 2), -a * sin(alpha / 180 * pi / 2), 0],
                [a * cos(alpha / 180 * pi / 2), a * sin(alpha / 180 * pi / 2), 0],
                [
                    a * cos(alpha / 180 * pi) / cos(alpha / 180 * pi / 2),
                    0,
                    a
                    * sqrt(
                        1 - cos(alpha / 180 * pi) ** 2 / cos(alpha / 180 * pi / 2) ** 2
                    ),
                ],
            ]
        )
        self.primitive = True
        if self.variant == "RHL1":
            eta = (1 + 4 * cos(alpha / 180 * pi)) / (2 + 4 * cos(alpha / 180 * pi))
            nu = 3 / 4 - eta / 2

            self.points = {
                "G": np.array([0, 0, 0]),
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

            self.path = [
                ["G", "L", "B1"],
                ["B", "Z", "G", "X"],
                ["Q", "F", "P1", "Z"],
                ["L", "P"],
            ]
        elif self.variant == "RHL2":
            eta = 1 / (2 * tan(self.alpha / 180 * pi / 2) ** 2)
            nu = 3 / 4 - eta / 2

            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([1 / 2, -1 / 2, 0]),
                "L": np.array([1 / 2, 0, 0]),
                "P": np.array([1 - nu, -nu, 1 - nu]),
                "P1": np.array([nu, nu - 1, nu - 1]),
                "Q": np.array([eta, eta, eta]),
                "Q1": np.array([1 - eta, -eta, -eta]),
                "Z": np.array([1 / 2, -1 / 2, 1 / 2]),
            }

            self.path = [["G", "P", "Z", "Q", "G", "F", "P1", "Q1", "L", "Z"]]

    @property
    def variant(self):
        r"""
        Two variants of the Lattice.

        :math:`\text{RHL}_1 \alpha < 90^{\circ}`,
        :math:`\text{RHL}_2 \alpha > 90^{\circ}`
        """
        if self.alpha < 90:
            return "RHL1"
        elif self.alpha > 90:
            return "RHL2"


# 12
class MCL(Lattice):
    r"""
    Monoclinic (MCL, mP)

    :math:`a, b \le c`, :math:`\alpha < 90^{\circ}`, :math:`\beta = \gamma = 90^{\circ}`.
    """

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        a, b, c = tuple(sorted([a, b, c]))
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.cell = np.array(
            [
                [self.a, 0, 0],
                [0, self.b, 0],
                [0, c * cos(self.alpha / 180 * pi), c * sin(self.alpha / 180 * pi)],
            ]
        )
        self.primitive = True

        eta = (1 - self.b * cos(self.alpha / 180 * pi) / self.c) / (
            2 * sin(self.alpha / 180 * pi) ** 2
        )
        nu = 1 / 2 - eta * self.c * cos(self.alpha / 180 * pi) / self.b

        self.points = {
            "G": np.array([0, 0, 0]),
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

        self.path = [
            ["G", "Y", "H", "C", "E", "M1", "A", "X", "H1"],
            ["M", "D", "Z"],
            ["Y", "D"],
        ]


# 13
class MCLC(Lattice):
    r"""
    C-centred monoclinic (MCLC, mS)

    :math:`a, b \le c`, :math:`\alpha < 90^{\circ}`, :math:`\beta = \gamma = 90^{\circ}`.
    """

    def __init__(self, a: float, b: float, c: float, alpha: float) -> None:
        if a > c:
            a, c = c, a
        if b > c:
            b, c = c, b
        if alpha > 90:
            raise ValueError("alpha has to be < 90")
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.cell = None
        self.make_conventional()
        self.primitive = False
        # Parameters
        if self.variant in ["MCLC1", "MCLC2"]:
            zeta = (2 - b * cos(alpha / 180 * pi) / c) / (
                4 * sin(alpha / 180 * pi) ** 2
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            psi = 3 / 4 - a**2 / (4 * b**2 * sin(alpha / 180 * pi) ** 2)
            phi = psi + (3 / 4 - psi) * b * cos(alpha / 180 * pi) / c
        elif self.variant in ["MCLC3", "MCLC4"]:
            mu = (1 + b**2 / a**2) / 4
            delta = b * c * cos(alpha / 180 * pi) / (2 * a**2)
            zeta = (
                mu
                - 1 / 4
                + (1 - b * cos(alpha / 180 * pi) / c) / (4 * sin(alpha / 180 * pi) ** 2)
            )
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            phi = 1 + zeta - 2 * mu
            psi = eta - 2 * delta
        elif self.variant == "MCLC5":
            zeta = (
                b**2 / a**2
                + (1 - b * cos(alpha / 180 * pi) / c) / sin(alpha / 180 * pi) ** 2
            ) / 4
            eta = 1 / 2 + 2 * zeta * c * cos(alpha / 180 * pi) / b
            mu = (
                eta / 2
                + b**2 / (4 * a**2)
                - b * c * cos(alpha / 180 * pi) / (2 * a**2)
            )
            nu = 2 * mu - zeta
            rho = 1 - zeta * a**2 / b**2
            omega = (
                (4 * nu - 1 - b**2 * sin(alpha / 180 * pi) ** 2 / a**2)
                * c
                / (2 * b * cos(alpha / 180 * pi))
            )
            delta = zeta * c * cos(alpha / 180 * pi) / b + omega / 2 - 1 / 4

        # Path
        if self.variant == "MCLC1":
            self.path = [
                ["G", "Y", "F", "L", "I"],
                ["I1", "Z", "F1"],
                ["Y", "X1"],
                ["X", "G", "N"],
                ["M", "G"],
            ]
            self.points = {
                "G": np.array([0, 0, 0]),
                "N": np.array([1 / 2, 0, 0]),
                "N1": np.array([0, -1 / 2, 0]),
                "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
                "F1": np.array([zeta, zeta, eta]),
                "F2": np.array([-zeta, -zeta, 1 - eta]),
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
        elif self.variant == "MCLC2":
            self.path = [["G", "Y", "F", "L", "I"], ["I1", "Z", "F1"], ["N", "G", "M"]]
            self.points = {
                "G": np.array([0, 0, 0]),
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
                "Y": np.array([1 / 2, 1 / 2, 0]),
                "Y1": np.array([-1 / 2, -1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }
        elif self.variant == "MCLC3":
            self.points = {
                "G": np.array([0, 0, 0]),
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
            self.path = [
                ["G", "Y", "F", "H", "Z", "I", "F1"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
        elif self.variant == "MCLC4":
            self.points = {
                "G": np.array([0, 0, 0]),
                "F": np.array([1 - phi, 1 - phi, 1 - psi]),
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
            self.path = [
                ["G", "Y", "F", "H", "Z", "I"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
        elif self.variant == "MCLC5":
            self.path = [
                ["G", "Y", "F", "L", "I"],
                ["I1", "Z", "H", "F1"],
                ["H1", "Y1", "X", "G", "N"],
                ["M", "G"],
            ]
            self.points = {
                "G": np.array([0, 0, 0]),
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

    def make_conventional(self):
        self.cell = np.array(
            [
                [self.a, 0, 0],
                [0, self.b, 0],
                [
                    0,
                    self.c * cos(self.alpha / 180 * pi),
                    self.c * sin(self.alpha / 180 * pi),
                ],
            ]
        )
        self.primitive = False

    def make_primitive(self):
        self.cell = np.array(
            [
                [self.a / 2, self.b / 2, 0],
                [-self.a / 2, self.b / 2, 0],
                [
                    0,
                    self.c * cos(self.alpha / 180 * pi),
                    self.c * sin(self.alpha / 180 * pi),
                ],
            ]
        )
        self.primitive = True

    @property
    def variant(self):
        r"""
        Five variant of the Lattice.

        :math:`\text{MCLC}_1: k_{\gamma} > 90^{\circ}`,
        :math:`\text{MCLC}_2: k_{\gamma} = 90^{\circ}`,
        :math:`\text{MCLC}_3: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`
        :math:`\text{MCLC}_4: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`
        :math:`\text{MCLC}_5: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`
        """

        if abs(self.k_gamma - 90) < 10e-5:
            return "MCLC2"
        elif self.k_gamma > 90:
            return "MCLC1"
        # TODO think about the ccriteria of accuracy
        elif self.k_gamma < 90:
            if (
                abs(
                    self.b * cos(self.alpha / 180 * pi) / self.c
                    + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                    - 1
                )
                < 10e-8
            ):
                return "MCLC4"
            elif (
                self.b * cos(self.alpha / 180 * pi) / self.c
                + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                < 1
            ):
                return "MCLC3"
            elif (
                self.b * cos(self.alpha / 180 * pi) / self.c
                + self.b**2 * sin(self.alpha / 180 * pi) ** 2 / self.a**2
                > 1
            ):
                return "MCLC5"


# 14
class TRI(Lattice):
    r"""Triclinic (TRI, aP)"""

    def __init__(
        self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.cell = np.array(
            [
                [a, 0, 0],
                [b * cos(gamma / 180 * pi), b * sin(gamma / 180 * pi), 0],
                [
                    c * cos(beta / 180 * pi),
                    c
                    / sin(gamma / 180 * pi)
                    * (
                        cos(alpha / 180 * pi)
                        - cos(beta / 180 * pi) * cos(gamma / 180 * pi)
                    ),
                    c
                    / sin(gamma / 180 * pi)
                    * sqrt(
                        sin(gamma / 180 * pi) ** 2
                        - cos(alpha / 180 * pi) ** 2
                        - cos(beta / 180 * pi) ** 2
                        + 2
                        * cos(alpha / 180 * pi)
                        * cos(beta / 180 * pi)
                        * cos(gamma / 180 * pi)
                    ),
                ],
            ]
        )
        self.primitive = True
        if self.variant in ["TRI1a", "TRI1b"]:
            self.points = {
                "G": np.array([0, 0, 0]),
                "L": np.array([1 / 2, 1 / 2, 0]),
                "M": np.array([0, 1 / 2, 1 / 2]),
                "N": np.array([1 / 2, 0, 1 / 2]),
                "R": np.array([1 / 2, 1 / 2, 1 / 2]),
                "X": np.array([1 / 2, 0, 0]),
                "Y": np.array([0, 1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }

            self.path = [
                ["X", "G", "Y"],
                ["L", "G", "Z"],
                ["N", "G", "M"],
                "R",
                "G",
            ]
        elif self.variant in ["TRI2a", "TRI2b"]:
            self.points = {
                "G": np.array([0, 0, 0]),
                "L": np.array([1 / 2, 1 / 2, 0]),
                "M": np.array([0, 1 / 2, 1 / 2]),
                "N": np.array([1 / 2, 0, 1 / 2]),
                "R": np.array([1 / 2, 1 / 2, 1 / 2]),
                "X": np.array([1 / 2, 0, 0]),
                "Y": np.array([0, 1 / 2, 0]),
                "Z": np.array([0, 0, 1 / 2]),
            }

            self.path = [
                ["X", "G", "Y"],
                ["L", "G", "Z"],
                ["N", "G", "M"],
                "R",
                "G",
            ]

    @property
    def variant(self):
        r"""
        Four variants of the Lattice.


        :math:`\text{TRI}_{1a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`
        :math:`\text{TRI}_{1b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`
        :math:`\text{TRI}_{2a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`
        :math:`\text{TRI}_{2b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`
        """
        if abs(self.k_gamma - 90) < 10e-5:
            if self.k_alpha > 90 and self.k_beta > 90:
                return "TRI2a"
            elif self.k_alpha < 90 and self.k_beta < 90:
                return "TRI2b"
        elif (min(self.k_gamma, self.k_beta, self.k_alpha)) > 90:
            return "TRI1a"
        elif (max(self.k_gamma, self.k_beta, self.k_alpha)) < 90:
            return "TRI1a"


# Examples
cub = CUB(pi)
fcc = FCC(pi)
bcc = BCC(pi)
tet = TET(pi, 2 * pi)
bct1 = BCT(2 * pi, pi)
bct2 = BCT(pi, 2 * pi)
orc = ORC(pi, 2 * pi, 3 * pi)
orcf1 = ORCF(0.9 * pi, 5 / 4 * pi, 5 / 3 * pi)
orcf2 = ORCF(1.1 * pi, 5 / 4 * pi, 5 / 3 * pi)
orcf3 = ORCF(pi, 5 / 4 * pi, 5 / 3 * pi)
orci = ORCI(pi, 2 * pi, 3 * pi)
orcc = ORCC(pi, 2 * pi, 3 * pi)
hex = HEX(pi, 2 * pi)
rhl1 = RHL(pi, 70)
rhl2 = RHL(pi, 110)
mcl = MCL(pi, 2 * pi, 3 * pi, alpha=80)
mclc1 = MCLC(1 * pi, 1.5 * pi, 2 * pi, 80)
mclc2 = MCLC(1.47721 * pi, 1.5 * pi, 2 * pi, 80)
mclc3 = MCLC(pi, pi / 2, pi, 80)
mclc4 = MCLC(1.06486355 * pi, pi, 1.2 * pi, 80)
mclc5 = MCLC(pi, pi, pi, 60)
tri1a = TRI(2 * pi, 3 * pi, 4 * pi, 60, 70, 80)
tri1b = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
tri2a = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
tri2b = TRI(pi, 2 * pi, 3 * pi, 100, 70, 65)
examples = [
    cub,
    fcc,
    bcc,
    tet,
    bct1,
    bct2,
    orc,
    orcf1,
    orcf2,
    orcf3,
    orci,
    orcc,
    hex,
    rhl1,
    rhl2,
    mcl,
    mclc1,
    mclc2,
    mclc3,
    mclc4,
    mclc5,
    tri1a,
    tri1b,
    tri2a,
    tri2b,
]

if __name__ == "__main__":
    from math import pi

    print(
        f"BCT1 {bct1.variant}",
        f"BCT2 {bct2.variant}",
        f"ORCF1 {orcf1.variant}",
        f"ORCF2 {orcf2.variant}",
        f"ORCF3 {orcf3.variant}",
        f"RHL1 {rhl1.variant}",
        f"RHL2 {rhl2.variant}",
        f"MCLC1 {mclc1.variant}",
        f"MCLC2 {mclc2.variant}",
        f"MCLC3 {mclc3.variant}",
        f"MCLC4 {mclc4.variant}",
        f"MCLC5 {mclc5.variant}",
        f"TRI1a {tri1a.variant}",
        f"TRI1b {tri1b.variant}",
        f"TRI2a {tri2a.variant}",
        f"TRI2b {tri2b.variant}",
        sep="\n",
    )

    for e in examples:
        l = tri1a
        print(l.variant)
        l.prepare_figure()
        l.plot("brillouin_kpath")
        l.show()
        break


# TODO FIX TRI Lattice
