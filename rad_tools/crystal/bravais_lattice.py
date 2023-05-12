r"""14 Bravais lattice"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from rad_tools.crystal.lattice import Lattice
from rad_tools import print_2D_array

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


class BravaisLattice(Lattice):
    PLOT_NAMES = {}

    def __init__(self, a, b, c) -> None:
        super().__init__(a, b, c)
        self.primitive = False
        self.points = {}
        self.path = []

    @property
    def reciprocal_cell(self):
        return np.array(
            [
                2 * pi / self.unit_cell_volume * np.cross(self.a2, self.a3),
                2 * pi / self.unit_cell_volume * np.cross(self.a3, self.a1),
                2 * pi / self.unit_cell_volume * np.cross(self.a1, self.a2),
            ]
        )

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
        r"""There is no variants of the Lattice"""
        return "Only one variant."

    def pepare_figure(self, background=True, focal_length=0.2) -> None:
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

    def plot_real_space(self, ax, vectors=True, colour=BLUE):
        if vectors:
            ax.text(
                self.cell[0][0] * 1.1,
                self.cell[0][1],
                self.cell[0][2],
                "$a_1$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.cell[1][0],
                self.cell[1][1] * 1.1,
                self.cell[1][2],
                "$a_2$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.cell[2][0],
                self.cell[2][1],
                self.cell[2][2] * 1.1,
                "$a_3$",
                fontsize=20,
                color=colour,
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

    def plot_reciprocal_space(self, ax, vectors=True, colour=RED):
        planes, edges, corners = deduct_zone([self.b1, self.b2, self.b3])

        if vectors:
            ax.text(
                self.b1[0] * 1.1,
                self.b1[1],
                self.b1[2],
                "$b_1$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.b2[0],
                self.b2[1] * 1.1,
                self.b2[2],
                "$b_2$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.b3[0],
                self.b3[1],
                self.b3[2] * 1.1,
                "$b_3$",
                fontsize=20,
                color=colour,
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

    def plot_kpath(self, ax, colour="black", **kwargs):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        for point in self.points:
            ax.scatter(
                *tuple(self.points[point] @ self.reciprocal_cell),
                s=64,
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

    def plot_wigner_seitz(self, ax, vectors=True, colour="black", **kwargs):
        reverse = False
        if not self.primitive:
            self.make_primitive()
            reverse = True
        planes, edges, corners = deduct_zone([self.a1, self.a2, self.a3])

        if vectors:
            ax.text(
                self.a1[0] * 1.1,
                self.a1[1],
                self.a1[2],
                "$a_1$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.a2[0],
                self.a2[1] * 1.1,
                self.a2[2],
                "$a_2$",
                fontsize=20,
                color=colour,
            )
            ax.text(
                self.a3[0],
                self.a3[1],
                self.a3[2] * 1.1,
                "$a_3$",
                fontsize=20,
                color=colour,
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


class CUB(BravaisLattice):
    r"""Cubic (CUB, cP)"""

    PLOT_NAMES = {"G": "$\\Gamma$", "M": "$M$", "R": "$R$", "X": "$X$"}

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


class FCC(BravaisLattice):
    r"""Face-centred cubic (FCC, cF)"""

    PLOT_NAMES = {
        "G": "$\\Gamma$",
        "K": "$K$",
        "L": "L",
        "U": "U",
        "W": "W",
        "X": "X",
    }

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


class BCC(BravaisLattice):
    r"""Body-centered cubic (BCC, cl)"""

    PLOT_NAMES = {"G": "$\\Gamma$", "H": "$H$", "P": "$P$", "N": "$N$"}

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


class TET(BravaisLattice):
    r"""Tetragonal (TET, tP)"""

    PLOT_NAMES = {
        "G": "$\\Gamma$",
        "A": "$A$",
        "M": "$M$",
        "R": "$R$",
        "X": "$X$",
        "Z": "$Z$",
    }

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


class BCT(BravaisLattice):
    r"""Body-centred tetragonal (BCT, tI)"""

    PLOT_NAMES = {
        "G": "$\\Gamma$",
        "M": "$M$",
        "N": "$N$",
        "P": "$P$",
        "X": "$X$",
        "Z": "$Z$",
        "Z1": "$Z_1$",
        "S": "$\\Sigma$",
        "S1": "$\\Sigma_1$",
        "Y": "$Y$",
        "Y1": "$Y_1$",
    }

    def __init__(self, a: float, c: float) -> None:
        if a == c:
            raise ValueError("Are you trying to create BCC Lattice (a == c)?")
        self.a = a
        self.c = c
        self.cell = np.diag([a, a, c])
        self.primitive = False
        if self.variant == "BCT1":
            eta = (1 + c**2 / a**2) / 4
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
            eta = (1 + a**2 / c**2) / 4
            zeta = a**2 / (2 * c**2)
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


if __name__ == "__main__":
    from math import pi

    focal_length = 0.2

    l = BCT(pi, 2 * pi)
    l.pepare_figure()
    l.plot(kind="primitive")
    l.plot(kind="kpath")
    l.plot(kind="brillouin_kpath")
    l.show()
    l.plot(kind="primitive")
