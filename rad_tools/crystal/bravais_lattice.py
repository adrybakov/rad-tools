r"""14 Bravais lattice"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from rad_tools.crystal.lattice import Lattice
from rad_tools import print_2D_array

__all__ = ["BravaisLattice", "CUB"]


class BravaisLattice(Lattice):
    def __init__(self, a, b, c) -> None:
        super().__init__(a, b, c)
        self.lattice_points = []

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

    def plot(self, kind="conventional"):
        rcParams["axes.linewidth"] = 0
        rcParams["xtick.color"] = "#B3B3B3"
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        ax.set_proj_type("persp", focal_length=0.2)
        ax.set_aspect("equal")
        ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
        ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
        ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 1)
        ax.set_xlabel("x", fontsize=15, alpha=0.5)
        ax.set_ylabel("y", fontsize=15, alpha=0.5)
        ax.set_zlabel("z", fontsize=15, alpha=0.5)
        ax.tick_params(axis="both", zorder=0, color="#B3B3B3")
        try:
            getattr(self, f"plot_{kind}")(ax)
        except AttributeError:
            d = dir(self)
            for i in d:
                print(i, i == f"plot_{kind}")
            raise ValueError(f"Plot kind '{kind}' does not exist!")
        plt.show()

    def plot_conventional(self, ax):
        for i in self.lattice_points:
            ax.scatter(*tuple(np.matmul(np.array(i), self.cell)), color="black", s=64)
        ax.text(
            self.cell[0][0] * 1.1,
            self.cell[0][1],
            self.cell[0][2],
            "$a_1$",
            fontsize=20,
        )
        ax.text(
            self.cell[1][0],
            self.cell[1][1] * 1.1,
            self.cell[1][2],
            "$a_2$",
            fontsize=20,
        )
        ax.text(
            self.cell[2][0],
            self.cell[2][1],
            self.cell[2][2] * 1.1,
            "$a_3$",
            fontsize=20,
        )
        for i in self.cell:
            ax.quiver(0, 0, 0, *tuple(i), arrow_length_ratio=0.2, color="black")

        def plot_line(line, shift):
            ax.plot(
                [shift[0], shift[0] + line[0]],
                [shift[1], shift[1] + line[1]],
                [shift[2], shift[2] + line[2]],
                color="black",
            )

        for i in range(0, 3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            plot_line(self.cell[i], np.zeros(3))
            plot_line(self.cell[i], self.cell[j])
            plot_line(self.cell[i], self.cell[k])
            plot_line(self.cell[i], self.cell[j] + self.cell[k])

    def plot_primitive(self, ax):
        self.make_primitive()
        self.plot_conventional(ax)
        self.make_conventional()

    def plot_brillouin(self, ax):
        self.make_primitive()
        planes = []
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    vector = i * self.b1 + j * self.b2 + k * self.b3
                    d_gamma = np.linalg.norm(vector / 2)
                    if i**2 + j**2 + k**2 != 0:
                        planes.append((np.array((i, j, k)), d_gamma))
        planes.sort(key=lambda x: x[1])
        for i in range(len(planes)):
            if planes[i][1] > planes[0][1]:
                break
        planes = planes[:i]
        for i in range(len(planes)):
            print(planes[i][0])
            planes[i] = planes[i][0]
        points = []
        for i in range(len(planes)):
            n_i = (
                planes[i][0] * self.b1 + planes[i][1] * self.b2 + planes[i][2] * self.b3
            ) / 2
            for j in range(i + 1, len(planes)):
                n_j = (
                    planes[j][0] * self.b1
                    + planes[j][1] * self.b2
                    + planes[j][2] * self.b3
                ) / 2
                for k in range(j + 1, len(planes)):
                    n_k = (
                        planes[k][0] * self.b1
                        + planes[k][1] * self.b2
                        + planes[k][2] * self.b3
                    ) / 2
                    A = np.array([n_i, n_j, n_k])
                    b = np.array(
                        [np.linalg.norm(n_i), np.linalg.norm(n_j), np.linalg.norm(n_k)]
                    )
                    try:
                        x = np.linalg.solve(A, b)
                        points.append((x, np.linalg.norm(x)))
                    except:
                        pass

        points.sort(key=lambda x: x[1])
        for i in range(len(points)):
            if points[i][1] > points[0][1]:
                break
        points = points[:i]
        for i in range(len(points)):
            print(points[i][0])
            points[i] = points[i][0]

        ax.text(
            self.b1[0] * 1.1,
            self.b1[1],
            self.b1[2],
            "$b_1$",
            fontsize=20,
        )
        ax.text(
            self.b2[0],
            self.b2[1] * 1.1,
            self.b2[2],
            "$b_2$",
            fontsize=20,
        )
        ax.text(
            self.b3[0],
            self.b3[1],
            self.b3[2] * 1.1,
            "$b_3$",
            fontsize=20,
        )
        for i in [self.b1, self.b2, self.b3]:
            ax.quiver(0, 0, 0, *tuple(i), arrow_length_ratio=0.2, color="black")
        self.make_conventional()

    def plot_kpath(self, ax):
        pass

    def plot_wigner_seitz(self, ax):
        pass


class CUB(BravaisLattice):
    r"""Cubic (CUB, cP)"""

    def __init__(self, a: float) -> None:
        self.cell = np.diag([a, a, a])
        self.primitive = True
        self.lattice_points = [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1],
        ]


if __name__ == "__main__":
    from math import pi

    l = CUB(pi)
    # l.plot()
    # l.plot(kind="primitive")
    l.plot(kind="brillouin")
    # l.plot(kind="kpath")
    # l.plot(kind="wigner_seitz")
