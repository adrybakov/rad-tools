r"""14 Bravais lattice"""

import numpy as np
import matplotlib.pyplot as plt
from rad_tools.crystal.lattice import Lattice

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

    def plot(self, type="conventional"):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        ax.set_proj_type("persp", focal_length=0.2)
        ax.set_aspect("equal")
        ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
        ax.set_xlabel("x", fontsize=15, alpha=0.5)
        ax.set_ylabel("y", fontsize=15, alpha=0.5)
        ax.set_zlabel("z", fontsize=15, alpha=0.5)
        types = {
            "conventional": self.plot_conventional,
            "primitive": self.plot_primitive,
            "brillouin": self.plot_brillouin,
            "kpath": self.plot_kpath,
            "wigner_seitz": self.plot_wigner_seitz,
        }
        types[type](ax)
        plt.show()

    def plot_conventional(self, ax):
        for i in self.lattice_points:
            ax.scatter(*tuple(np.matmul(np.array(i), self.cell)), color="black", s=64)
        ax.text(
            self.cell[0][0] * 1.1,
            self.cell[0][1],
            self.cell[0][2],
            "$a_1$",
            fontsize=15,
        )
        ax.text(
            self.cell[1][0],
            self.cell[1][1] * 1.1,
            self.cell[1][2],
            "$a_2$",
            fontsize=15,
        )
        ax.text(
            self.cell[2][0],
            self.cell[2][1],
            self.cell[2][2] * 1.1,
            "$a_3$",
            fontsize=15,
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
        pass

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
    l = CUB(23)
    l.plot()
    l.plot(type="primitive")
