# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
from argparse import ArgumentParser

import matplotlib.pyplot as plt

import radtools as rad

ROOT_DIR = "."
OUTPUT_PATH = os.path.join(
    ROOT_DIR, "docs", "source", "user-guide", "library", "bravais-lattices"
)


def plot():
    elev = {
        "CUB": [28, 28, 28],
        "FCC": [23, 28, 46],
        "BCC": [11, 30, 46],
        "TET": [30, 30, 30],
        "BCT1": [30, 37, 26],
        "BCT2": [36, 40, 41],
        "ORC": [35, 36, 20],
        "ORCF1": [21, 24, 44],
        "ORCF2": [15, 25, 38],
        "ORCF3": [25, 27, 23],
        "ORCI": [35, 32, 30],
        "ORCC": [22, 39, 33],
        "HEX": [19, 35, 32],
        "RHL1": [-41, 35, 19],
        "RHL2": [14, 35, 30],
        "MCL": [12, 25, 11],
        "MCLC1": [21, 29, 45],
        "MCLC2": [11, 30, 13],
        "MCLC3": [23, 39, 41],
        "MCLC4": [19, 25, 5],
        "MCLC5": [-5, 20, 13],
        "TRI1a": [45, 31, 9],
        "TRI2a": [10, 39, 30],
        "TRI1b": [30, 12, 21],
        "TRI2b": [40, 17, 19],
    }
    azim = {
        "CUB": [23, 23, 23],
        "FCC": [28, 23, 19],
        "BCC": [25, -180, 19],
        "TET": [23, 30, 30],
        "BCT1": [28, 72, 59],
        "BCT2": [28, 85, 59],
        "ORC": [34, 35, 30],
        "ORCF1": [49, 38, 28],
        "ORCF2": [36, 28, 14],
        "ORCF3": [62, 36, 38],
        "ORCI": [23, -12, 12],
        "ORCC": [57, 15, -40],
        "HEX": [20, 23, 10],
        "RHL1": [-13, 52, -19],
        "RHL2": [-85, 52, -29],
        "MCL": [25, 40, 37],
        "MCLC1": [53, 83, 67],
        "MCLC2": [64, 62, 29],
        "MCLC3": [34, 68, 61],
        "MCLC4": [45, 70, 66],
        "MCLC5": [59, 56, 66],
        "TRI1a": [-55, -20, 18],
        "TRI2a": [32, 44, 62],
        "TRI1b": [42, 11, 39],
        "TRI2b": [1, 54, 13],
    }

    wtps = {
        "CUB": [["brillouin-kpath"], ["primitive"], ["wigner-seitz"]],
        "FCC": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "BCC": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "TET": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "BCT1": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "BCT2": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "ORC": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "ORCF1": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "ORCF2": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "ORCF3": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "ORCI": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "ORCC": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "HEX": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "RHL1": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "RHL2": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "MCL": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "MCLC1": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "MCLC2": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "MCLC3": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "MCLC4": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "MCLC5": [
            ["brillouin-kpath"],
            ["primitive", "conventional"],
            ["wigner-seitz"],
        ],
        "TRI1a": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "TRI2a": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "TRI1b": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
        "TRI2b": [
            ["brillouin-kpath"],
            ["primitive"],
            ["wigner-seitz"],
        ],
    }
    names = {
        "CUB": ["brillouin", "real", "wigner-seitz"],
        "FCC": ["brillouin", "real", "wigner-seitz"],
        "BCC": ["brillouin", "real", "wigner-seitz"],
        "TET": ["brillouin", "real", "wigner-seitz"],
        "BCT1": ["brillouin", "real", "wigner-seitz"],
        "BCT2": ["brillouin", "real", "wigner-seitz"],
        "ORC": ["brillouin", "real", "wigner-seitz"],
        "ORCF1": ["brillouin", "real", "wigner-seitz"],
        "ORCF2": ["brillouin", "real", "wigner-seitz"],
        "ORCF3": ["brillouin", "real", "wigner-seitz"],
        "ORCI": ["brillouin", "real", "wigner-seitz"],
        "ORCC": ["brillouin", "real", "wigner-seitz"],
        "HEX": ["brillouin", "real", "wigner-seitz"],
        "RHL1": ["brillouin", "real", "wigner-seitz"],
        "RHL2": ["brillouin", "real", "wigner-seitz"],
        "MCL": ["brillouin", "real", "wigner-seitz"],
        "MCLC1": ["brillouin", "real", "wigner-seitz"],
        "MCLC2": ["brillouin", "real", "wigner-seitz"],
        "MCLC3": ["brillouin", "real", "wigner-seitz"],
        "MCLC4": ["brillouin", "real", "wigner-seitz"],
        "MCLC5": ["brillouin", "real", "wigner-seitz"],
        "TRI1a": ["brillouin", "real", "wigner-seitz"],
        "TRI2a": ["brillouin", "real", "wigner-seitz"],
        "TRI1b": ["brillouin", "real", "wigner-seitz"],
        "TRI2b": ["brillouin", "real", "wigner-seitz"],
    }

    for i, name in enumerate(names):
        output_subname = (name.translate(str.maketrans("", "", "12345ab"))).lower()
        l = rad.lattice_example(name)
        for j, wtp in enumerate(wtps[name]):
            py_file = open(
                os.path.join(
                    OUTPUT_PATH, output_subname, f"{name.lower()}_{names[name][j]}.py"
                ),
                "w",
                encoding="utf-8",
            )
            py_file.write(
                f'import radtools as rad\n\nl = rad.lattice_example("{name}")\nbackend = rad.PlotlyBackend()\n'
            )
            backend = rad.PlotlyBackend()
            for data in wtp:
                if len(wtp) == 1:
                    backend.plot(l, kind=data)
                    py_file.write(f'backend.plot(l, kind="{data}")\n')
                else:
                    if data == "conventional":
                        backend.plot(l, kind=data, color="black", label=data)
                        py_file.write(
                            f'backend.plot(l, kind="{data}", label="{data}", color="black")\n'
                        )
                    else:
                        backend.plot(l, kind=data, label=data)
                        py_file.write(
                            f'backend.plot(l, kind="{data}", label="{data}")\n'
                        )
            backend.save(
                os.path.join(
                    OUTPUT_PATH, output_subname, f"{name.lower()}_{names[name][j]}.html"
                ),
                kwargs_write_html=dict(full_html=False, include_plotlyjs=False),
                kwargs_update_layout=dict(
                    width=600,
                    height=600,
                ),
            )
            plt.close()
            py_file.write(
                "# Save an image:\n"
                + f'backend.save("{name.lower()}_{names[name][j]}.png")\n'
            )
            py_file.write("# Interactive plot:\n" + f"backend.show()\n")
            del backend
            py_file.close()
        print(f"{name} done")


if __name__ == "__main__":
    plot()
