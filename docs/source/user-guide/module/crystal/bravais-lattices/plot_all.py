from argparse import ArgumentParser
from os.path import join
import radtools as rad
import matplotlib.pyplot as plt


def plot(output_path="."):
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
    }

    for i, name in enumerate(names):
        output_subname = (
            name.replace("1", "")
            .replace("2", "")
            .replace("3", "")
            .replace("4", "")
            .replace("5", "")
            .replace("1a", "")
            .replace("2a", "")
            .replace("1b", "")
            .replace("2b", "")
        )
        l = rad.lattice_example(name)
        for j, wtp in enumerate(wtps[name]):
            py_file = open(
                join(
                    output_path, output_subname, f"{name.lower()}_{names[name][j]}.py"
                ),
                "w",
            )
            py_file.write(
                f'import radtools as rad\n\nl = rad.lattice_example(f"{name}")\n'
            )
            for data in wtp:
                if len(wtp) == 1:
                    l.plot(data)
                    py_file.write(f'l.plot("{data}")\n')
                else:
                    if data == "conventional":
                        l.plot(data, colour="black", label=data)
                        py_file.write(
                            f"l.plot(\n"
                            + f'    "{data}",\n'
                            + f'    label="{data}",\n'
                            + '    colour="black"\n'
                            + ")\n"
                        )
                    else:
                        l.plot(data, label=data)
                        py_file.write(
                            f"l.plot(\n"
                            + f'    "{data}",\n'
                            + f'    label="{data}",\n'
                            + ")\n"
                        )
                    py_file.write("l.legend()\n")
                    l.legend()
            l.savefig(
                join(
                    output_path, output_subname, f"{name.lower()}_{names[name][j]}.png"
                ),
                elev=elev[name][j],
                azim=azim[name][j],
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()
            py_file.write(
                "# Save an image:\n"
                + "l.savefig(\n"
                + f'    "{name.lower()}_{names[name][j]}.png",\n'
                + f"    elev={elev[name][j]},\n"
                + f"    azim={azim[name][j]},\n"
                + "    dpi=300,\n"
                + '   bbox_inches="tight",\n'
                + ")\n"
            )
            py_file.write(
                "# Interactive plot:\n"
                + f"l.show(elev={elev[name][j]}, azim={azim[name][j]})\n"
            )
            l.ax = None
            l.fig = None
            py_file.close()
        print(f"{name} done")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-op",
        "--output-path",
        metavar="path",
        type=str,
        default=".",
        help="Folder for output.",
    )

    args = parser.parse_args()

    plot(args.output_path)
