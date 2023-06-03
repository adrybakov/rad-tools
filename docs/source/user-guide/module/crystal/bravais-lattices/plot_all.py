from argparse import ArgumentParser
from os.path import join
import radtools as rad


def plot(output_path="."):
    lattices = rad.lattice_example()
    elev = {
        "CUB": [28, 28, 28],
        "FCC": [23, 28, 46],
        "BCC": [11, 30, 46],
        "TET": [30, 30, 30],
    }
    azim = {
        "CUB": [23, 23, 23],
        "FCC": [28, 23, 19],
        "BCC": [25, -180, 19],
        "TET": [23, 30, 30],
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
    }
    names = {
        "CUB": ["brillouin", "real", "wigner-seitz"],
        "FCC": ["brillouin", "real", "wigner-seitz"],
        "BCC": ["brillouin", "real", "wigner-seitz"],
        "TET": ["brillouin", "real", "wigner-seitz"],
    }

    for i, name in enumerate(lattices[:4]):
        l = rad.lattice_example(name)
        for j, wtp in enumerate(wtps[name]):
            py_file = open(
                join(output_path, name.lower(), f"plot_{names[name][j]}.py"), "w"
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
                join(output_path, name.lower(), f"{name.lower()}_{names[name][j]}.png"),
                elev=elev[name][j],
                azim=azim[name][j],
                dpi=300,
                bbox_inches="tight",
            )
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
