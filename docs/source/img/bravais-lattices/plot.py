from argparse import ArgumentParser
from os.path import join
import radtools as rad


def plot(output_path="."):
    lattices = rad.lattice_example()
    elev = {"CUB": [28, 28, 28], "FCC": [23, 28, 46], "BCC": [11, 30, 46]}
    azim = {"CUB": [23, 23, 23], "FCC": [28, 23, 19], "BCC": [25, -180, 19]}
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
    }
    names = {
        "CUB": ["brillouin", "primitive", "wigner-seitz"],
        "FCC": ["brillouin", "both", "wigner-seitz"],
        "BCC": ["brillouin", "both", "wigner-seitz"],
    }

    for i, name in enumerate(lattices[:3]):
        l = rad.lattice_example(name)
        for j, wtp in enumerate(wtps[name]):
            for data in wtp:
                if len(wtp) == 1:
                    l.plot(data)
                else:
                    if data == "conventional":
                        l.plot(data, colour="black", label=data)
                    else:
                        l.plot(data, label=data)
            l.savefig(
                join(output_path, f"{name}_{names[name][j]}"),
                elev=elev[name][j],
                azim=azim[name][j],
                dpi=300,
                bbox_inches="tight",
            )
            l.ax = None
            l.fig = None
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
