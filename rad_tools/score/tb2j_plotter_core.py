from argparse import ArgumentParser
from math import sqrt
from os import makedirs
from os.path import abspath, join

import numpy as np
from matplotlib import pyplot as plt

from rad_tools.io.internal import read_template
from rad_tools.io.tb2j import read_exchange_model
from rad_tools.routines import OK, RESET, atom_mark_to_latex, rot_angle


def manager(
    input_filename,
    output_path=".",
    output_name="exchange",
    what_to_plot="iso",
    draw_cells=False,
    min_distance=None,
    max_distance=None,
    distance=None,
    template_file=None,
    R_vector=None,
    scale_atoms=1,
    scale_data=1,
    title=None,
    force_symmetry=False,
    verbose=False,
):
    r"""
    Main function call of the tb2j-plotter.py script.

    Full documentation on the behaviour is available
    :ref:`here <tb2j-plotter>`. Parameters of the function directly
    corresponds to the arguments of the script.

    If you want to have the behaviour of the tb2j-plotter.py script
    but in a format of a function call use this function.
    """

    if force_symmetry and template_file is None:
        raise ValueError("--force-symmetry option requires a template file.")

    if distance is not None:
        min_distance = distance
        max_distance = distance

    try:
        makedirs(output_path)
    except FileExistsError:
        pass

    if R_vector is not None:
        R_vector = np.array(R_vector[: len(R_vector) // 3 * 3], dtype=int).reshape(
            (len(R_vector) // 3, 3)
        )
        R_vector = list(map(tuple, R_vector.tolist()))

    model = read_exchange_model(input_filename, quiet=not verbose)
    if template_file is not None:
        template = read_template(template_file)
    else:
        template = None

    if force_symmetry:
        model.force_symmetry(template=template)

    model.filter(
        min_distance=min_distance,
        max_distance=max_distance,
        R_vector=R_vector,
        template=template,
    )

    dummy = True
    ha = "center"
    x_min, y_min, z_min, x_max, y_max, z_max = model.space_dimensions
    X = max(abs(x_min), abs(x_max))
    Y = max(abs(y_min), abs(y_max))
    if X == 0 and Y == 0:
        X = Y = 1
    linewidth = 1
    fontsize = 11
    if what_to_plot == "all":
        wtps = ["iso", "distance", "dmi"]
    else:
        wtps = [what_to_plot]
    for wtp in wtps:
        if X < Y:
            fig = plt.figure(figsize=(6.4 * X / Y, 4.8))
        else:
            fig = plt.figure(figsize=(6.4, 4.8 * Y / X))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.set_aspect("equal")
        ax.set_xlabel("x, Angstroms")
        ax.set_ylabel("y, Angstroms")

        for atom1, atom2, R in model:
            bond = model[(atom1, atom2, R)]
            dis = model.get_distance(atom1, atom2, R)
            x1, y1, z1 = model.get_atom_coordinates(atom1)
            x2, y2, z2 = model.get_atom_coordinates(atom2, R)
            xm = (x1 + x2) / 2
            ym = (y1 + y2) / 2
            zm = (z1 + z2) / 2

            ax.scatter(x1, y1, s=100 * fig.dpi / 72.0, c="white")
            ax.scatter(x2, y2, s=100 * fig.dpi / 72.0, c="white")

            ax.text(
                x2,
                y2,
                atom_mark_to_latex(atom2),
                va="center",
                ha="center",
                fontsize=1.5 * fontsize * scale_atoms,
            )
            ax.text(
                x1,
                y1,
                atom_mark_to_latex(atom1),
                va="center",
                ha="center",
                fontsize=1.5 * fontsize * scale_atoms,
                color="#A04F4D",
            )
            if wtp == "iso":
                ax.text(
                    xm,
                    ym,
                    str(round(bond.iso, 4)),
                    va="bottom",
                    ha=ha,
                    rotation_mode="anchor",
                    rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                    fontsize=fontsize * scale_data,
                )
            elif wtp == "distance":
                ax.text(
                    xm,
                    ym,
                    str(round(dis, 4)),
                    va="bottom",
                    ha=ha,
                    rotation_mode="anchor",
                    rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                    fontsize=fontsize * scale_data,
                )
            elif wtp == "dmi":
                ax.text(
                    xm,
                    ym,
                    str(round(sqrt(np.sum(bond.dmi**2)), 4)),
                    va="bottom",
                    ha=ha,
                    rotation_mode="anchor",
                    rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                    fontsize=fontsize * scale_data,
                )

        xlims = (ax.get_xlim()[0] - 0.5, ax.get_xlim()[1] + 0.5)
        ylims = (ax.get_ylim()[0] - 0.5, ax.get_ylim()[1] + 0.5)

        if draw_cells:
            cells = model.cell_list
            a_x, a_y, a_z = tuple(model.cell[0])
            b_x, b_y, b_z = tuple(model.cell[1])
            Rx_min = 0
            Rx_max = 0
            Ry_min = 0
            Ry_max = 0
            Rz_min = 0
            Rz_max = 0
            for Rx, Ry, Rz in cells:
                Rx_min = min(Rx, Rx_min)
                Rx_max = max(Rx, Rx_max)
                Ry_min = min(Ry, Ry_min)
                Ry_max = max(Ry, Ry_max)
                Rz_min = min(Rz, Rz_min)
                Rz_max = max(Rz, Rz_max)

            for i in range(Rx_min, Rx_max + 2):
                for j in range(Ry_min, Ry_max + 2):
                    for k in range(Rz_min, Rz_max + 2):
                        ax.plot(
                            [Rx_min * a_x + j * b_x, (Rx_max + 1) * a_x + j * b_x],
                            [Rx_min * a_y + j * b_y, (Rx_max + 1) * a_y + j * b_y],
                            linewidth=linewidth,
                            color="#DFDFDF",
                        )
                        ax.plot(
                            [i * a_x + Ry_min * b_x, i * a_x + (Ry_max + 1) * b_x],
                            [i * a_y + Ry_min * b_y, i * a_y + (Ry_max + 1) * b_y],
                            linewidth=linewidth,
                            color="#DFDFDF",
                        )
            ax.plot(
                np.array([0, a_x, a_x + b_x, b_x, 0]),
                np.array([0, a_y, a_y + b_y, b_y, 0]),
                linewidth=linewidth,
                color="#CA7371",
            )

        ax.set_xlim(*xlims)
        ax.set_ylim(*ylims)

        if title is not None:
            ax.set_title(title, fontsize=1.5 * fontsize)

        png_path = join(output_path, f"{output_name}.{wtp}.png")
        plt.savefig(png_path, dpi=400, bbox_inches="tight")
        print(f"{OK}2D plot with {wtp} is in {abspath(png_path)}{RESET}")


def create_parser():
    parser = ArgumentParser(description="Script for visualization of TB2J results")

    parser.add_argument(
        "-if",
        "--input-filename",
        metavar="filename",
        type=str,
        required=True,
        help="Relative or absolute path to the 'exchange.out' file,"
        + "including the name and extension of the file itself.",
    )
    parser.add_argument(
        "-op",
        "--output-path",
        metavar="path",
        type=str,
        default=".",
        help="Relative or absolute path to the folder for saving outputs.",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        metavar="filename",
        type=str,
        default="exchange",
        help="Seedname for the output files.",
    )
    parser.add_argument(
        "-wtp",
        "--what-to-plot",
        metavar="value",
        type=str,
        choices=["all", "iso", "distance", "dmi"],
        default="all",
        help="Type of data for display.",
    )
    parser.add_argument(
        "-dc",
        "--draw-cells",
        action="store_true",
        default=False,
        help="Whenever to draw the cells.",
    )
    parser.add_argument(
        "-R",
        "--R-vector",
        metavar="R1a R1b R1c",
        type=int,
        nargs="*",
        default=None,
        help="R vectors for filtering the model.",
    )
    parser.add_argument(
        "-maxd",
        "--max-distance",
        metavar="distance",
        type=float,
        default=None,
        help="(<=) Maximum distance.",
    )
    parser.add_argument(
        "-mind",
        "--min-distance",
        metavar="distance",
        type=float,
        default=None,
        help="(>=) Minimum distance.",
    )
    parser.add_argument(
        "-d",
        "--distance",
        metavar="distance",
        type=float,
        default=None,
        help="(=) Exact distance.",
    )
    parser.add_argument(
        "-tf",
        "--template-file",
        metavar="filename",
        type=str,
        default=None,
        help="Relative or absolute path to the template file, "
        + "including the name and extension of the file.",
    )
    parser.add_argument(
        "-sa",
        "--scale-atoms",
        metavar="factor",
        default=1,
        type=float,
        help="Scale for the size of atom marks.",
    )
    parser.add_argument(
        "-sd",
        "--scale-data",
        metavar="factor",
        default=1,
        type=float,
        help="Scale for the size of data text.",
    )
    parser.add_argument(
        "-t",
        "--title",
        metavar="title",
        default=None,
        type=str,
        help="Title for the plots.",
    )
    parser.add_argument(
        "-fs",
        "--force-symmetry",
        default=False,
        action="store_true",
        help="Force the Exchange model to have the symmetry of the template.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )

    return parser
