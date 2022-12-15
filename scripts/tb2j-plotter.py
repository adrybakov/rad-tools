#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join, abspath
from os import makedirs
from math import sqrt

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from rad_tools.io.tb2j import read_exchange_model
from rad_tools.io.internal import read_template
from rad_tools.routines import atom_mark_to_latex, rot_angle, OK, RESET


def plot_2d(filename, out_dir='.',
            out_name='exchange',
            wtp='iso',
            draw_cells=False,
            min_distance=None,
            max_distance=None,
            template=None,
            R_vector=None,
            double_bonds=False,
            scale_atoms=1,
            scale_data=1,
            title=None):

    mode_name = "2d"
    messages = {'iso': "isotropic exchange",
                "distance": "distances",
                "dmi": "dmi"}

    model = read_exchange_model(filename)
    template = read_template(template)
    model.filter(min_distance=min_distance,
                 max_distance=max_distance,
                 R_vector=R_vector,
                 template=template)

    dummy = True
    ha = 'right'
    if not double_bonds:
        model = model.remove_double_bonds()
        dummy = False
        ha = 'center'
    x_min, y_min, z_min, x_max, y_max, z_max = model.space_dimensions
    X = max(abs(x_min), abs(x_max))
    Y = max(abs(y_min), abs(y_max))
    if X == 0 and Y == 0:
        X = Y = 1
    linewidth = 1
    fontsize = 11

    if X < Y:
        fig = plt.figure(figsize=(6.4*X/Y, 4.8))
    else:
        fig = plt.figure(figsize=(6.4, 4.8*Y/X))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_aspect("equal")
    ax.set_xlabel('x, Angstroms')
    ax.set_ylabel('y, Angstroms')

    for atom1 in model.bonds:
        for atom2 in model.bonds[atom1]:
            for R in model.bonds[atom1][atom2]:
                bond = model.bonds[atom1][atom2][R]
                x1, y1, z1 = model.get_atom_coordinates(atom1)
                x2, y2, z2 = model.get_atom_coordinates(atom2, R)
                xm = (x1 + x2) / 2
                ym = (y1 + y2) / 2
                zm = (z1 + z2) / 2

                ax.scatter(x1, y1, s=100 * fig.dpi/72., c='white')
                ax.scatter(x2, y2, s=100 * fig.dpi/72., c='white')

                ax.text(x1, y1, atom_mark_to_latex(atom1),
                        va='center', ha='center',
                        fontsize=1.5 * fontsize * scale_atoms)
                ax.text(x2, y2, atom_mark_to_latex(atom2),
                        va='center', ha='center',
                        fontsize=1.5 * fontsize * scale_atoms)
                if wtp == 'iso':
                    ax.text(xm, ym, str(round(bond.iso, 4)),
                            va='bottom', ha=ha,
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=fontsize * scale_data)
                elif wtp == 'distance':
                    ax.text(xm, ym, str(round(bond.dis, 4)),
                            va='bottom', ha=ha,
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=fontsize * scale_data)
                elif wtp == 'dmi':
                    ax.text(xm, ym, str(round(sqrt(np.sum(bond.dmi**2)), 4)),
                            va='bottom', ha=ha,
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=fontsize * scale_data)

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
                    ax.plot([Rx_min * a_x + j * b_x,
                             (Rx_max + 1) * a_x + j * b_x],
                            [Rx_min * a_y + j * b_y,
                             (Rx_max + 1) * a_y + j * b_y],
                            linewidth=linewidth, color="#DFDFDF")
                    ax.plot([i * a_x + Ry_min * b_x,
                             i * a_x + (Ry_max + 1) * b_x],
                            [i * a_y + Ry_min * b_y,
                            i * a_y + (Ry_max + 1) * b_y],
                            linewidth=linewidth, color="#DFDFDF")
        ax.plot(np.array([0, a_x, a_x + b_x, b_x, 0]),
                np.array([0, a_y, a_y + b_y, b_y, 0]),
                linewidth=linewidth, color="#CA7371")

    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)

    if title is not None:
        ax.set_title(title, fontsize=1.5 * fontsize)

    png_path = join(out_dir, f'{out_name}.{wtp}.png')
    plt.savefig(png_path, dpi=400, bbox_inches="tight")
    print(f'{OK}{mode_name} plot with {wtp} is in {abspath(png_path)}{RESET}')


if __name__ == '__main__':
    parser = ArgumentParser(
        description="Script for visualisation of TB2J results",
        epilog="""
               For the full description of arguments see the docs:
               https://rad-tools.adrybakov.com/en/stable/user-guide/tb2j-plotter.html
               """)

    plot_data_type = ['iso', 'distance', 'dmi']

    parser.add_argument("-f", "--filename",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the *exchange.out* file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder for saving outputs.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='exchange',
                        help="""
                        Seedname for the output files.
                        """
                        )
    parser.add_argument("-wtp", "--what-to-plot",
                        type=str,
                        choices=['all'] + plot_data_type,
                        default='all',
                        help="""
                        Type of data for display.
                        """
                        )
    parser.add_argument("-dc", "--draw-cells",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to draw the cells.
                        """
                        )
    parser.add_argument("-R", "--R-vector",
                        type=int,
                        nargs="*",
                        default=None,
                        help="""
                        R vectors for filtering the model.
                        """
                        )
    parser.add_argument("-maxd", "--max-distance",
                        type=float,
                        default=None,
                        help="""
                        (<=) Maximum distance.
                        """
                        )
    parser.add_argument("-mind", "--min-distance",
                        type=float,
                        default=None,
                        help="""
                        (>=) Minimum distance.
                        """
                        )
    parser.add_argument("-d", "--distance",
                        type=float,
                        default=None,
                        help="""
                        (=) Exact distance.
                        """
                        )
    parser.add_argument("-tf", "--template",
                        type=str,
                        default=None,
                        help="""
                        Relative or absolute path to the template file,
                        including the name and extention of the file.
                        """)
    parser.add_argument("-db", "--double-bonds",
                        default=False,
                        action="store_true",
                        help="""
                        Whenever to keep both bonds.
                        """)
    parser.add_argument("-sa", "--scale-atoms",
                        default=1,
                        type=float,
                        help="""
                        Scale for the size of atom marks.
                        """)
    parser.add_argument("-sd", "--scale-data",
                        default=1,
                        type=float,
                        help="""
                        Scale for the size of data text.
                        """)
    parser.add_argument("-t", "--title",
                        default=None,
                        type=str,
                        help="""
                        Title for the plots.
                        """)

    args = parser.parse_args()

    if args.distance is not None:
        args.min_distance = args.max_distance = args.distance

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

    if args.R_vector is not None:
        args.R_vector = np.array(args.R_vector[:len(args.R_vector) // 3 * 3],
                                 dtype=int).reshape((len(args.R_vector)//3, 3))
        args.R_vector = list(map(tuple, args.R_vector.tolist()))

    if args.what_to_plot != 'all':
        plot_data_type = [args.what_to_plot]

    for data_type in plot_data_type:
        plot_2d(filename=args.filename,
                out_dir=args.output_dir,
                out_name=args.output_name,
                wtp=data_type,
                draw_cells=args.draw_cells,
                min_distance=args.min_distance,
                max_distance=args.max_distance,
                template=args.template,
                R_vector=args.R_vector,
                double_bonds=args.double_bonds,
                scale_atoms=args.scale_atoms,
                scale_data=args.scale_data,
                title=args.title)
