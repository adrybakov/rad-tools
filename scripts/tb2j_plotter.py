#! /usr/local/bin/python3
from argparse import ArgumentParser
from os.path import join, split
from math import atan

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModelTB2J
from rad_tools.routines import check_make_dir, atom_mark_to_latex, rot_angle


def plot_2d(filename, out_dir, wtp='iso', draw_cells=False,
            min_distance=None, max_distance=None,
            template=None, R_vector=None,
            double_bonds=False):

    model = ExchangeModelTB2J(filename)
    model = model.filter(min_distance=min_distance,
                         max_distance=max_distance,
                         R_vector=R_vector,
                         template=template)
    dummy = False
    if not double_bonds:
        model = model.remove_double_bonds()
        dummy = True
    X, Y, Z = model.get_space_dimensions()
    if X == 0 and Y == 0:
        X = Y = 1
    fontsize = 10 * 1.1 * Y / 5
    plt.rcParams.update({'font.size': fontsize})
    mpl.rcParams.update({'axes.linewidth': 1.1 * Y / 5})
    mpl.rcParams.update({'xtick.major.width': 1.1 * Y / 5})
    mpl.rcParams.update({'ytick.major.width': 1.1 * Y / 5})
    mpl.rcParams.update({'xtick.major.size': 5 * 1.1 * Y / 5})
    mpl.rcParams.update({'ytick.major.size': 5 * 1.1 * Y / 5})
    mpl.rcParams.update({'xtick.major.pad': 5 * 1.1 * Y / 5})
    mpl.rcParams.update({'ytick.major.pad': 5 * 1.1 * Y / 5})
    mpl.rcParams.update({'axes.labelpad': 1.1 * Y / 5})

    fig = plt.figure(figsize=(1.1 * X, 1.1 * Y))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlabel('x, Angstroms')
    ax.set_ylabel('y, Angstroms')
    label_ax = fig.add_axes([0.85, 0.1, 0.1, 0.8])
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(0, 1)
    label_ax.axis('off')

    for atom1 in model.bonds:
        for atom2 in model.bonds[atom1]:
            for R in model.bonds[atom1][atom2]:
                bond = model.bonds[atom1][atom2][R]
                x1, y1, z1, x2, y2, z2 = model.get_atom_coordinates(atom1,
                                                                    atom2,
                                                                    R)
                xm = (x1 + x2) / 2
                ym = (y1 + y2) / 2
                zm = (z1 + z2) / 2

                ax.scatter(x1, y1, s=100 * fig.dpi/72., c='white')
                ax.scatter(x2, y2, s=100 * fig.dpi/72., c='white')

                ax.text(x1, y1, atom_mark_to_latex(atom1),
                        va='center', ha='center',
                        fontsize=15)
                ax.text(x2, y2, atom_mark_to_latex(atom2),
                        va='center', ha='center',
                        fontsize=15)
                if wtp == 'iso':
                    ax.text(xm, ym, bond.iso,
                            va='bottom', ha='right',
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=10)
                elif wtp == 'distance':
                    ax.text(xm, ym, bond.dis,
                            va='bottom', ha='right',
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1),
                            fontsize=10)

    if draw_cells:
        cells = model.get_cells()
        a_x, a_y, a_z = tuple(model.cell[0])
        b_x, b_y, b_z = tuple(model.cell[1])
        for Rx, Ry, Rz in cells:
            Rx = Rx * (a_x + b_x)
            Ry = Ry * (a_y + b_y)
            ax.plot(np.array([0, a_x, a_x + b_x, b_x, 0]) + Rx,
                    np.array([0, a_y, a_y + b_y, b_y, 0]) + Ry,
                    linewidth=1, color="#BCBF5A")

    plt.savefig(join(out_dir, f'{split(filename)[1]}.png'), dpi=400)


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for visualisation of "
                            "TB2J results")

    parser.add_argument("-f", "--file",
                        type=str,
                        required=True,
                        help="TB2J *.out file")
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="Directory for saving outputs "
                        "if there will be one. "
                        "Could be non-existing, "
                        "but the parent directory have to exist")
    parser.add_argument("-wtp", "--what-to-plot",
                        type=str,
                        choices=['iso', 'distance'],
                        default='iso',
                        help="What to plot? Isotropic exchange by default")
    parser.add_argument("-dc", "--draw-cells",
                        action="store_true",
                        default=False,
                        help="Whenever to draw the supercells. "
                        "False by default")
    parser.add_argument("-R", "--R-vector",
                        type=int,
                        nargs="*",
                        default=None,
                        help="R vectors for filtering the model. "
                        "Specify 3n integers, they will be split into vectors. "
                        "(First 3 -> R_1, ...) "
                        "If 3n+1 or 3n+2 integers are specified then "
                        "the last 1 or 2 numbers will be ignored")
    parser.add_argument("-maxd", "--max-distance",
                        type=float,
                        default=None,
                        help="Maximum distance for the neighbors to be shown "
                        "(<=)"
                        )
    parser.add_argument("-mind", "--min-distance",
                        type=float,
                        default=None,
                        help="Minimum distance for the neighbors to be shown "
                        "(>=)"
                        )
    parser.add_argument("-d", "--distance",
                        type=float,
                        default=None,
                        help="Exact distance for the neighbors to be shown "
                        "(=)"
                        )
    parser.add_argument("-t", "--template",
                        type=str,
                        default=None,
                        help="Template for filtering the Exchange, "
                        "it have to be a plain text file "
                        "which will be passed to the Template class")
    parser.add_argument("-db", "--double-bonds",
                        default=False,
                        action="store_true",
                        help="Whenever to show both directions "
                        "of the bond (if they are in TB2J file at first place)")

    args = parser.parse_args()

    if args.distance is not None:
        args.min_distance = args.max_distance = args.distance

    check_make_dir(args.output_dir)
    if args.R_vector is not None:
        args.R_vector = np.array(args.R_vector[:len(args.R_vector) // 3 * 3],
                                 dtype=int).reshape((len(args.R_vector)//3, 3))
        args.R_vector = list(map(tuple, args.R_vector.tolist()))

    plot_2d(filename=args.file,
            out_dir=args.output_dir,
            wtp=args.what_to_plot,
            draw_cells=args.draw_cells,
            min_distance=args.min_distance,
            max_distance=args.max_distance,
            template=args.template,
            R_vector=args.R_vector,
            double_bonds=args.double_bonds)
