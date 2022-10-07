#! /usr/local/bin/python3
from argparse import ArgumentParser
from os.path import join, split, abspath
from math import atan, sqrt

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModelTB2J
from rad_tools.routines import check_make_dir, atom_mark_to_latex, rot_angle,\
    OK, RESET


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
            atoms=None,
            title=None):

    mode_name = "2d"
    messages = {'iso': "isotropic exchange",
                "distance": "distances"}

    model = ExchangeModelTB2J(filename)
    model = model.filter(min_distance=min_distance,
                         max_distance=max_distance,
                         R_vector=R_vector,
                         template=template)
    if atoms is None:
        atoms = model.magnetic_atoms

    dummy = True
    ha = 'right'
    if not double_bonds:
        model = model.remove_double_bonds()
        dummy = False
        ha = 'center'
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
                        fontsize=1.5 * fontsize * scale_atoms)
                ax.text(x2, y2, atom_mark_to_latex(atom2),
                        va='center', ha='center',
                        fontsize=1.5 * fontsize * scale_atoms)
                if wtp == 'iso':
                    ax.text(xm, ym, round(bond.iso, 4),
                            va='bottom', ha=ha,
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=fontsize * scale_data)
                elif wtp == 'distance':
                    ax.text(xm, ym, round(bond.dis, 4),
                            va='bottom', ha=ha,
                            rotation_mode='anchor',
                            rotation=rot_angle(x2 - x1, y2 - y1, dummy=dummy),
                            fontsize=fontsize * scale_data)

    if draw_cells:
        cells = model.get_cells()
        a_x, a_y, a_z = tuple(model.cell[0])
        b_x, b_y, b_z = tuple(model.cell[1])
        for Rx, Ry, Rz in cells:
            X_shift = Rx * a_x + Ry * b_x
            Y_shift = Rx * a_y + Ry * b_y
            ax.plot(np.array([0, a_x, a_x + b_x, b_x, 0]) + X_shift,
                    np.array([0, a_y, a_y + b_y, b_y, 0]) + Y_shift,
                    linewidth=1, color="#BCBF5A")

    if title is not None:
        ax.set_title(title, fontsize=1.5 * fontsize)

    png_path = join(out_dir, f'{out_name}.{wtp}.png')
    plt.savefig(png_path, dpi=400)
    print(f'{OK}{mode_name} plot with {wtp} is in {abspath(png_path)}{RESET}')


def get_molecule_center(atoms, exclude_atoms):
    filtered_atoms = {}
    for atom in atoms:
        keep = True
        for exclude_atom in exclude_atoms:
            if exclude_atom in atom:
                keep = False
                break
        if keep:
            filtered_atoms[atom] = atoms[atom]
    atoms = filtered_atoms
    x_min, y_min, z_min = None, None, None
    x_max, y_max, z_max = None, None, None
    for atom in atoms:
        x, y, z = atoms[atom]

        if x_min is None:
            x_min = x
        else:
            x_min = min(x, x_min)
        if y_min is None:
            y_min = y
        else:
            y_min = min(y, y_min)
        if z_min is None:
            z_min = z
        else:
            z_min = min(z, z_min)

        if x_max is None:
            x_max = x
        else:
            x_max = max(x, x_max)
        if y_max is None:
            y_max = y
        else:
            y_max = max(y, y_max)
        if z_max is None:
            z_max = z
        else:
            z_max = max(z, z_max)
            return (x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2


def get_distance(x1, y1, z1, x2, y2, z2):
    return sqrt((x2 - x1)**2 + (y2-y1)**2 + (z2-z1)**2)


def plot_molecule(filename, out_dir='.',
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
                  atoms=None,
                  title=None):

    mpl.rcParams.update(mpl.rcParamsDefault)
    mode_name = "molecule"
    messages = {'iso': "isotropic exchange",
                "distance": "distances"}

    model = ExchangeModelTB2J(filename)
    model = model.filter(min_distance=min_distance,
                         max_distance=max_distance,
                         R_vector=R_vector,
                         template=template)

    if atoms is None:
        atoms = model.magnetic_atoms

    fig, ax = plt.subplots()

    ax.set_xlabel('Distance to the molecule center, Angstroms')
    if wtp == 'iso':
        ax.set_ylabel('Isotropic exchange parameter, meV')
    elif wtp == 'distance':
        ax.set_ylabel('Bond length, Angstroms')

    data = [[], []]
    x, y, z = get_molecule_center(model._atoms, atoms)
    print(f"{OK}Molecule center is detected at "
          f"({round(x, 4)}, {round(y, 4)}, {round(z, 4)}){RESET}")

    for atom1 in model.bonds:
        for atom2 in model.bonds[atom1]:
            for R in model.bonds[atom1][atom2]:
                xa, ya, za = model.get_bond_coordinate(atom1, atom2, R)
                d = get_distance(x, y, z, xa, ya, za)
                bond = model.bonds[atom1][atom2][R]
                data[0].append(d)
                if wtp == 'iso':
                    data[1].append(bond.iso)
                elif wtp == 'distance':
                    data[1].append(bond.dis)
    data = np.array(data).T.tolist()

    data.sort(key=lambda x: x[0])
    data = np.array(data).T.tolist()
    ax.plot(data[0], data[1], "o-", color="black")
    if title is not None:
        ax.set_title(title, fontsize=15)

    png_path = join(out_dir, f'{out_name}.{mode_name}.{wtp}.png')
    plt.savefig(png_path, dpi=400)
    print(f'{OK}{mode_name} plot with {wtp} is in {abspath(png_path)}{RESET}')


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for visualisation of TB2J results",
                            epilog="""
                            See the docs: 
                            https://rad-tools.adrybakov.com/en/latest/tb2j_plotter.html
                            """)

    plot_data_type = ['iso', 'distance']
    plot_mode = {"2d": plot_2d, "molecule": plot_molecule}

    parser.add_argument("-f", "--file",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the *exchange.out* file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-m", "--mode",
                        type=str,
                        default="2d",
                        choices=["all", "2d", "molecule"],
                        help="""
                        Mode of plotting.

    Two modes are supported: structure with the view from above 
    and the plots with *value* over distance between bond and 
    the center of the molecule.
                        """
                        )
    parser.add_argument("-a", "--atoms",
                        type=str,
                        default=None,
                        nargs="*",
                        help="""
                        Atoms from the substrate

    Marks of atoms from the substracte (Same as in TB2J). 
    You can specify only names. For example instead of "Cr12" one can provide 
    "Cr" and then all Cr atoms will be thouth as a substrate ones.
                        """
                        )
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder for saving outputs.

    If the folder does not exist then it is created from the specified path.
    The creation is applied recursevly to the path, starting from the right
    until the existing folder is reached.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='exchange',
                        help="""
                        Seedname for the output files.

    Output files will have the following name structure:
    output-name.display_data_type.png
                        """
                        )
    parser.add_argument("-wtp", "--what-to-plot",
                        type=str,
                        choices=['all'] + plot_data_type,
                        default='all',
                        help="""
                        Type of data for display.

    Specifying the data for display at the graph. 
    Everything is displayed by default, each value in a separate picture. 
    Currently available for display: Isotropic exchange parameter, distance.
                        """
                        )
    parser.add_argument("-dc", "--draw-cells",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to draw the supercell`s shape.

    If specified then the shape of all supercells 
    presented in the model (after filtering) is drawn.
                        """
                        )
    parser.add_argument("-R", "--R-vector",
                        type=int,
                        nargs="*",
                        default=None,
                        help="""
                        R vectors for filtering the model.

    In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to). 
    Atom 1 is always located in (0, 0, 0) supercell, while atom 2 is located in 
    R = (i, j, k) supercell. This parameter tells the script to keep only the 
    bonds for which atom 2 is located in one of specified R supercells. 
    In order to specify supercells provide a set of integers separated 
    by spaces. They are grouped by three starting from the left to form a set 
    of R vectors. If the last group will contain 1 or 2 integers they will be 
    ignored.
                        """
                        )
    parser.add_argument("-maxd", "--max-distance",
                        type=float,
                        default=None,
                        help="""
                        (<=) Maximum distance.

    All the bonds with the distance beetwen atom 1 and atom 2 
    greater then maximum distance are excluded from the model.
                        """
                        )
    parser.add_argument("-mind", "--min-distance",
                        type=float,
                        default=None,
                        help="""
                        (>=) Minimum distance.

    All the bonds with the distance beetwen atom 1 and atom 2 
    lower then minimum distance are excluded from the model.
                        """
                        )
    parser.add_argument("-d", "--distance",
                        type=float,
                        default=None,
                        help="""
                        (=) Exact distance.

    Only the bonds with the exact distance remains in the model.
    Note: there is no point in specifying maximum or minimum distance when 
    this parameter is specified.
                        """
                        )
    parser.add_argument("-t", "--template",
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

    In TB2J file there are two bonds for the pair of atom 1 and atom 2: 
    from 1 to 2 and from 2 to 1 (when R = (0, 0, 0)). Isotropic and 
    anisotropic exchange and distance usially are exactly the same. 
    DMI vector have the same module and opposite directions. 
    If this parameter is specifyied then both bonds are displayed. 
    Otherwise bonds are combined in one by taking the average beetween
    exchange parameters (Note that it forces DMI to be equal to zero).
                        """)
    parser.add_argument("-sa", "--scale-atoms",
                        default=1,
                        type=float,
                        help="""
                        Scale for the size of atom marks.

    Use it if you want to display atom marks bigger or smaller. 
    Have to be positive.
                        """)
    parser.add_argument("-sd", "--scale-data",
                        default=1,
                        type=float,
                        help="""
                        Scale for the size of data text.

    Use it if you want to display data text marks bigger or smaller. 
    Have to be positive.
                        """)
    parser.add_argument("--title",
                        default=None,
                        type=str,
                        help="""
                        Title for the plots

    Title will be displayed in the picture.
                        """)

    args = parser.parse_args()

    if args.distance is not None:
        args.min_distance = args.max_distance = args.distance

    check_make_dir(args.output_dir)
    if args.R_vector is not None:
        args.R_vector = np.array(args.R_vector[:len(args.R_vector) // 3 * 3],
                                 dtype=int).reshape((len(args.R_vector)//3, 3))
        args.R_vector = list(map(tuple, args.R_vector.tolist()))

    if args.what_to_plot != 'all':
        plot_data_type = [args.what_to_plot]
    if args.mode != 'all':
        plot_mode = [plot_mode[key] for key in [args.mode]]
    else:
        plot_mode = [plot_mode[key] for key in plot_mode]

    for mode in plot_mode:
        for data_type in plot_data_type:
            mode(filename=args.file,
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
                 atoms=args.atoms,
                 title=args.title)
