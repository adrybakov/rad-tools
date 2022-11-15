#! /usr/local/bin/python3

from argparse import ArgumentParser
from os import makedirs
from os.path import basename, abspath, join

from rad_tools.routines import OK, RESET

import matplotlib.pyplot as plt
import numpy as np

LABELS = {
    "s": ["s (up)", "s (down)"],
    "p": ["$p_z$ (up)", "$p_z$ (down)",
          "$p_y$ (up)", "$p_y$ (down)",
          "$p_x$ (up)", "$p_x$ (down)"],
    "d": ["$d_{z^2}$ (up)", "$d_{z^2}$",
          "$d_{zx}$ (up)", "$d_{zx}$ (down)",
          "$d_{zy}$ (up)", "$d_{zy}$ (down)",
          "$d_{x^2 - y^2}$ (up)", "$d_{x^2 - y^2}$ (down)",
          "$d_{xy}$ (up)", "$d_{xy}$ (down)"]
}

COLORS = {
    "s": ["blue", "red"],
    "p": ["#FF5C5C", "#6B0000",
          "#37FF5C", "#00A61E",
          "#2591FF", "#0061C4"],
    "d": ["#FF297A", "#6B0029",
          "#FFF821", "#898500",
          "#19FF19", "#00A600",
          "#10F9FF", "#00BFC4",
          "#5A08FF", "#4B00E1"]
}


def read_dos(filename, efermi=0):
    dos = np.loadtxt(filename, skiprows=1).T
    with open(filename, "r") as file:
        header = file.readline()
    dos[0] -= efermi
    if "up" in header and "dw" in header:
        for i in range(0, (len(dos) - 1) // 2):
            dos[i * 2 + 2] *= -1
    return dos


def decompose_filename(filename):
    atom_type = filename.split("(")[1].split(")")[0]
    atom_number = filename.split("#")[1].split("(")[0]
    orbital_type = filename.split("(")[2].split(")")[0]
    orbital_number = filename.split("#")[2].split("(")[0]
    return atom_type, orbital_type, atom_number, orbital_number


def manager(filename, out_dir, out_name, window=None, interactive=False):
    atom_type, orbital_type, atom_number, orbital_number = decompose_filename(
        basename(filename))
    dos = read_dos(filename)

    fig, ax = plt.subplots(figsize=(12, 4))

    ax.set_xlabel("Energy, eV")
    ax.set_ylabel("DOS, states/eV")
    ax.set_title(f"Atom: {atom_type}, orbitals: {orbital_type}")
    ax.vlines(0, 0, 1, transform=ax.get_xaxis_transform(),
              color='grey', linewidths=0.5,
              linestyles='dashed')

    ax.fill_between(dos[0], dos[2], dos[1],
                    color='grey', alpha=0.5, lw=0)
    # plot up
    for i in range(0, (len(dos) - 3) // 2):
        ax.plot(dos[0], dos[3 + 2 * i],
                '-',
                c=COLORS[orbital_type][2 * i],
                label=LABELS[orbital_type][2 * i],
                linewidth=0.7)
    # plot down
    for i in range(0, (len(dos) - 3) // 2):
        ax.plot(dos[0], dos[4 + 2 * i],
                '-',
                c=COLORS[orbital_type][2 * i + 1],
                label=LABELS[orbital_type][2 * i + 1],
                linewidth=0.7)

    ax.legend(loc='best', ncol=2)

    if window is not None:
        ax.set_xlim(*tuple(window))

    if interactive:
        plt.show()
    else:
        png_path = join(
            out_dir, f"{out_name}.#{atom_number}({atom_type})_wfc#{orbital_number}({orbital_type}).png")
        plt.savefig(png_path, dpi=400)
        print(f'{OK}Plot is in {abspath(png_path)}{RESET}')


if __name__ == "__main__":
    parser = ArgumentParser(description="Script for visualisation of projected density of states.",
                            epilog="""
               For the full description of arguments see the docs: 
               https://rad-tools.adrybakov.com/en/stable/user-guide/rad-dos-plotter.html
               """)
    parser.add_argument("-f", "--filename",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the dos file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder
                        for saving outputs.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='pdos',
                        help="""
                        Seedname for the output files.
                        """
                        )
    parser.add_argument("-w", "--window",
                        type=float,
                        nargs=2,
                        default=None,
                        help="""
                        Energy window.
                        """)
    parser.add_argument("-i", "--interactive",
                        action="store_true",
                        default=False,
                        help="""
                        Interactive mode flag.
                        """)

    args = parser.parse_args()

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

        manager(filename=args.filename,
                out_dir=args.output_dir,
                out_name=args.output_name,
                window=args.window,
                interactive=args.interactive)
