#! /usr/local/bin/python3

from argparse import ArgumentParser
from os import makedirs, walk
from os.path import basename, abspath, join

from rad_tools.routines import OK, RESET, YELLOW

import matplotlib.pyplot as plt
import numpy as np

LABELS = {
    "s": ["s (up)", "s (down)"],
    "p": ["$p_z$ (up)", "$p_z$ (down)",
          "$p_y$ (up)", "$p_y$ (down)",
          "$p_x$ (up)", "$p_x$ (down)"],
    "d": ["$d_{z^2}$ (up)", "$d_{z^2}$ (down)",
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

CASES = ["collinear, spin-unpolarized",
         "collinear, spin-polarized",
         "non-collinear, non spin-orbit"]


def decompose_filenames(files):
    atoms = {}
    wfcs = {}
    for file in files:
        atom_number = int(file.split(".pdos_atm#")[1].split("(")[0])
        atom_symbol = file.split("(")[1].split(")")[0]
        wfc_number = int(file.split(")_wfc#")[1].split("(")[0])
        wfc_symbol = file.split(")_wfc#")[1].split("(")[1].split(")")[0]
        if atom_symbol not in atoms:
            atoms[atom_symbol] = []
        atoms[atom_symbol].append(atom_number)

        if (atom_number, atom_symbol) not in wfcs:
            wfcs[(atom_number, atom_symbol)] = []

        wfcs[(atom_number, atom_symbol)].append((wfc_number, wfc_symbol))

    for atom in atoms:
        atoms[atom].sort()
    for key in wfcs:
        wfcs[key].sort(key=lambda x: x[0])

    return atoms, wfcs


def detect_case(filename):
    with open(filename) as file:
        header = file.readline().lower()
    case = None
    if ("dos" in header
            and "pdos" in header):
        case = 0
    if ("dos" in header
        and "pdosup" in header
            and "pdosdw" in header):
        case = 2
    if ("dosup" in header
        and "dosdw" in header
        and "pdosup" in header
            and "pdosdw" in header):
        case = 1
    if case is None:
        raise RuntimeError("Unable to detect case, analized header:\n" +
                           f"{header}")
    return case


def analyse_input_folder(input_dir, filpdos=None):
    files = []
    for (dirpath, dirnames, filenames) in walk(input_dir):
        files.extend(filenames)
        break
    new_files = []
    if filpdos is not None:
        for file in files:
            if filpdos in file:
                new_files.append(file)
        files = new_files
        filpdos = [filpdos]
        files = [files]
    else:
        filpdos = set()
        for file in files:
            if ".pdos_tot" in file:
                filpdos.add(file.split(".pdos_tot")[0])
                new_files.append(file)
            elif ".pdos_atm#" in file:
                filpdos.add(file.split(".pdos_atm#")[0])
                new_files.append(file)
        filpdos = list(filpdos)

        print(f"Following DOS seednames (filpdos) are detected:{YELLOW}")
        for item in filpdos:
            print(f"   * {item}")
        print(RESET)

        if len(filpdos) != 1:
            files = [[] for _ in range(0, len(filpdos))]
            for file in new_files:
                for i, seedname in enumerate(filpdos):
                    if seedname in file:
                        files[i].append(file)
        else:
            files = [new_files]
    return filpdos, files


def plot_pdos_tot(case,
                  input_dir,
                  seedname,
                  output_dir=".",
                  window=None,
                  efermi=0.):
    dos = np.loadtxt(join(input_dir, f"{seedname}.pdos_tot"),
                     skiprows=1).T
    dos[0] -= efermi

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.set_xlabel("Energy, eV", fontsize=18)
    ax.set_ylabel("DOS, states/eV", fontsize=18)
    ax.set_title(f"DOS and PDOS", fontsize=18)
    ax.vlines(0, 0, 1, transform=ax.get_xaxis_transform(),
              color='grey', linewidths=0.5,
              linestyles='dashed')
    if window is not None:
        dos = dos.T
        i_min, i_max = 0, 0
        for i in range(0, len(dos)):
            if dos[i][0] < window[0]:
                i_min = i
            if dos[i][0] < window[1]:
                i_max = i
        dos = dos[i_min:i_max].T
        ax.set_xlim(*tuple(window))

    if case == 0:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='grey',
                        alpha=0.5, label="DOS")
        ax.plot(dos[0], dos[2], "-",
                lw=1, color="orange",
                alpha=0.5, label="PDOS")
        ncol = 1
    if case == 1:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='blue',
                        alpha=0.5, label="DOS (up)")
        ax.plot(dos[0], dos[3], "-",
                lw=1, color="orange",
                alpha=0.5, label="PDOS (up)")
        ax.fill_between(dos[0], 0, -dos[2],
                        lw=0, color='red',
                        alpha=0.5,  label="DOS (down)")
        ax.plot(dos[0], -dos[4], "-",
                lw=1, color="darkgreen",
                alpha=0.5, label="PDOS (down)")
        ncol = 2
    if case == 2:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='grey',
                        alpha=0.5, label="DOS")
        ax.plot(dos[0], dos[2], "-",
                lw=1, color="orange",
                alpha=0.5, label="PDOS (up)")
        ax.plot(dos[0], -dos[3], "-",
                lw=1, color="darkgreen",
                alpha=0.5, label="PDOS (down)")
        ncol = 1

    ax.legend(loc="best", ncol=ncol)

    png_path = join(
        output_dir, seedname, f"{seedname}.pdos_tot.png")
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    print(f"{YELLOW}Total dos plot for seedname {seedname} is in {abspath(png_path)}{RESET}")
    plt.close()


def manager(input_dir,
            filpdos=None,
            output_dir=".",
            window=None,
            efermi=0.,
            separate=False):

    try:
        makedirs(output_dir)
    except FileExistsError:
        pass

    filpdos, files = analyse_input_folder(input_dir, filpdos)
    for s_i, seedname in enumerate(filpdos):
        print(f"({s_i + 1}/{len(filpdos)}) " +
              f"Start to work with {seedname} seedname")
        case = detect_case(join(input_dir, f"{seedname}.pdos_tot"))
        print(f"{CASES[case]} case detected.")

        try:
            makedirs(join(output_dir, seedname))
        except FileExistsError:
            pass

        atoms, wfcs = decompose_filenames(
            [file for file in files[s_i] if ".pdos_tot" not in file])

        plot_pdos_tot(case,
                      input_dir,
                      seedname,
                      output_dir=output_dir,
                      window=window,
                      efermi=efermi)

        print(f"{OK}" +
              f"Finish to work with {seedname} seedname, results are in \n" +
              f"  {abspath(join(output_dir, seedname))}" +
              f"{RESET}")

    # ax.fill_between(dos[0], dos[2], dos[1],
    #                 color='grey', alpha=0.5, lw=0)
    # if window is not None:
    #     ax.set_xlim(*tuple(window))


def get_parser():
    parser = ArgumentParser(
        description="Script for visualisation of density of states.")
    parser.add_argument("-id", "--input-dir",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the folder with dos files.
                        """
                        )
    parser.add_argument("-f", "--filpdos",
                        type=str,
                        default=None,
                        help="""
                        Prefix for output files containing PDOS(E). 
                        As specified in the QE projwfc.x input file.
                        """)
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder
                        for saving outputs.
                        """
                        )
    parser.add_argument("-w", "--window",
                        type=float,
                        nargs=2,
                        default=None,
                        help="""
                        Energy window for the plots.  
                        By default whole range present in the files is plotted
                        """)
    parser.add_argument("-ef", "--efermi",
                        type=float,
                        default=0.,
                        help="""
                        Fermi energy, zero by default.
                        """)
    parser.add_argument("-s", "--separate",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to plot projected DOS for each atom 
                        of the same type separately.
                        """)

    return parser
