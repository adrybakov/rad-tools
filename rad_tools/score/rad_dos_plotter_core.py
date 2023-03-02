#! /usr/local/bin/python3

from argparse import ArgumentParser
from os import makedirs, walk
from os.path import abspath, join
from copy import deepcopy

from rad_tools.routines import OK, RESET, YELLOW

import matplotlib.pyplot as plt
import numpy as np

CASES = ["collinear, spin-unpolarized",
         "collinear, spin-polarized",
         "non-collinear, non spin-orbit"]


def analyse_input_folder(input_dir, filpdos=None):
    r"""
    Analyse input folder, detects seednames for the dos output files.

    Parameters
    ----------
    input_dir : str
        Directory with DOS files.
    filpdos : str, default None
        Seddname for the DOS file. If not provided then all seednames
        from the folder are detected and processed.

    Returns
    -------
    filpdos : list
        List of seednames found in ``input_dir``.
    filenames : list
        List of lists of filenames for each seedname in ``filpdos``.
    """

    files = []
    for (dirpath, dirnames, filenames) in walk(input_dir):
        files.extend(filenames)
        break

    filenames = []
    if filpdos is not None:
        for filename in files:
            if filpdos in filename:
                filenames.append(filename)
        files = filenames
        filpdos = [filpdos]
        filenames = [filenames]
    else:
        filpdos = set()
        for filename in files:
            if ".pdos_tot" in filename:
                filpdos.add(filename.split(".pdos_tot")[0])
                filenames.append(filename)
            elif ".pdos_atm#" in filename:
                filpdos.add(filename.split(".pdos_atm#")[0])
                filenames.append(filename)
        filpdos = list(filpdos)

        print(f"Following DOS seednames (filpdos) are detected:{YELLOW}")
        for item in filpdos:
            print(f"   * {item}")
        print(RESET)

        if len(filpdos) != 1:
            tmp = [[] for _ in range(0, len(filpdos))]
            for filename in filenames:
                for i, seedname in enumerate(filpdos):
                    if seedname in filename:
                        tmp[i].append(filename)
            filenames = tmp
        else:
            filenames = [filenames]
    return filpdos, filenames


def detect_case(filename):
    r"""
    Detects case of the DOS calculations.

    Parameters
    ----------
    filename: str
        Path to the file in which the calculation case have to be detected

    Returns
    -------
    case : int
        Case index.
    """

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


def decompose_filenames(filenames):
    r"""
    Decompose filenames and extract information about atoms and projections.

    Parameters
    ----------
    filenames : str
        List of filenames with projected DOS (not pdos_total).

    Returns
    -------
    atoms : dist
        Distionary with the atom labels and their numbers.

        .. code-block:: python

            atoms = {atom1: [n_1, n_2, n_3, ...]}
    wfcs : dict
        Dictionary of projection functions and their numbers.

        .. code-block:: python

            wfcs = {(atom1, atom1_number): (wfc_label, wfc_number)}
    """

    atoms = {}
    wfcs = {}
    for filename in filenames:

        # Detect names and numbers
        atom_number = int(filename.split(".pdos_atm#")[1].split("(")[0])
        atom_symbol = filename.split("(")[1].split(")")[0]
        wfc_number = int(filename.split(")_wfc#")[1].split("(")[0])
        wfc_symbol = filename.split(")_wfc#")[1].split("(")[1].split(")")[0]

        if atom_symbol not in atoms:
            atoms[atom_symbol] = []
        atoms[atom_symbol].append(atom_number)

        if (atom_symbol, atom_number) not in wfcs:
            wfcs[(atom_symbol, atom_number)] = []
        wfcs[(atom_symbol, atom_number)].append((wfc_symbol, wfc_number))

    # Sort entries
    for atom in atoms:
        atoms[atom] = list(set(atoms[atom]))
        atoms[atom].sort()
    for key in wfcs:
        wfcs[key].sort(key=lambda x: x[1])

    return atoms, wfcs


def combine_by_atoms(input_dir,
                     seedname,
                     output_dir,
                     atoms,
                     wfcs,
                     case):
    r"""
    Combine PDOS by atom type.

    Parameters
    ----------
    input_dir : str
    seedname : str
    output_dir : str
    atoms : dict
    wfcs : dict
    case : int
    """
    dos_by_atom = {}
    print(f"Start to combine PDOS by atoms")
    for atom in atoms:
        print(f"    {len(atoms[atom])} atoms of {atom} detected")
        for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
            with open(join(input_dir,
                           f"{seedname}.pdos_atm#{atoms[atom][0]}({atom})_wfc#{wfc_number}({wfc_symbol})")) as file:
                header = file.readline().replace("\n", "")
            for i in atoms[atom]:
                if i == atoms[atom][0]:
                    dos = np.loadtxt(join(input_dir,
                                          f"{seedname}.pdos_atm#{i}({atom})_wfc#{wfc_number}({wfc_symbol})"),
                                     skiprows=1).T
                else:
                    dos += np.loadtxt(join(input_dir,
                                           f"{seedname}.pdos_atm#{i}({atom})_wfc#{wfc_number}({wfc_symbol})"),
                                      skiprows=1).T
            dos[0] = dos[0] / len(atoms[atom])
            energy = dos[0]
            if case == 0:
                dos_by_atom[f"{atom}, {wfc_symbol}"] = dos[1]
            elif case == 1 or case == 2:
                dos_by_atom[f"{atom}, {wfc_symbol}"] = [dos[1], dos[2]]

            np.savetxt(join(output_dir,
                            seedname,
                            "summed-by-atom",
                            f"{atom}_wfc#{wfc_number}({wfc_symbol})"),
                       dos.T,
                       header=header,
                       comments="",
                       fmt="%8.3f  "+"%.3E  "*(dos.shape[0]-1))
    print(f"Files with PDOS combined by atom are in " +
          f"{abspath(join(output_dir,seedname,'summed-by-atom'))}")


def filter_window(dos, window, efermi=0.):
    r"""
    Filter dos with respect to the window and make a shift to Fermi energy.

    Parameters
    ----------
    dos : array
        Density of states array.
    window : tuple, default None
        Limits for the energy.
    efermi : float, default 0
        Fermi energy.
    """

    dos[0] -= efermi
    if window is not None:
        i_min, i_max = 0, 0
        for i in range(0, len(dos[0])):
            if dos[0][i] < window[0]:
                i_min = i
            if dos[0][i] < window[1]:
                i_max = i
        dos = dos.T[i_min:i_max].T
    return dos


def plot_pdos_tot(dos, output_name, case, window=None, efermi=0.):

    dos = filter_window(dos, window=window, efermi=efermi)

    ax = plt.subplots(figsize=(8, 4))[1]

    ax.set_xlabel("Energy, eV", fontsize=18)
    ax.set_ylabel("DOS, states/eV", fontsize=18)
    if efermi != 0:
        ax.set_title(f"DOS and PDOS (0 is Fermi energy)", fontsize=18)
    else:
        ax.set_title(
            f"DOS and PDOS (0 is 0)", fontsize=18)
    ax.vlines(0, 0, 1, transform=ax.get_xaxis_transform(),
              color='grey', linewidths=0.5,
              linestyles='dashed')

    if case == 0:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='grey',
                        alpha=0.3, label="DOS")
        ax.plot(dos[0], dos[2], "-",
                lw=1, color="orange",
                alpha=0.7, label="PDOS")
        ncol = 1
    if case == 1:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='blue',
                        alpha=0.3, label="DOS (up)")
        ax.plot(dos[0], dos[3], "-",
                lw=1, color="blue",
                alpha=0.7, label="PDOS (up)")
        ax.fill_between(dos[0], 0, -dos[2],
                        lw=0, color='red',
                        alpha=0.3,  label="DOS (down)")
        ax.plot(dos[0], -dos[4], "-",
                lw=1, color="red",
                alpha=0.7, label="PDOS (down)")
        ncol = 2
    if case == 2:
        ax.fill_between(dos[0], 0, dos[1],
                        lw=0, color='grey',
                        alpha=0.3, label="DOS")
        ax.plot(dos[0], dos[2], "-",
                lw=1, color="blue",
                alpha=0.7, label="PDOS (up)")
        ax.plot(dos[0], -dos[3], "-",
                lw=1, color="red",
                alpha=0.7, label="PDOS (down)")
        ncol = 1

    ax.legend(loc="best", ncol=ncol)

    png_path = f"{output_name}.png"
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    print(f"Total dos plot  is in {abspath(png_path)}")
    plt.close()


def plot_projected(l, dos, output_name, title, case,
                   window=None, efermi=0.):
    labels = {"s": ["s"],
              "p": ["$p_z$", "$p_y$", "$p_x$"],
              "d": ["$d_{z^2}$", "$d_{zx}$",
                    "$d_{zy}$", "$d_{x^2 - y^2}$", "$d_{xy}$"]}
    factor = {"s": 1,
              "p": 3,
              "d": 5}
    dos = filter_window(dos, window=window, efermi=efermi)
    if window is None:
        window = (np.amin(dos[0]), np.amax(dos[0]))
    fig, axs = plt.subplots(factor[l], 1, figsize=(9, factor[l]*2))
    if l == "s":
        axs = [axs]
    fig.subplots_adjust(hspace=0)
    for i in range(0, factor[l]):
        ax = axs[i]
        ax.set_ylabel('DOS, states/eV', fontsize=12)
        if i == factor[l] - 1:
            ax.set_xlabel('E, ev', fontsize=15)
        else:
            ax.axes.get_xaxis().set_visible(False)
        ax.set_xlim(*tuple(window))
        ax.vlines(0, 0, 1, transform=ax.get_xaxis_transform(),
                  color='grey', linewidths=0.5, linestyles='dashed')
        if title is not None and i == 0:
            ax.set_title(title)

        if case == 0:
            ax.fill_between(dos[0], 0, dos[1],
                            lw=0, color='black',
                            alpha=0.3, label="p")
            ax.plot(dos[0], dos[2+i],
                    "-",
                    lw=0.5,
                    color="black",
                    alpha=0.8,
                    label=labels[l][i])
        elif case in [1, 2]:
            ax.fill_between(dos[0], 0, dos[1],
                            lw=0, color='blue',
                            alpha=0.3, label=f"{l} (up)")
            ax.plot(dos[0], dos[3 + 2*i],
                    "-",
                    lw=0.5,
                    color="blue",
                    alpha=0.8,
                    label=f"{labels[l][i]} (up)")
            ax.fill_between(dos[0], 0, -dos[2],
                            lw=0, color='red',
                            alpha=0.3, label=f"{l} (down)")
            ax.plot(dos[0], -dos[4 + 2*i],
                    "-",
                    lw=0.5,
                    color="red",
                    alpha=0.8,
                    label=f"{labels[l][i]} (down)")
        ax.legend(loc=(1.025, 0.2),
                  bbox_transform=ax.transAxes)

    plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")


def plot_projected_relative(l, dos, output_name, title, case,
                            window=None, efermi=0.,
                            normalize=False):
    labels = {"s": ["s"],
              "p": ["$p_z$", "$p_y$", "$p_x$"],
              "d": ["$d_{z^2}$", "$d_{zx}$",
                    "$d_{zy}$", "$d_{x^2 - y^2}$", "$d_{xy}$"]}
    factor = {"s": 1,
              "p": 3,
              "d": 5}
    colors = ["#3C6FD8", "#55A822", "#8F30CE", "#B67E23", "#C22661"]
    dos = filter_window(dos, window=window, efermi=efermi)
    if window is None:
        window = (np.amin(dos[0]), np.amax(dos[0]))
    fig, axs = plt.subplots(2, 1, figsize=(9, 4))
    ax1 = axs[0]
    ax2 = axs[1]
    fig.subplots_adjust(hspace=0)
    if normalize:
        ax1.set_ylabel('PDOS / LDOS', fontsize=12)
        ax2.set_ylabel('PDOS / LDOS', fontsize=12)
    else:
        ax1.set_ylabel('DOS, states/eV', fontsize=12)
        ax2.set_ylabel('DOS, states/eV', fontsize=12)
    ax1.axes.get_xaxis().set_visible(False)
    ax2.set_xlabel('E, ev', fontsize=15)
    ax1.set_xlim(*tuple(window))
    ax2.set_xlim(*tuple(window))
    ax1.vlines(0, 0, 1, transform=ax1.get_xaxis_transform(),
               color='grey', linewidths=0.5, linestyles='dashed')
    ax2.vlines(0, 0, 1, transform=ax2.get_xaxis_transform(),
               color='grey', linewidths=0.5, linestyles='dashed')
    if title is not None:
        ax1.set_title(title)

    y_up = np.zeros(dos[0].shape)
    y_prev_up = np.zeros(dos[0].shape)
    y_down = np.zeros(dos[0].shape)
    y_prev_down = np.zeros(dos[0].shape)
    if case == 0:
        y_up = np.zeros(dos[0].shape)
        y_prev_up = np.zeros(dos[0].shape)
        y_down = np.zeros(dos[0].shape)
        y_prev_down = np.zeros(dos[0].shape)
    else:
        y = np.zeros(dos[0].shape)
        y_prev = np.zeros(dos[0].shape)

    for i in range(0, factor[l]):
        if case == 0:
            y += dos[2+i]
            if normalize:
                ax1.fill_between(dos[0], y_prev/dos[1], y/dos[1],
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
            else:
                ax1.fill_between(dos[0], y_prev, y,
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
            y_prev = deepcopy(y)
        elif case in [1, 2]:
            y_up += dos[3+2*i]
            y_down += dos[4+2*i]
            if normalize:
                ax1.fill_between(dos[0], y_prev_up/dos[1], y_up/dos[1],
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
                ax2.fill_between(dos[0], y_prev_down/dos[2], y_down/dos[2],
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
            else:
                ax1.fill_between(dos[0], y_prev_up, y_up,
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
                ax2.fill_between(dos[0], y_prev_down, y_down,
                                 lw=0, color=colors[i],
                                 label=f"{labels[l][i]}")
            y_prev_up = deepcopy(y_up)
            y_prev_down = deepcopy(y_down)

    if normalize:
        ax1.set_ylim(0, 1)
        ax2.set_ylim(1, 0)
    else:
        ax1.set_ylim(0, None)
        ax2.set_ylim(0, None)
        ax2.invert_yaxis()
    ax1.legend(loc=(1.025, 0.2),
               bbox_transform=ax1.transAxes, title="up", reverse=True)
    ax2.legend(loc=(1.025, 0.2),
               bbox_transform=ax2.transAxes, title="down")

    plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")


def manager(input_dir,
            filpdos=None,
            output_dir=".",
            window=None,
            efermi=0.,
            separate=False,
            relative=False,
            normalize=False):

    try:
        makedirs(output_dir)
    except FileExistsError:
        pass

    filpdos, files = analyse_input_folder(input_dir, filpdos)

    for s_i, seedname in enumerate(filpdos):

        # Preparations
        print(f"{YELLOW}({s_i + 1}/{len(filpdos)}) " +
              f"Start to work with {seedname} seedname{RESET}")

        case = detect_case(join(input_dir, f"{seedname}.pdos_tot"))
        print(f"{CASES[case]} case detected.")

        try:
            makedirs(join(output_dir, seedname))
        except FileExistsError:
            pass
        try:
            makedirs(join(output_dir, seedname, "summed-by-atom"))
        except FileExistsError:
            pass
        if separate:
            try:
                makedirs(join(output_dir, seedname, "individual-plots"))
            except FileExistsError:
                pass

        atoms, wfcs = decompose_filenames(
            [file for file in files[s_i] if ".pdos_tot" not in file])

        # Plot PDOS vs DOS
        total_dos = np.loadtxt(join(input_dir, f"{seedname}.pdos_tot"),
                               skiprows=1).T
        plot_pdos_tot(total_dos,
                      join(output_dir, seedname, "pdos-vs-dos"),
                      case,
                      window=window,
                      efermi=efermi)

        # Summ PDOS for the same atom types
        combine_by_atoms(input_dir,
                         seedname,
                         output_dir,
                         atoms,
                         wfcs,
                         case)

        # Plot summed PDOS
        for atom in atoms:
            for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                path = join(input_dir,
                            seedname,
                            "summed-by-atom",
                            f"{atom}_wfc#{wfc_number}({wfc_symbol})")
                dos = np.loadtxt(path, skiprows=1).T
                if efermi == 0:
                    title = f"PDOS for {atom} ({wfc_symbol}) (0 is 0)"
                else:
                    title = f"PDOS for {atom} ({wfc_symbol}) (0 is Fermi energy)"
                if relative:
                    plot_projected_relative(wfc_symbol,
                                            dos,
                                            path,
                                            title,
                                            case,
                                            window=window,
                                            efermi=efermi,
                                            normalize=normalize)
                else:
                    plot_projected(wfc_symbol,
                                   dos,
                                   path,
                                   title,
                                   case,
                                   window=window,
                                   efermi=efermi)
        # Plot all individual plots
        if separate:
            for atom in atoms:
                for a_number in atoms[atom]:
                    for wfc_symbol, wfc_number in wfcs[(atom, a_number)]:
                        dos = np.loadtxt(join(input_dir,
                                              f"{seedname}.pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})"),
                                         skiprows=1).T
                        if efermi == 0:
                            title = f"PDOS for {atom} #{a_number} ({wfc_symbol}) (0 is 0)"
                        else:
                            title = f"PDOS for {atom} #{a_number} ({wfc_symbol}) (0 is Fermi energy)"
                        plot_projected(wfc_symbol,
                                       dos,
                                       join(output_dir, seedname,
                                            "individual-plots",
                                            f"pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})"),
                                       title,
                                       case,
                                       window=window,
                                       efermi=efermi)

        print(f"{OK}" +
              f"Finish to work with {seedname} seedname, results are in " +
              f"  {abspath(join(output_dir, seedname))}" +
              f"{RESET}")


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
    parser.add_argument("-r", "--relative",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to use relative style.
                        """)
    parser.add_argument("-n", "--normalize",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to use normalize relative style.
                        """)

    return parser
