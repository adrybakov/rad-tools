#! /usr/local/bin/python3

import re
from argparse import ArgumentParser
from copy import deepcopy
from os import makedirs, walk
from os.path import abspath, join

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import __version__ as matplotlib_version

from rad_tools.routines import OK, RESET, YELLOW

CASES = [
    "collinear, spin-unpolarized",
    "collinear, spin-polarized",
    "non-collinear, non spin-orbit",
]


class DOS:
    def __init__(self, filename=None) -> None:
        self.filename = filename
        self.dos = None
        self.case = None

    @property
    def energy(self):
        if self.dos is not None:
            return self.dos[0]
        return None

    @property
    def dos(self):
        if self.dos is not None:
            return self.dos[1:]
        return None

    def __add__(self, other):
        if not isinstance(other, DOS):
            raise TypeError

    def read_dos(self):
        self.dos = np.loadtxt(self.filename, skiprows=1).T


def analyse_input_folder(input_path, filpdos=None):
    r"""
    Analyze input folder, detects seednames for the dos output files.

    Parameters
    ----------
    input_path : str
        Directory with DOS files.
    filpdos : str, default None
        Seedname for the DOS file. If not provided then all seednames
        from the folder are detected and processed.

    Returns
    -------
    filpdos : list
        List of seednames found in ``input_path``.
    filenames : list
        List of lists of filenames for each seedname in ``filpdos``.
    """

    pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf_0-9j.]*\\)"

    # Get list of files in the folder
    files = []
    for dirpath, dirnames, filenames in walk(input_path):
        files.extend(filenames)
        break

    filenames = []
    # Search for provided filpdos
    if filpdos is not None:
        for filename in files:
            if re.fullmatch(f"{re.escape(filpdos)}.pdos_tot", filename) or re.fullmatch(
                f"{re.escape(filpdos)}{pattern}",
                filename,
            ):
                filenames.append(filename)
        files = filenames
        filpdos = [filpdos]
        filenames = [filenames]
    # Try to detect all filpdos
    else:
        filpdos = set()
        for filename in files:
            if ".pdos_tot" in filename and ".pdos_tot" == filename[-9:]:
                filpdos.add(filename[:-9])
                filenames.append(filename)
            elif re.match(f".*{pattern}$", filename):
                filpdos.add(re.split(f"{pattern}$", filename)[0])
                filenames.append(filename)
        filpdos = list(filpdos)

        print(f"Following DOS seednames (filpdos) are detected:{YELLOW}")
        for item in filpdos:
            print(f"   * {item}")
        print(RESET)

        # Group filpdos with corresponding filenames
        if len(filpdos) != 1:
            tmp = [[] for _ in range(0, len(filpdos))]
            for filename in filenames:
                for i, seedname in enumerate(filpdos):
                    if re.fullmatch(
                        f"{re.escape(seedname)}.pdos_tot", filename
                    ) or re.fullmatch(
                        f"{re.escape(seedname)}{pattern}",
                        filename,
                    ):
                        tmp[i].append(filename)
            filenames = tmp
        else:
            filenames = [filenames]
    return filpdos, filenames


def load_dos(filename):
    dos = np.loadtxt(filename, skiprows=1).T
    with open(filename, "r") as file:
        header = file.readline()
    if "ik" in header:
        nkpoints = int(dos[0][-1])
        ndpoints = len(dos[0])
        dos = dos[1:].reshape((len(dos) - 1, nkpoints, ndpoints // nkpoints))
        dos = np.sum(dos, axis=(1))
        dos /= nkpoints
    return dos


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
        header = file.readline().lower().split()
    case = None
    if "dos(e)" in header and "pdos(e)" in header:
        case = 0
    if "dos(e)" in header and "pdosup(e)" in header and "pdosdw(e)" in header:
        case = 2
    if (
        "dosup(e)" in header
        and "dosdw(e)" in header
        and "pdosup(e)" in header
        and "pdosdw(e)" in header
    ):
        case = 1
    if case is None:
        raise RuntimeError("Unable to detect case, analysed header:\n" + f"{header}")
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
        Dictionary with the atom labels and their numbers.

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


def combine_by_atoms(input_path, seedname, output_path, atoms, wfcs, case, suffix):
    r"""
    Combine PDOS by atom type.

    Parameters
    ----------
    input_path : str
    seedname : str
    output_path : str
    atoms : dict
    wfcs : dict
    case : int
    """
    dos_by_atom = {}
    print(f"Start to combine PDOS by atoms")
    for atom in atoms:
        print(f"    {len(atoms[atom])} atoms of {atom} detected")
        for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
            with open(
                join(
                    input_path,
                    f"{seedname}.pdos_atm#{atoms[atom][0]}({atom})_wfc#{wfc_number}({wfc_symbol})",
                )
            ) as file:
                header = file.readline().replace("\n", "")
                if "ik" in header:
                    header = header.replace("ik", "")
            for i in atoms[atom]:
                if i == atoms[atom][0]:
                    dos = load_dos(
                        join(
                            input_path,
                            f"{seedname}.pdos_atm#{i}({atom})_wfc#{wfc_number}({wfc_symbol})",
                        )
                    )
                else:
                    dos += load_dos(
                        join(
                            input_path,
                            f"{seedname}.pdos_atm#{i}({atom})_wfc#{wfc_number}({wfc_symbol})",
                        )
                    )
            dos[0] = dos[0] / len(atoms[atom])
            energy = dos[0]
            if case == 0:
                dos_by_atom[f"{atom}, {wfc_symbol}"] = dos[1]
            elif case == 1 or case == 2:
                dos_by_atom[f"{atom}, {wfc_symbol}"] = [dos[1], dos[2]]

            np.savetxt(
                join(
                    output_path,
                    f"{seedname}{suffix}",
                    "summed-by-atom",
                    f"{atom}_wfc#{wfc_number}({wfc_symbol})",
                ),
                dos.T,
                header=header,
                comments="",
                fmt="%8.3f  " + "%.3E  " * (dos.shape[0] - 1),
            )
    print(
        f"Files with PDOS combined by atom are in "
        + f"{abspath(join(output_path, f'{seedname}{suffix}', 'summed-by-atom'))}"
    )


def filter_window(dos, window, efermi=0.0):
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


def plot_pdos_tot(dos, output_name, case, window=None, efermi=0.0):
    dos = filter_window(dos, window=window, efermi=efermi)
    if window is None:
        window = (np.amin(dos[0]), np.amax(dos[0]))

    ax = plt.subplots(figsize=(8, 4))[1]

    ax.set_xlabel("Energy, eV", fontsize=18)
    ax.set_ylabel("DOS, states/eV", fontsize=18)
    if efermi != 0:
        ax.set_title(f"DOS and PDOS (0 is Fermi energy)", fontsize=18)
    else:
        ax.set_title(f"DOS and PDOS (0 is 0)", fontsize=18)
    ax.vlines(
        0,
        0,
        1,
        transform=ax.get_xaxis_transform(),
        color="grey",
        linewidths=0.5,
        linestyles="dashed",
    )

    if case == 0:
        ax.fill_between(dos[0], 0, dos[1], lw=0, color="grey", alpha=0.3, label="DOS")
        ax.plot(dos[0], dos[2], "-", lw=1, color="black", alpha=0.7, label="PDOS")
        ncol = 1
    if case == 1:
        ax.fill_between(
            dos[0], 0, dos[1], lw=0, color="blue", alpha=0.3, label="DOS (up)"
        )
        ax.plot(dos[0], dos[3], "-", lw=1, color="blue", alpha=0.7, label="PDOS (up)")
        ax.fill_between(
            dos[0], 0, -dos[2], lw=0, color="red", alpha=0.3, label="DOS (down)"
        )
        ax.plot(dos[0], -dos[4], "-", lw=1, color="red", alpha=0.7, label="PDOS (down)")
        ncol = 2
    if case == 2:
        ax.fill_between(dos[0], 0, dos[1], lw=0, color="grey", alpha=0.3, label="DOS")
        ax.fill_between(dos[0], -dos[1], 0, lw=0, color="grey", alpha=0.3, label="-DOS")
        ax.plot(dos[0], dos[2], "-", lw=1, color="blue", alpha=0.7, label="PDOS (up)")
        ax.plot(dos[0], -dos[3], "-", lw=1, color="red", alpha=0.7, label="PDOS (down)")
        ncol = 1

    ax.set_xlim(*tuple(window))
    ax.legend(loc="best", ncol=ncol)

    png_path = f"{output_name}.png"
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    print(f"Total dos plot  is in {abspath(png_path)}")
    plt.close()


def plot_projected(
    l, dos, output_name, title, case, window=None, efermi=0.0, labels=None, factor=None
):
    if l in ["s", "p", "d"]:
        labels = {
            "s": ["s"],
            "p": ["$p_z$", "$p_y$", "$p_x$"],
            "d": ["$d_{z^2}$", "$d_{zx}$", "$d_{zy}$", "$d_{x^2 - y^2}$", "$d_{xy}$"],
        }
        labels = labels[l]

        factor = {"s": 1, "p": 3, "d": 5}
        factor = factor[l]
    elif labels is None or factor is None:
        raise ValueError
    dos = filter_window(dos, window=window, efermi=efermi)
    if window is None:
        window = (np.amin(dos[0]), np.amax(dos[0]))
    fig, axs = plt.subplots(factor, 1, figsize=(9, factor * 2))
    if len(labels) == 1:
        axs = [axs]
    fig.subplots_adjust(hspace=0)
    for i in range(0, factor):
        ax = axs[i]
        ax.set_ylabel("DOS, states/eV", fontsize=12)
        if i == factor - 1:
            ax.set_xlabel("E, ev", fontsize=15)
        else:
            ax.axes.get_xaxis().set_visible(False)
        ax.set_xlim(*tuple(window))
        ax.vlines(
            0,
            0,
            1,
            transform=ax.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )
        if title is not None and i == 0:
            ax.set_title(title)

        if case == 0:
            ax.fill_between(
                dos[0], 0, dos[1], lw=0, color="black", alpha=0.3, label="p"
            )
            ax.plot(
                dos[0],
                dos[2 + i],
                "-",
                lw=0.5,
                color="black",
                alpha=0.8,
                label=labels[i],
            )
        elif case in [1, 2]:
            ax.fill_between(
                dos[0], 0, dos[1], lw=0, color="blue", alpha=0.3, label=f"{l} (up)"
            )
            ax.plot(
                dos[0],
                dos[3 + 2 * i],
                "-",
                lw=0.5,
                color="blue",
                alpha=0.8,
                label=f"{labels[i]} (up)",
            )
            ax.fill_between(
                dos[0], 0, -dos[2], lw=0, color="red", alpha=0.3, label=f"{l} (down)"
            )
            ax.plot(
                dos[0],
                -dos[4 + 2 * i],
                "-",
                lw=0.5,
                color="red",
                alpha=0.8,
                label=f"{labels[i]} (down)",
            )
        ax.legend(loc=(1.025, 0.2), bbox_transform=ax.transAxes)

    plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")
    plt.close()


def plot_projected_relative(
    l,
    dos,
    output_name,
    title,
    case,
    window=None,
    efermi=0.0,
    normalize=False,
    labels=None,
    factor=None,
):
    if l in ["s", "p", "d"]:
        colors = {
            "s": ["#9737FF"],
            "p": ["#0000FF", "#00FF00", "#FF0000"],
            "d": ["#00FF00", "#FF00FF", "#00FFFF", "#3E3847", "#FFD600"],
        }
        labels = {
            "s": ["s"],
            "p": ["$p_z$", "$p_y$", "$p_x$"],
            "d": ["$d_{z^2}$", "$d_{zx}$", "$d_{zy}$", "$d_{x^2 - y^2}$", "$d_{xy}$"],
        }
        colors = colors[l]
        labels = labels[l]

        factor = {"s": 1, "p": 3, "d": 5}
        factor = factor[l]
    elif labels is None or factor is None:
        raise ValueError
    else:
        colors = [
            "#0000FF",
            "#FF0000",
            "#00FF00",
            "#FF00FF",
            "#00FFFF",
            "#3E3847",
            "#FFD600",
            "#366B35",
            "#FF6F00",
        ]
    dos = filter_window(dos, window=window, efermi=efermi)
    if window is None:
        window = (np.amin(dos[0]), np.amax(dos[0]))
    if case == 0:
        fig, axs = plt.subplots(figsize=(9, 4))
        ax1 = axs
    elif case in [1, 2]:
        fig, axs = plt.subplots(2, 1, figsize=(9, 4))
        ax1 = axs[0]
        ax2 = axs[1]
    fig.subplots_adjust(hspace=0)
    if normalize:
        ax1.set_ylabel("PDOS / LDOS", fontsize=12)
        if case in [1, 2]:
            ax2.set_ylabel("PDOS / LDOS", fontsize=12)
    else:
        ax1.set_ylabel("DOS, states/eV", fontsize=12)
        if case in [1, 2]:
            ax2.set_ylabel("DOS, states/eV", fontsize=12)

    ax1.set_xlim(*tuple(window))
    ax1.vlines(
        0,
        0,
        1,
        transform=ax1.get_xaxis_transform(),
        color="grey",
        linewidths=0.5,
        linestyles="dashed",
    )

    if case in [1, 2]:
        ax2.set_xlabel("E, ev", fontsize=15)
        ax2.set_xlim(*tuple(window))
        ax2.vlines(
            0,
            0,
            1,
            transform=ax2.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )
        ax1.axes.get_xaxis().set_visible(False)
    else:
        ax1.set_xlabel("E, ev", fontsize=15)

    if title is not None:
        ax1.set_title(title)

    if case in [1, 2]:
        y_up = np.zeros(dos[0].shape)
        y_prev_up = np.zeros(dos[0].shape)
        y_down = np.zeros(dos[0].shape)
        y_prev_down = np.zeros(dos[0].shape)
    else:
        y = np.zeros(dos[0].shape)
        y_prev = np.zeros(dos[0].shape)

    for i in range(0, factor):
        if case == 0:
            y += dos[2 + i]
            if normalize:
                ax1.fill_between(
                    dos[0],
                    y_prev / dos[1],
                    y / dos[1],
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
            else:
                ax1.fill_between(
                    dos[0],
                    y_prev,
                    y,
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
            y_prev = deepcopy(y)
        elif case in [1, 2]:
            y_up += dos[3 + 2 * i]
            y_down += dos[4 + 2 * i]
            if normalize:
                ax1.fill_between(
                    dos[0],
                    y_prev_up / dos[1],
                    y_up / dos[1],
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
                ax2.fill_between(
                    dos[0],
                    y_prev_down / dos[2],
                    y_down / dos[2],
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
            else:
                ax1.fill_between(
                    dos[0],
                    y_prev_up,
                    y_up,
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
                ax2.fill_between(
                    dos[0],
                    y_prev_down,
                    y_down,
                    lw=0,
                    color=colors[i % factor],
                    label=f"{labels[i]}",
                )
            y_prev_up = deepcopy(y_up)
            y_prev_down = deepcopy(y_down)

    if normalize:
        ax1.set_ylim(0, 1)
        ax2.set_ylim(1, 0)
    else:
        ax1.set_ylim(0, None)
        ax2.set_ylim(0, None)
        ax2.invert_yaxis()
    if matplotlib_version >= "3.7.0":
        ax1.legend(
            loc=(1.025, 0.2), bbox_transform=ax1.transAxes, title="up", reverse=True
        )
    else:
        ax1.legend(loc=(1.025, 0.2), bbox_transform=ax1.transAxes, title="up")
    ax2.legend(loc=(1.025, 0.2), bbox_transform=ax2.transAxes, title="down")

    plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")
    plt.close()


def prepare_directories(output_path, seedname, separate, suffix):
    try:
        makedirs(join(output_path, seedname))
    except FileExistsError:
        pass
    try:
        makedirs(join(output_path, seedname, "summed-by-atom"))
    except FileExistsError:
        pass
    try:
        makedirs(join(output_path, seedname, "summed-by-atom", "wfc-in-atoms"))
    except FileExistsError:
        pass
    if separate:
        try:
            makedirs(join(output_path, seedname, "individual-plots"))
        except FileExistsError:
            pass
        try:
            makedirs(
                join(
                    output_path,
                    f"{seedname}{suffix}",
                    "individual-plots",
                    "wfc-in-atoms",
                )
            )
        except FileExistsError:
            pass


def manager(
    input_path,
    filpdos=None,
    output_path=".",
    window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    verbose=False,
):
    suffix = ""
    if relative or normalize:
        suffix += "-"
    if relative:
        suffix += "r"
    if normalize:
        suffix += "n"

    try:
        makedirs(output_path)
    except FileExistsError:
        pass

    filpdos, files = analyse_input_folder(input_path, filpdos)

    for s_i, seedname in enumerate(filpdos):
        # Preparations
        print(
            f"{YELLOW}({s_i + 1}/{len(filpdos)}) "
            + f"Start to work with {seedname} seedname{RESET}"
        )

        case = detect_case(join(input_path, f"{seedname}.pdos_tot"))
        print(f"{CASES[case]} case detected.")

        prepare_directories(output_path, f"{seedname}{suffix}", separate, suffix)

        atoms, wfcs = decompose_filenames(
            [file for file in files[s_i] if ".pdos_tot" != file[-9:]]
        )

        # Plot PDOS vs DOS
        total_dos = load_dos(join(input_path, f"{seedname}.pdos_tot"))
        plot_pdos_tot(
            total_dos,
            join(output_path, f"{seedname}{suffix}", "pdos-vs-dos"),
            case,
            window=window,
            efermi=efermi,
        )

        # Summ PDOS for the same atom types
        combine_by_atoms(input_path, seedname, output_path, atoms, wfcs, case, suffix)

        # Plot summed PDOS
        for atom in atoms:
            for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                path = join(
                    input_path,
                    f"{seedname}{suffix}",
                    "summed-by-atom",
                    f"{atom}_wfc#{wfc_number}({wfc_symbol})",
                )
                dos = load_dos(path)
                if efermi == 0:
                    title = f"PDOS for {atom} ({wfc_symbol}) (0 is 0)"
                else:
                    title = f"PDOS for {atom} ({wfc_symbol}) (0 is Fermi energy)"
                if relative:
                    plot_projected_relative(
                        wfc_symbol,
                        dos,
                        path,
                        title,
                        case,
                        window=window,
                        efermi=efermi,
                        normalize=normalize,
                    )
                else:
                    plot_projected(
                        wfc_symbol, dos, path, title, case, window=window, efermi=efermi
                    )

        # Plot all individual plots
        if separate:
            for atom in atoms:
                for a_number in atoms[atom]:
                    for wfc_symbol, wfc_number in wfcs[(atom, a_number)]:
                        dos = load_dos(
                            join(
                                input_path,
                                f"{seedname}.pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})",
                            )
                        )
                        if efermi == 0:
                            title = (
                                f"PDOS for {atom} #{a_number} ({wfc_symbol}) (0 is 0)"
                            )
                        else:
                            title = f"PDOS for {atom} #{a_number} ({wfc_symbol}) (0 is Fermi energy)"
                        if relative:
                            plot_projected_relative(
                                wfc_symbol,
                                dos,
                                join(
                                    output_path,
                                    f"{seedname}{suffix}",
                                    "individual-plots",
                                    f"pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})",
                                ),
                                title,
                                case,
                                window=window,
                                efermi=efermi,
                                normalize=normalize,
                            )
                        else:
                            plot_projected(
                                wfc_symbol,
                                dos,
                                join(
                                    output_path,
                                    f"{seedname}{suffix}",
                                    "individual-plots",
                                    f"pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})",
                                ),
                                title,
                                case,
                                window=window,
                                efermi=efermi,
                            )

        # Plot orbital contribution in the whole atom
        for atom in atoms:
            custom_dos = [load_dos(join(input_path, f"{seedname}.pdos_tot"))[0]]
            wfc_symbols = []
            for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                wfc_symbols.append(f"#{wfc_number}, {wfc_symbol}")
                wfs_dos = load_dos(
                    join(
                        input_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        f"{atom}_wfc#{wfc_number}({wfc_symbol})",
                    )
                )
                if case == 0:
                    if len(custom_dos) == 1:
                        custom_dos.append(deepcopy(wfs_dos[1]))
                    else:
                        custom_dos[1] += wfs_dos[1]
                    custom_dos.append(wfs_dos[1])
                elif case in [1, 2]:
                    if len(custom_dos) == 1:
                        custom_dos.append(deepcopy(wfs_dos[1]))
                        custom_dos.append(deepcopy(wfs_dos[2]))
                    else:
                        custom_dos[1] += wfs_dos[1]
                        custom_dos[2] += wfs_dos[2]
                    custom_dos.append(wfs_dos[1])
                    custom_dos.append(wfs_dos[2])
            if efermi == 0:
                title = f"wfc contribution to {atom} (0 is 0)"
            else:
                title = f"wfc contribution to {atom} (0 is Fermi energy)"

            if relative:
                plot_projected_relative(
                    atom,
                    np.array(custom_dos),
                    join(
                        output_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        "wfc-in-atoms",
                        atom,
                    ),
                    title,
                    case,
                    window=window,
                    efermi=efermi,
                    normalize=normalize,
                    labels=wfc_symbols,
                    factor=len(wfc_symbols),
                )
            else:
                plot_projected(
                    atom,
                    np.array(custom_dos),
                    join(
                        output_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        "wfc-in-atoms",
                        atom,
                    ),
                    title,
                    case,
                    window=window,
                    efermi=efermi,
                    labels=wfc_symbols,
                    factor=len(wfc_symbols),
                )

        # Plot orbital contribution in the whole atom (individual atoms)
        if separate:
            for atom in atoms:
                for a_number in atoms[atom]:
                    custom_dos = [load_dos(join(input_path, f"{seedname}.pdos_tot"))[0]]
                    wfc_symbols = []
                    for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                        wfc_symbols.append(f"#{wfc_number}, {wfc_symbol}")
                        wfs_dos = load_dos(
                            join(
                                input_path,
                                f"{seedname}.pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})",
                            )
                        )
                        if case == 0:
                            if len(custom_dos) == 1:
                                custom_dos.append(deepcopy(wfs_dos[1]))
                            else:
                                custom_dos[1] += wfs_dos[1]
                            custom_dos.append(wfs_dos[1])
                        elif case in [1, 2]:
                            if len(custom_dos) == 1:
                                custom_dos.append(deepcopy(wfs_dos[1]))
                                custom_dos.append(deepcopy(wfs_dos[2]))
                            else:
                                custom_dos[1] += wfs_dos[1]
                                custom_dos[2] += wfs_dos[2]
                            custom_dos.append(wfs_dos[1])
                            custom_dos.append(wfs_dos[2])
                    if efermi == 0:
                        title = f"wfc contribution to #{a_number} {atom} (0 is 0)"
                    else:
                        title = f"wfc contribution to #{a_number} {atom} (0 is Fermi energy)"

                    if relative:
                        plot_projected_relative(
                            f"#{a_number} {atom}",
                            np.array(custom_dos),
                            join(
                                output_path,
                                f"{seedname}{suffix}",
                                "individual-plots",
                                "wfc-in-atoms",
                                f"#{a_number}_{atom}",
                            ),
                            title,
                            case,
                            window=window,
                            efermi=efermi,
                            normalize=normalize,
                            labels=wfc_symbols,
                            factor=len(wfc_symbols),
                        )
                    else:
                        plot_projected(
                            f"#{a_number} {atom}",
                            np.array(custom_dos),
                            join(
                                output_path,
                                f"{seedname}{suffix}",
                                "individual-plots",
                                "wfc-in-atoms",
                                f"#{a_number} {atom}",
                            ),
                            title,
                            case,
                            window=window,
                            efermi=efermi,
                            labels=wfc_symbols,
                            factor=len(wfc_symbols),
                        )

        # Plot atom contribution to total
        custom_dos = load_dos(join(input_path, f"{seedname}.pdos_tot"))
        if case == 0:
            custom_dos = [custom_dos[0], custom_dos[2]]
        elif case == 1:
            custom_dos = [custom_dos[0], custom_dos[3], custom_dos[4]]
        elif case == 2:
            custom_dos = [custom_dos[0], custom_dos[2], custom_dos[3]]
        labels = []
        for atom in atoms:
            labels.append(atom)
            if case == 0:
                custom_dos.append(None)
            elif case in [1, 2]:
                custom_dos.append(None)
                custom_dos.append(None)
            for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                wfs_dos = load_dos(
                    join(
                        input_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        f"{atom}_wfc#{wfc_number}({wfc_symbol})",
                    )
                )
                if case == 0:
                    if custom_dos[-1] is None:
                        custom_dos[-1] = deepcopy(wfs_dos[1])
                    else:
                        custom_dos[-1] += wfs_dos[1]
                elif case in [1, 2]:
                    if custom_dos[-1] is None:
                        custom_dos[-2] = deepcopy(wfs_dos[1])
                        custom_dos[-1] = deepcopy(wfs_dos[2])
                    else:
                        custom_dos[-2] += wfs_dos[1]
                        custom_dos[-1] += wfs_dos[2]
            if efermi == 0:
                title = f"Atom contribution to total PDOS (0 is 0)"
            else:
                title = f"Atom contribution to total PDOS (0 is Fermi energy)"

            if relative:
                plot_projected_relative(
                    "Total PDOS",
                    np.array(custom_dos),
                    join(
                        output_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        "Atom-contributions",
                    ),
                    title,
                    case,
                    window=window,
                    efermi=efermi,
                    normalize=normalize,
                    labels=labels,
                    factor=len(labels),
                )
            else:
                plot_projected(
                    "Total PDOS",
                    np.array(custom_dos),
                    join(
                        output_path,
                        f"{seedname}{suffix}",
                        "summed-by-atom",
                        "Atom-contributions",
                    ),
                    title,
                    case,
                    window=window,
                    efermi=efermi,
                    labels=labels,
                    factor=len(labels),
                )

        # Plot atom contribution to total (individual atoms)
        if separate:
            custom_dos = load_dos(join(input_path, f"{seedname}.pdos_tot"))
            if case == 0:
                custom_dos = [custom_dos[0], custom_dos[2]]
            elif case == 1:
                custom_dos = [custom_dos[0], custom_dos[3], custom_dos[4]]
            elif case == 2:
                custom_dos = [custom_dos[0], custom_dos[2], custom_dos[3]]
            labels = []
            for atom in atoms:
                for a_number in atoms[atom]:
                    labels.append(f"#{a_number} {atom}")
                    if case == 0:
                        custom_dos.append(None)
                    elif case in [1, 2]:
                        custom_dos.append(None)
                        custom_dos.append(None)
                    for wfc_symbol, wfc_number in wfcs[(atom, atoms[atom][0])]:
                        wfs_dos = load_dos(
                            join(
                                input_path,
                                f"{seedname}.pdos_atm#{a_number}({atom})_wfc#{wfc_number}({wfc_symbol})",
                            )
                        )
                        if case == 0:
                            if custom_dos[-1] is None:
                                custom_dos[-1] = deepcopy(wfs_dos[1])
                            else:
                                custom_dos[-1] += wfs_dos[1]
                        elif case in [1, 2]:
                            if custom_dos[-1] is None:
                                custom_dos[-2] = deepcopy(wfs_dos[1])
                                custom_dos[-1] = deepcopy(wfs_dos[2])
                            else:
                                custom_dos[-2] += wfs_dos[1]
                                custom_dos[-1] += wfs_dos[2]
                    if efermi == 0:
                        title = f"Atom contribution to total PDOS (0 is 0)"
                    else:
                        title = f"Atom contribution to total PDOS (0 is Fermi energy)"

                    if relative:
                        plot_projected_relative(
                            "Total PDOS",
                            np.array(custom_dos),
                            join(
                                output_path,
                                f"{seedname}{suffix}",
                                "individual-plots",
                                "Atom-contributions",
                            ),
                            title,
                            case,
                            window=window,
                            efermi=efermi,
                            normalize=normalize,
                            labels=labels,
                            factor=len(labels),
                        )
                    else:
                        plot_projected(
                            "Total PDOS",
                            np.array(custom_dos),
                            join(
                                output_path,
                                f"{seedname}{suffix}",
                                "individual-plots",
                                "Atom-contributions",
                            ),
                            title,
                            case,
                            window=window,
                            efermi=efermi,
                            labels=labels,
                            factor=len(labels),
                        )

        print(
            f"{OK}"
            + f"Finish to work with {seedname} seedname, results are in "
            + f"  {abspath(join(output_path, f'{seedname}{suffix}'))}"
            + f"{RESET}"
        )


def create_parser():
    parser = ArgumentParser(
        description="Script for visualisation of density of states."
    )
    parser.add_argument(
        "-ip",
        "--input-path",
        metavar="path",
        type=str,
        required=True,
        help="Relative or absulute path to the folder with dos files.",
    )
    parser.add_argument(
        "-f",
        "--filpdos",
        metavar="filpdos",
        type=str,
        default=None,
        help="Prefix for output files containing PDOS(E). "
        + "As specified in the QE projwfc.x input file.",
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
        "-w",
        "--window",
        metavar=("min", "max"),
        type=float,
        nargs=2,
        default=None,
        help="Energy window for the plots. "
        + "By default whole range present in the files is plotted.",
    )
    parser.add_argument(
        "-ef",
        "--efermi",
        metavar="energy",
        type=float,
        default=0.0,
        help="Fermi energy, zero by default.",
    )
    parser.add_argument(
        "-s",
        "--separate",
        action="store_true",
        default=False,
        help="Whenever to plot projected DOS for each atom of the same type separately.",
    )
    parser.add_argument(
        "-r",
        "--relative",
        action="store_true",
        default=False,
        help="Whenever to use relative style.",
    )
    parser.add_argument(
        "-n",
        "--normalize",
        action="store_true",
        default=False,
        help="Whenever to use normalize relative style.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )

    return parser
