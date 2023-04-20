#! /usr/local/bin/python3

import re
from argparse import ArgumentParser
from copy import deepcopy
from os import makedirs, walk
from os.path import abspath, join

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import __version__ as matplotlib_version

from rad_tools.dos import DOS
from rad_tools.routines import cprint


def detect_seednames(input_path):
    r"""
    Analyze input folder, detects seednames for the dos output files.

    Parameters
    ----------
    input_path : str
        Directory with DOS files.

    Returns
    -------
    seednames : list
        List of seednames found in ``input_path``.
    """

    pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf_0-9j.]*\\)"

    # Get list of files in the folder
    files = []
    for dirpath, dirnames, filenames in walk(input_path):
        files.extend(filenames)
        break

    seednames = set()
    for file in files:
        if ".pdos_tot" in file and ".pdos_tot" == file[-9:]:
            seednames.add(file[:-9])
        elif re.match(f".*{pattern}$", file):
            seednames.add(re.split(f"{pattern}$", file)[0])
    seednames = list(seednames)

    return seednames


def manager(
    input_path,
    seedname=None,
    output_path=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    verbose=False,
    interactive=False,
):
    makedirs(output_path, exist_ok=True)
    suffix = ""
    if relative or normalize:
        suffix += "-"
    if relative:
        suffix += "r"
    if normalize:
        suffix += "n"

    if seedname is None:
        seednames = detect_seednames(input_path)
        print(f"Following DOS seednames (filpdos) are detected:")
        for item in seednames:
            cprint(f"   * {item}", colour="yellow")
    else:
        seednames = [seedname]

    for s_i, seedname in enumerate(seednames):
        output_root = join(output_path, f"{seedname}{suffix}")
        # Preparations
        cprint(
            f"({s_i + 1}/{len(seednames)}) Start to work with {seedname} seedname",
            colour="yellow",
        )
        makedirs(output_root, exist_ok=True)

        dos = DOS(seedname, input_path, efermi=efermi, window=energy_window)
        print(f"{dos.case_name} case detected.")

        # Plot PDOS vs DOS
        dos.plot_pdos_tot(
            output_name=join(output_root, "pdos-vs-dos"),
            interactive=interactive,
            ylim=dos_window,
        )

        # Plot PDOS for each atom/wfc (summed for atom numbers)
        makedirs(join(output_root, "summed"), exist_ok=True)
        data = {}
        for atom, atom_number, wfc, wfc_number in dos:
            if atom not in data:
                data[atom] = []
            data[atom].append((wfc, wfc_number))
        for atom in data:
            for wfcs in data[atom]:
                if efermi == 0:
                    title = f"PDOS for {atom} ({wfcs}) (0 is 0)"
                else:
                    title = f"PDOS for {atom} ({wfcs}) (0 is Fermi energy)"
                pdos = dos.pdos(
                    atom=atom,
                    atom_numbers=None,
                    wfc=wfcs[0],
                    wfc_numbers=wfcs[1],
                    squeeze=True,
                )
                pdos.plot_projected(
                    output_name=join(
                        output_root, "summed", f"{atom}_{wfcs[0]}_{wfcs[1]}"
                    ),
                    ylim=None,
                    title=title,
                    relative=relative,
                    normalize=normalize,
                    interactive=interactive,
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
        "-s",
        "--seedname",
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
        "-ew",
        "--energy-window",
        metavar=("min", "max"),
        type=float,
        nargs=2,
        default=None,
        help="Energy window for the plots. "
        + "By default whole range present in the files is plotted.",
    )
    parser.add_argument(
        "-dw",
        "--dos-window",
        metavar=("min", "max"),
        type=float,
        nargs=2,
        default=None,
        help="DOS window for the plots. "
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
        "-sep",
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
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        default=False,
        help="Interactive plotting.",
    )

    return parser
