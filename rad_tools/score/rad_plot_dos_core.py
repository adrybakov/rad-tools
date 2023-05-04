#! /usr/local/bin/python3

import re
from argparse import ArgumentParser
from os import makedirs, walk
from os.path import abspath, join

from tqdm import tqdm

from rad_tools.dos import DOSQE, plot_projected, PDOS
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
    save_pickle=False,
    save_txt=False,
):
    makedirs(output_path, exist_ok=True)

    suffix = ""
    if relative:
        suffix += "-relative"
    if normalize:
        suffix += "-normalized"
    if separate:
        suffix += "-separate"

    # Detect seednames if not provided.
    if seedname is None:
        seednames = detect_seednames(input_path)
        print(f"Following DOS seednames (filpdos) are detected:")
        for item in seednames:
            cprint(f"   * {item}", colour="yellow")
    else:
        seednames = [seedname]

    # Work with each seedname.
    for s_i, seedname in enumerate(seednames):
        cprint(
            f"({s_i + 1}/{len(seednames)}) Start to work with {seedname} seedname",
            colour="yellow",
        )
        # Preparations
        output_root = join(output_path, f"{seedname}{suffix}")
        makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_path, energy_window=energy_window)
        print(f"{dos.casename} case detected.")
        for atom in dos.atoms:
            print(f"    {len(dos.atom_numbers(atom))} of {atom} detected")

        # Plot PDOS vs DOS
        cprint("Total DOS vs total PDOS", colour="yellow")
        dos.plot_pdos_tot(
            output_name=join(output_root, "pdos-vs-dos"),
            interactive=interactive,
            efermi=efermi,
            xlim=energy_window,
            ylim=dos_window,
            save_pickle=save_pickle,
        )
        print(f"  Result is in {join(output_root, 'pdos-vs-dos')}")

        # Plot PDOS for each atom/wfc
        cprint("Orbital-resolved PDOS:", colour="yellow")
        local_output = join(output_root, "orbital-orbital-resolved")
        makedirs(local_output, exist_ok=True)
        data = {}
        for atom, atom_number, wfc, wfc_number in dos:
            if atom not in data:
                data[atom] = []
            data[atom].append((wfc, wfc_number))

        for atom in data:
            for wfc, wfc_number in data[atom]:
                if separate:
                    atom_numbers = dos.atom_numbers(atom)
                else:
                    atom_numbers = [None]

                for atom_number in tqdm(
                    atom_numbers, desc=f"  {atom} {wfc} #{wfc_number}"
                ):
                    if separate:
                        atom_name = f"{atom}{atom_number}"
                    else:
                        atom_name = atom

                    if efermi == 0:
                        title = f"PDOS for {atom_name} ({wfc} #{wfc_number}) (0 is 0)"
                    else:
                        title = f"PDOS for {atom_name} ({wfc} #{wfc_number}) (0 is Fermi energy)"

                    pdos = dos.pdos(
                        atom=atom,
                        wfc=wfc,
                        wfc_number=wfc_number,
                        atom_numbers=atom_number,
                    )
                    if save_txt:
                        pdos.dump_txt(
                            join(local_output, f"{atom_name}_{wfc}_{wfc_number}.txt")
                        )
                    plot_projected(
                        pdos=pdos,
                        efermi=efermi,
                        output_name=join(
                            local_output, f"{atom_name}_{wfc}_{wfc_number}"
                        ),
                        title=title,
                        xlim=energy_window,
                        ylim=dos_window,
                        relative=relative,
                        normalize=normalize,
                        interactive=interactive,
                        save_pickle=save_pickle,
                    )
        print(f"  Results are in {abspath(local_output)}")

        # Plot wfc contribution into each atom
        cprint("Orbital's contribution for each atom.", colour="yellow")
        local_output = join(output_root, "atom-orbital-resolved")
        makedirs(local_output, exist_ok=True)

        for atom in dos.atoms:
            if separate:
                atom_numbers = dos.atom_numbers(atom)
            else:
                atom_numbers = [None]
            for atom_number in tqdm(atom_numbers, desc=f"  {atom}"):
                if separate:
                    atom_name = f"{atom}{atom_number}"
                else:
                    atom_name = atom
                projectors = []
                pdos = []
                for wfc, wfc_number in dos.wfcs(atom, atom_number):
                    projectors.append(f"{wfc} #{wfc_number}")
                    pdos.append(dos.pdos(atom, wfc, wfc_number, atom_number).ldos)

                pdos = PDOS(
                    energy=dos.energy,
                    pdos=pdos,
                    projectors_group=atom_name,
                    projectors=projectors,
                    spin_pol=dos.case in [2, 3],
                )
                if efermi == 0:
                    title = f"PDOS for {atom_name} (0 is 0)"
                else:
                    title = f"PDOS for {atom_name} (0 is Fermi energy)"
                if save_txt:
                    pdos.dump_txt(join(local_output, f"{atom_name}.txt"))
                plot_projected(
                    pdos=pdos,
                    efermi=efermi,
                    output_name=join(local_output, atom_name),
                    title=title,
                    xlim=energy_window,
                    ylim=dos_window,
                    relative=relative,
                    normalize=normalize,
                    interactive=interactive,
                    save_pickle=save_pickle,
                )
        print(f"  Results are in {abspath(local_output)}")

        # Plot atom contributions into total PDOS
        cprint("Atom's contributions into total PDOS:", colour="yellow")
        projectors = []
        pdos = []
        for atom in dos.atoms:
            if separate:
                atom_numbers = dos.atom_numbers(atom)
            else:
                atom_numbers = [None]
            for atom_number in atom_numbers:
                if separate:
                    atom_name = f"{atom}{atom_number}"
                else:
                    atom_name = atom
                projectors.append(atom_name)
                for i, (wfc, wfc_number) in enumerate(dos.wfcs(atom, atom_number)):
                    if i == 0:
                        ldos = dos.pdos(atom, wfc, wfc_number, atom_number).ldos
                    else:
                        ldos += dos.pdos(atom, wfc, wfc_number, atom_number).ldos
                pdos.append(ldos)

        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.case in [2, 3],
        )
        if efermi == 0:
            title = f"Atom contribution in PDOS (0 is 0)"
        else:
            title = f"Atom contribution in PDOS (0 is Fermi energy)"
        if save_txt:
            pdos.dump_txt(join(output_root, "atomic-contributions.txt"))
        plot_projected(
            pdos=pdos,
            efermi=efermi,
            output_name=join(output_root, "atomic-contributions"),
            title=title,
            xlim=energy_window,
            ylim=dos_window,
            relative=relative,
            normalize=normalize,
            interactive=interactive,
            save_pickle=save_pickle,
        )
        print(f"  Result is in {abspath(local_output)}")


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
    parser.add_argument(
        "-sp",
        "--save-pickle",
        action="store_true",
        default=False,
        help="Whenever to save figures as .pickle files.",
    )
    parser.add_argument(
        "-st",
        "--save-txt",
        action="store_true",
        default=False,
        help="Whenever to save some data as txt files.",
    )

    return parser
