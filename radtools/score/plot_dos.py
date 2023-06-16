# TODO Fix background total + non-colinear-nonspin-orbit

import re
from argparse import ArgumentParser
from os import makedirs, walk
from os.path import abspath, isfile, join

from termcolor import cprint
from tqdm import tqdm

from radtools.dos.dos import DOSQE, PATTERN
from radtools.dos.pdos import COLOURS, PDOS, plot_projected


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

    # Get list of files in the folder
    files = []
    for dirpath, dirnames, filenames in walk(input_path):
        files.extend(filenames)
        break

    seednames = set()
    for file in files:
        if ".pdos_tot" in file and ".pdos_tot" == file[-9:]:
            seednames.add(file[:-9])
        elif re.match(f".*{PATTERN}$", file):
            seednames.add(re.split(f"{PATTERN}$", file)[0])
    seednames = list(seednames)

    return seednames


def plot_orbital_resolved(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Orbital-resolved PDOS:", "green")
    local_output = join(output_root, "orbital-resolved")
    makedirs(local_output, exist_ok=True)
    data = {}
    for atom, atom_number, wfc, wfc_number in dos:
        if atom not in data:
            data[atom] = []
        data[atom].append((wfc, wfc_number))

    # Avoid repetitions
    for atom in data:
        data[atom] = list(set(data[atom]))

    for atom in data:
        for wfc, wfc_number in data[atom]:
            if separate:
                atom_numbers = dos.atom_numbers(atom)
            else:
                atom_numbers = [None]

            for atom_number in tqdm(atom_numbers, desc=f"  {atom} {wfc} #{wfc_number}"):
                if separate:
                    atom_name = f"{atom}#{atom_number}"
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
                    background_total=background_total,
                )
                if background_total:
                    pdos.projectors_group = "Total PDOS"
                if save_txt:
                    pdos.dump_txt(
                        join(local_output, f"{atom_name}_{wfc}#{wfc_number}.txt")
                    )
                plot_projected(
                    pdos=pdos,
                    efermi=efermi,
                    output_name=join(local_output, f"{atom_name}_{wfc}#{wfc_number}"),
                    title=title,
                    xlim=energy_window,
                    ylim=dos_window,
                    relative=relative,
                    normalize=normalize,
                    interactive=interactive,
                    save_pickle=save_pickle,
                    colours=colours,
                    legend_fontsize=legend_fontsize,
                    axes_labels_fontsize=axes_labels_fontsize,
                    title_fontsize=title_fontsize,
                )
    cprint(f"Results are in {abspath(local_output)}", "blue")


def plot_atom_resolved(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Orbital's contribution for each atom.", "green")
    local_output = join(output_root, "atom-resolved")
    makedirs(local_output, exist_ok=True)

    for atom in dos.atoms:
        if separate:
            atom_numbers = dos.atom_numbers(atom)
        else:
            atom_numbers = [None]
        for atom_number in tqdm(atom_numbers, desc=f"  {atom}"):
            if separate:
                atom_name = f"{atom}#{atom_number}"
            else:
                atom_name = atom
            projectors = []
            pdos = []
            for wfc, wfc_number in dos.wfcs(atom, atom_number):
                projectors.append(f"{wfc} #{wfc_number}")
                pdos.append(dos.pdos(atom, wfc, wfc_number, atom_number).ldos)

            if background_total:
                pdos = PDOS(
                    energy=dos.energy,
                    pdos=pdos,
                    ldos=dos.total_pdos(),
                    projectors_group=atom_name,
                    projectors=projectors,
                    spin_pol=dos.case in [2, 3],
                )
                pdos.projectors_group = "Total PDOS"
            else:
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
                colours=colours,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )
    cprint(f"Results are in {abspath(local_output)}", "blue")


def plot_atom_to_total(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Atom's contributions into total PDOS:", "green")
    projectors = []
    pdos = []
    for atom in dos.atoms:
        if separate:
            atom_numbers = dos.atom_numbers(atom)
        else:
            atom_numbers = [None]
        for atom_number in atom_numbers:
            if separate:
                atom_name = f"{atom}#{atom_number}"
            else:
                atom_name = atom
            projectors.append(atom_name)
            for i, (wfc, wfc_number) in enumerate(dos.wfcs(atom, atom_number)):
                if i == 0:
                    ldos = dos.pdos(atom, wfc, wfc_number, atom_number).ldos
                else:
                    ldos += dos.pdos(atom, wfc, wfc_number, atom_number).ldos
            pdos.append(ldos)
    if background_total:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            ldos=dos.total_pdos(),
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.case in [2, 3],
        )
        pdos.projectors_group = "Total PDOS"
    else:
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
        colours=colours,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
    )
    cprint(f"Result is in {abspath(output_root)}", "blue")


def plot_custom(
    dos,
    custom,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Plotting custom plot", "green")
    print("Input is understood as:")
    projectors = []
    pdos = []
    if (
        labels is not None
        and len(labels) != len(custom)
        and len(labels) != len(custom) + 1
    ):
        raise ValueError(
            f"Got {len(labels)} labels, but {len(custom)} PDOS, have to be the same or n custom, n+1 labels."
        )
    if labels is not None and len(labels) == len(custom) + 1:
        if labels[0].lower() == "none":
            total_label = None
        elif labels[0].lower() == "default":
            total_label = "default"
        else:
            total_label = labels[0]
        labels = labels[1:]
    else:
        total_label = "default"
    for i_e, entry in enumerate(custom):
        if labels is not None:
            projectors.append(labels[i_e])
        else:
            projectors.append(entry)
        cprint(f'"{entry}":', "green")
        entry = entry.replace(" ", "").replace(")", "")
        subentries = entry.split(";")
        pdos_element = None
        for subentry in subentries:
            atom_part = subentry.split("(")[0]
            atom = atom_part.split("#")[0]
            if "#" in subentry:
                atom_numbers = list(map(int, atom_part.split("#")[1:]))
            else:
                atom_numbers = dos.atom_numbers(atom)

            cprint(
                f"  * PDOS is summed among the following {atom} atoms:",
                "green",
                end="\n      ",
            )
            for i, number in enumerate(atom_numbers):
                print(f"{atom}#{number}", end="")
                if i != len(atom_numbers) - 1:
                    print(end=", ")
                else:
                    print()

            if "(" in subentry:
                wfc_parts = subentry.split("(")[1].split(",")
                wfcs = []
                for wfc_part in wfc_parts:
                    wfc = wfc_part.split("#")[0]
                    if "#" in wfc_part:
                        wfc_numbers = list(map(int, wfc_part.split("#")[1:]))
                        for number in wfc_numbers:
                            wfcs.append((wfc, number))
                    else:
                        wfc_list = dos.wfcs(atom)
                        for name, number in wfc_list:
                            if name == wfc:
                                wfcs.append((wfc, number))
            else:
                wfcs = dos.wfcs(atom)

            print(
                f"  * For each {atom} atom PDOS is summed among the following projections:",
                end="\n      ",
            )
            for i, (name, number) in enumerate(wfcs):
                print(f"{name}#{number}", end="")
                if i != len(wfcs) - 1:
                    print(end=", ")
                else:
                    print()

            for name, number in wfcs:
                if pdos_element is None:
                    pdos_element = dos.pdos(
                        atom=atom,
                        wfc=name,
                        wfc_number=number,
                        atom_numbers=atom_numbers,
                    ).ldos
                else:
                    pdos_element += dos.pdos(
                        atom=atom,
                        wfc=name,
                        wfc_number=number,
                        atom_numbers=atom_numbers,
                    ).ldos

        pdos.append(pdos_element)

    if background_total:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            ldos=dos.total_pdos(),
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.case in [2, 3],
        )
    else:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            projectors_group="Total (sum)",
            projectors=projectors,
            spin_pol=dos.case in [2, 3],
        )

    if isfile(join(output_root, "custom.png")):
        i = 1
        while isfile(join(output_root, f"custom{i}.png")):
            i += 1
        output_name = f"custom{i}"
    else:
        output_name = "custom"

    if save_txt:
        pdos.dump_txt(join(output_root, f"{output_name}.txt"))
    plot_projected(
        pdos=pdos,
        efermi=efermi,
        output_name=join(output_root, output_name),
        xlim=energy_window,
        ylim=dos_window,
        relative=relative,
        normalize=normalize,
        interactive=interactive,
        save_pickle=save_pickle,
        colours=colours,
        total_label=total_label,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
    )
    cprint(f"Result is in {abspath(join(output_root, f'{output_name}.png'))}", "blue")


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
    background_total=False,
    custom=None,
    colours=None,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    r"""
    :ref:`rad-plot-dos` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-plot-dos>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    # Create the output directory if it does not exist
    makedirs(output_path, exist_ok=True)

    if colours is None:
        colours = COLOURS

    suffix = ""
    if relative:
        suffix += "-relative"
    if normalize:
        suffix += "-normalized"
    if separate:
        suffix += "-separate"
    if background_total:
        suffix += "-vstotal"

    # Detect seednames if not provided.
    if seedname is None:
        seednames = detect_seednames(input_path)
        print(f"Following DOS seednames (filpdos) are detected:")
        for item in seednames:
            cprint(f"   * {item}", "green", attrs=["bold"])
    else:
        seednames = [seedname]

    # Work with each seedname.
    for s_i, seedname in enumerate(seednames):
        cprint(
            f"({s_i + 1}/{len(seednames)}) Start to work with {seedname} seedname",
            "yellow",
            attrs=["bold"],
        )
        # Preparations
        output_root = join(output_path, f"{seedname}{suffix}")
        makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_path, energy_window=energy_window, efermi=efermi)
        print(f"{dos.casename} case detected.")
        for atom in dos.atoms:
            print(f"  {len(dos.atom_numbers(atom))} {atom} detected")

        # Plot predefined plots
        if custom is None:
            # Plot PDOS vs DOS
            cprint("Total DOS vs total PDOS", "green")
            dos.plot_pdos_tot(
                output_name=join(output_root, "pdos-vs-dos"),
                interactive=interactive,
                efermi=efermi,
                xlim=energy_window,
                ylim=dos_window,
                save_pickle=save_pickle,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )
            cprint(f"Result is in {join(output_root, 'pdos-vs-dos')}", "blue")

            # Plot PDOS for each atom/wfc
            plot_orbital_resolved(
                dos=dos,
                output_root=output_root,
                energy_window=energy_window,
                dos_window=dos_window,
                efermi=efermi,
                separate=separate,
                relative=relative,
                normalize=normalize,
                interactive=interactive,
                save_pickle=save_pickle,
                save_txt=save_txt,
                background_total=background_total,
                colours=colours,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )

            # Plot wfc contribution into each atom
            plot_atom_resolved(
                dos=dos,
                output_root=output_root,
                energy_window=energy_window,
                dos_window=dos_window,
                efermi=efermi,
                separate=separate,
                relative=relative,
                normalize=normalize,
                interactive=interactive,
                save_pickle=save_pickle,
                save_txt=save_txt,
                background_total=background_total,
                colours=colours,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )

            # Plot atom contributions into total PDOS
            plot_atom_to_total(
                dos=dos,
                output_root=output_root,
                energy_window=energy_window,
                dos_window=dos_window,
                efermi=efermi,
                separate=separate,
                relative=relative,
                normalize=normalize,
                interactive=interactive,
                save_pickle=save_pickle,
                save_txt=save_txt,
                background_total=background_total,
                colours=colours,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )
        # Plot custom plot
        else:
            plot_custom(
                dos=dos,
                custom=custom,
                output_root=output_root,
                energy_window=energy_window,
                dos_window=dos_window,
                efermi=efermi,
                relative=relative,
                normalize=normalize,
                interactive=interactive,
                save_pickle=save_pickle,
                save_txt=save_txt,
                background_total=background_total,
                colours=colours,
                labels=labels,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
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
        default=".",
        help="Relative or absulute path to the folder with dos files.",
    )
    parser.add_argument(
        "-s",
        "--seedname",
        metavar="name",
        type=str,
        default=None,
        help="Prefix for input files with PDOS(E).",
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
        help="Whether to plot projected DOS for each atom of the same type separately.",
    )
    parser.add_argument(
        "-r",
        "--relative",
        action="store_true",
        default=False,
        help="Whether to use relative style.",
    )
    parser.add_argument(
        "-n",
        "--normalize",
        action="store_true",
        default=False,
        help="Whether to normalized PDOS values to 1.",
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
        help="Whether to save figures as .pickle files.",
    )
    parser.add_argument(
        "-st",
        "--save-txt",
        action="store_true",
        default=False,
        help="Whether to save some data as txt files.",
    )
    parser.add_argument(
        "-bt",
        "--background-total",
        action="store_true",
        default=False,
        help="Whether to use total PDOS as the background for all plots.",
    )
    parser.add_argument(
        "--custom",
        type=str,
        metavar="description",
        default=None,
        nargs="*",
        help="Custom PDOS plot. See docs for info.",
    )
    parser.add_argument(
        "-cls",
        "--colours",
        type=str,
        metavar="colours",
        default=None,
        nargs="*",
        help="Colours for the relative and custom plots.",
    )
    parser.add_argument(
        "-lbs",
        "--labels",
        type=str,
        metavar="labels",
        default=None,
        nargs="*",
        help="Labels for the custom plots.",
    )
    parser.add_argument(
        "-alfs",
        "--axes-labels-fontsize",
        type=int,
        default=14,
        metavar="fontsize",
        help="Fontsize of the labes of the axes.",
    )
    parser.add_argument(
        "-lfs",
        "--legend-fontsize",
        type=int,
        default=12,
        metavar="fontsize",
        help="Fontsize of the legend.",
    )
    parser.add_argument(
        "-tfs",
        "--title-fontsize",
        type=int,
        default=18,
        metavar="fontsize",
        help="Fontsize of the title.",
    )

    return parser
