from argparse import ArgumentParser
from os import makedirs
from os.path import abspath, isfile, join

from termcolor import cprint

from radtools.dos.dos import DOSQE, detect_seednames
from radtools.dos.pdos import PDOS
from radtools.dos.plotting import plot_fatbands, COLOURS, plot_fatbands


def plot_custom(
    dos,
    custom,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    interactive=False,
    save_pickle=False,
    save_txt=False,
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

    pdos = PDOS(
        energy=dos.energy,
        pdos=pdos,
        ldos=dos.total_pdos(),
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

    plot_fatbands(
        pdos=pdos,
        efermi=efermi,
        output_name=join(output_root, output_name),
        ylim=energy_window,
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
    efermi=0.0,
    separate=False,
    verbose=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    custom=None,
    colours=None,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    r"""
    :ref:`rad-plot-fatbands` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-plot-fatbands>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    # Create the output directory if it does not exist
    makedirs(output_path, exist_ok=True)

    if colours is None:
        colours = COLOURS

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
        output_root = join(output_path, f"{seedname}")
        if output_root != "":
            makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_path, energy_window=energy_window, efermi=efermi)
        print(f"{dos.casename} case detected.")
        for atom in dos.atoms:
            print(f"  {len(dos.atom_numbers(atom))} {atom} detected")

        plot_custom(
            dos=dos,
            custom=custom,
            output_root=output_root,
            energy_window=energy_window,
            efermi=efermi,
            interactive=interactive,
            save_pickle=save_pickle,
            save_txt=save_txt,
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
