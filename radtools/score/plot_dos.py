from argparse import ArgumentParser
from os import makedirs
from os.path import abspath, join

from termcolor import cprint

from radtools.dos.dos import DOSQE, detect_seednames

from radtools.dos.plotting import COLOURS

from radtools.dos.pdos_plotting import (
    plot_custom_pdos,
    plot_orbital_resolved,
    plot_atom_resolved,
    plot_atom_to_total,
)


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
        if output_root != "":
            makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_path, energy_window=energy_window, efermi=efermi)
        print(f"{dos.casename} case detected.")
        for atom in dos.atoms:
            print(f"  {len(dos.atom_numbers(atom))} {atom} detected")

        pdos_parameters = {
            "output_root": output_root,
            "energy_window": energy_window,
            "dos_window": dos_window,
            "efermi": efermi,
            "relative": relative,
            "normalize": normalize,
            "interactive": interactive,
            "save_pickle": save_pickle,
            "save_txt": save_txt,
            "background_total": background_total,
            "colours": colours,
            "legend_fontsize": legend_fontsize,
            "axes_labels_fontsize": axes_labels_fontsize,
            "title_fontsize": title_fontsize,
        }

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
            cprint(f"Result is in {abspath(join(output_root, 'pdos-vs-dos'))}", "blue")

            # Plot PDOS for each atom/wfc
            plot_orbital_resolved(dos=dos, **pdos_parameters, separate=separate)

            # Plot wfc contribution into each atom
            plot_atom_resolved(dos=dos, **pdos_parameters, separate=separate)

            # Plot atom contributions into total PDOS
            plot_atom_to_total(dos=dos, **pdos_parameters, separate=separate)
        # Plot custom plot
        else:
            plot_custom_pdos(dos=dos, custom=custom, **pdos_parameters, labels=labels)


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
