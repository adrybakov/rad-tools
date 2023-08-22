from argparse import ArgumentParser
import os

from termcolor import cprint

from radtools.dos.dos import DOSQE, detect_seednames
from radtools.dos.pdos_plotting import (
    plot_atom_resolved,
    plot_atom_to_total,
    plot_custom_pdos,
    plot_orbital_resolved,
)
from radtools.dos.plotting import COLOURS


def manager(
    input_folder=".",
    seedname=None,
    output_name="",
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

    Parameters
    ----------
    input_folder : str, default "."
        Relative or absolute path to the folder with PDOS files.

        Console argument: ``-if`` / ``--input-folder``

        Metavar: "path"

        .. versionchanged:: 0.8.0 Renamed from ``input_path``
    seedname : str, optional
        Prefix for input files with PDOS(E).

        In the case of Quantum Espresso-produced seedname is the same
        as specified in the QE projwfc.x input file (filpdos).

        If it is not provided the script tries to
        detect it automatically in the
        ``rad-plot-dos_`input-folder`` folder.

        Console argument: ``-s`` / ``--seedname``

        Metavar: "name"

        .. versionchanged:: 0.5.21 from "filpdos" to "seedname".
    output_name : str, optional
        Relative or absolute path to the folder for saving outputs.

        Console argument: ``-on`` / ``--output-name``

        Metavar: "path"
    energy_window : tuple of 2 float, optional
        Energy window for the plots.

        By default the whole energy range present in the files is plotted.

        Console argument: ``-ew`` / ``--energy-window``

        Metavar: ("min", "max")

        .. versionchanged:: 0.5.21 Renamed from "window" to "energy-window".
    dos_window : tuple of 2 float, optional
        DOS window for the plots.

        By default the whole states/eV range is plotted.

        Console argument: ``-dw`` / ``--dos-window``

        Metavar: ("min", "max")

        .. versionadded:: 0.5.21
    efermi : float, default 0.0
        Fermi energy.

        Zero is shifted to Fermi energy.

        Console argument: ``-ef`` / ``--efermi``

        Metavar: "energy"
    separate : bool, default False
        Whether to plot projected DOS for each atom of the same type separately.

        Console argument: ``-sep`` / ``--separate``
    relative : bool, default False
        Whether to use relative style.

        Console argument: ``-r`` / ``--relative``
    normalize : bool, default False
        Whether to normalized PDOS values to 1.

        Console argument: ``-n`` / ``--normalize``

        (with respect to LDOS of each plot or to total PDOS if
        :ref:`rad-plot-dos_background-total` is used).
    verbose : bool, default False
        Verbose output, propagates to the called methods.

        Console argument: ``-v`` / ``--verbose``
    interactive : bool, default False
        Interactive plotting.

        Console argument: ``-i`` / ``--interactive``
    save_pickle : bool, default False
        Whether to save figures as .pickle files.

        Console argument: ``-sp`` / ``--save-pickle``

        .. versionadded:: 0.5.21
    save_txt : bool, default False
        Whether to save some data as txt files.

        It does not affect "pdos-vs-dos.png",
        because these data are accessible directly from PDOS input files.

        Console argument: ``-st`` / ``--save-txt``

        .. versionadded:: 0.5.21
    background_total : bool, default False
        Whether to use total PDOS as the background for all plots.

        Console argument: ``-bt`` / ``--background-total``

        .. versionadded:: 0.5.21
    custom : list of str, optional
        Custom PDOS plot. See :ref:`rad-plot-dos_custom-plots` for info.

        Console argument: ``--custom``

        Metavar: "description"

        .. versionadded:: 0.7.5
    colours : list of str, optional
        Colours for the relative and custom plots.

        Values are passed directly to the matplotlib as strings,
        therefore any valid value is allowed. Examples: "red" or "#FF0000".
        When ``custom`` is used the order of colours is the same as for
        the values of the ``custom``.

        Console argument: ``-cls`` / ``--colours``

        .. versionadded:: 0.7.5
    labels : list of str, optional
        Labels for the custom plots.

        Amount of labels have to be the same as the amount of ``custom`` strings, or one more.
        If one more, then first one is interpreted as the label for the background
        (Use "None" to switch it off). If the amount of argument is one more  and the first one is None,
        then the label for the total PDOS is switched off and the total PDOS itself is not plotted.

        Console argument: ``-lbs`` / ``--labels``

        .. versionadded:: 0.7.6
    legend_fontsize : int, default 12
        Fontsize of the legend.

        Console argument: ``-lfs`` / ``--legend-fontsize``

        Metavar: "fontsize"

        .. versionadded:: 0.7.8
    axes_labels_fontsize : int, default 14
        Fontsize of the labes of the axes.

        Console argument: ``-alfs`` / ``--axes-labels-fontsize``

        Metavar: "fontsize"

        .. versionadded:: 0.7.8
    title_fontsize : int, default 18
        Fontsize of the title.

        Console argument: ``-tfs`` / ``--title-fontsize``

        Metavar: "fontsize"

        .. versionadded:: 0.7.8

    """

    out_head, out_tail = os.path.split(output_name)
    if len(out_head) == 0:
        out_head = input_folder

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
        seednames = detect_seednames(input_folder)
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
        output_root = os.path.join(out_head, f"{out_tail}{seedname}{suffix}")
        if output_root != "":
            os.makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_folder, energy_window=energy_window, efermi=efermi)
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
                output_name=os.path.join(output_root, "pdos-vs-dos"),
                interactive=interactive,
                efermi=efermi,
                xlim=energy_window,
                ylim=dos_window,
                save_pickle=save_pickle,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )
            cprint(
                f"Result is in {os.path.abspath(os.path.join(output_root, 'pdos-vs-dos'))}",
                "blue",
            )

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

    parser = ArgumentParser()
    parser.add_argument(
        "-if",
        "--input-folder",
        default=".",
        metavar="path",
        type=str,
        help='Relative or absolute path to the folder with PDOS files.',
    )
    parser.add_argument(
        "-s",
        "--seedname",
        default=None,
        metavar="name",
        type=str,
        help='Prefix for input files with PDOS(E).',
    )
    parser.add_argument(
        "-on",
        "--output-name",
        default="",
        metavar="path",
        type=str,
        help='Relative or absolute path to the folder for saving outputs.',
    )
    parser.add_argument(
        "-ew",
        "--energy-window",
        default=None,
        metavar=("min", "max"),
        type=float,
        nargs=2,
        help='Energy window for the plots.',
    )
    parser.add_argument(
        "-dw",
        "--dos-window",
        default=None,
        metavar=("min", "max"),
        type=float,
        nargs=2,
        help='DOS window for the plots.',
    )
    parser.add_argument(
        "-ef",
        "--efermi",
        default=0.0,
        metavar="energy",
        type=float,
        help='Fermi energy.',
    )
    parser.add_argument(
        "-sep",
        "--separate",
        default=False,
        action="store_true",
        help='Whether to plot projected DOS for each atom of the same type separately.',
    )
    parser.add_argument(
        "-r",
        "--relative",
        default=False,
        action="store_true",
        help='Whether to use relative style.',
    )
    parser.add_argument(
        "-n",
        "--normalize",
        default=False,
        action="store_true",
        help='Whether to normalized PDOS values to 1.',
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help='Verbose output, propagates to the called methods.',
    )
    parser.add_argument(
        "-i",
        "--interactive",
        default=False,
        action="store_true",
        help='Interactive plotting.',
    )
    parser.add_argument(
        "-sp",
        "--save-pickle",
        default=False,
        action="store_true",
        help='Whether to save figures as .pickle files.',
    )
    parser.add_argument(
        "-st",
        "--save-txt",
        default=False,
        action="store_true",
        help='Whether to save some data as txt files.',
    )
    parser.add_argument(
        "-bt",
        "--background-total",
        default=False,
        action="store_true",
        help='Whether to use total PDOS as the background for all plots.',
    )
    parser.add_argument(
        "--custom",
        default=None,
        metavar="description",
        type=str,
        nargs="*",
        help='Custom PDOS plot. See :ref:`rad-plot-dos_custom-plots` for info.',
    )
    parser.add_argument(
        "-cls",
        "--colours",
        default=None,
        type=str,
        nargs="*",
        help='Colours for the relative and custom plots.',
    )
    parser.add_argument(
        "-lbs",
        "--labels",
        default=None,
        type=str,
        nargs="*",
        help='Labels for the custom plots.',
    )
    parser.add_argument(
        "-lfs",
        "--legend-fontsize",
        default=12,
        metavar="fontsize",
        type=int,
        help='Fontsize of the legend.',
    )
    parser.add_argument(
        "-alfs",
        "--axes-labels-fontsize",
        default=14,
        metavar="fontsize",
        type=int,
        help='Fontsize of the labes of the axes.',
    )
    parser.add_argument(
        "-tfs",
        "--title-fontsize",
        default=18,
        metavar="fontsize",
        type=int,
        help='Fontsize of the title.',
    )

    return parser
