from argparse import ArgumentParser
import os

from termcolor import cprint

from radtools.dos.dos import DOSQE, detect_seednames
from radtools.dos.fatbands_plotting import plot_custom_fatbands
from radtools.dos.plotting import COLOURS


def manager(
    input_folder=".",
    seedname=None,
    output_name="",
    energy_window=None,
    k_window=None,
    efermi=0.0,
    verbose=False,
    interactive=False,
    separate=False,
    save_pickle=False,
    save_txt=False,
    custom=None,
    colours=None,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
    k_points=None,
):
    r"""
    :ref:`rad-plot-fatbands` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-plot-fatbands>`.
    Parameters of the function directly
    correspond to the arguments of the script.

    Parameters
    ----------
    input_folder : str, default "."
        Relative or absolute path to the folder with PDOS files.

        Console argument: ``-if`` / ``--input-folder``

        Metavar: "path"
    seedname : str, optional
        Prefix for input files with PDOS(E).

        In the case of Quantum Espresso-produced seedname is the same
        as specified in the QE projwfc.x input file (filpdos).

        If it is not provided the script tries to
        detect it automatically in the
        ``rad-plot-fatbands_input-folder`` folder.

        Console argument: ``-s`` / ``--seedname``

        Metavar: "name"
    output_name : str, optional
        Relative or absolute path to the folder for saving outputs.

        Console argument: ``-on`` / ``--output-name``

        Metavar: "path"
    energy_window : tuple of 2 float, optional
        Energy window for the plots.

        By default the whole energy range present in the files is plotted.

        Console argument: ``-ew`` / ``--energy-window``

        Metavar: ("min", "max")
    k_window : tuple of 2 float, optional
        K-point window for the plots.

        By default the whole states/eV range is plotted.

        Console argument: ``-kw`` / ``--k-window``

        Metavar: ("min", "max")
    efermi : float, default 0.0
        Fermi energy.

        Zero is shifted to Fermi energy.

        Console argument: ``-ef`` / ``--efermi``

        Metavar: "energy"
    verbose : bool, default False
        Verbose output, propagates to the called methods.

        Console argument: ``-v`` / ``--verbose``
    interactive : bool, default False
        Interactive plotting.

        Console argument: ``-i`` / ``--interactive``
    separate : bool, default False
        Whether to plot projected DOS for each atom of the same type separately.

        Console argument: ``-sep`` / ``--separate``
    save_pickle : bool, default False
        Whether to save figures as .pickle files.

        Console argument: ``-sp`` / ``--save-pickle``
    save_txt : bool, default False
        Whether to save some data as txt files.

        Console argument: ``-st`` / ``--save-txt``
    custom : list of str, optional
        Custom PDOS plot.

        Console argument: ``--custom``

        Metavar: "description"
    colours : list of str, optional
        Colours for the relative and custom plots.

        Values are passed directly to the matplotlib as strings,
        therefore any valid value is allowed. Examples: "red" or "#FF0000".
        When ``custom`` is used the order of colours is the same as for
        the values of the ``custom``.

        Console argument: ``-cls`` / ``--colours``
    labels : list of str, optional
        Labels for the custom plots.

        Amount of labels have to be the same as the amount of ``custom`` strings, or one more.
        If one more, then first one is interpreted as the label for the background
        (Use "None" to switch it off). If the amount of argument is one more  and the first one is None,
        then the label for the total PDOS is switched off and the total PDOS itself is not plotted.

        Console argument: ``-lbs`` / ``--labels``
    legend_fontsize : int, default 12
        Fontsize of the legend.

        Console argument: ``-lfs`` / ``--legend-fontsize``

        Metavar: "fontsize"
    axes_labels_fontsize : int, default 14
        Fontsize of the labes of the axes.

        Console argument: ``-alfs`` / ``--axes-labels-fontsize``

        Metavar: "fontsize"
    title_fontsize : int, default 18
        Fontsize of the title.

        Console argument: ``-tfs`` / ``--title-fontsize``

        Metavar: "fontsize"
    k_points : list of str, optional
        Plot coordinates of high symmetry points.

        Console argument: ``-kp`` / ``--k-points``

        Metavar: "G 0.23 K 0.64 Y 1.2 ..."
    """

    out_head, out_tail = os.path.split(output_name)
    if len(out_head) == 0:
        out_head = input_folder

    if colours is None:
        colours = COLOURS

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
        output_root = os.path.join(out_head, f"{out_tail}{seedname}")
        if output_root != "":
            os.makedirs(output_root, exist_ok=True)

        # Load DOS data.
        dos = DOSQE(seedname, input_folder, energy_window=energy_window, efermi=efermi)
        print(f"{dos.casename} case detected.")
        for atom in dos.atoms:
            print(f"  {len(dos.atom_numbers(atom))} {atom} detected")

        plot_custom_fatbands(
            dos=dos,
            custom=custom,
            output_root=output_root,
            energy_window=energy_window,
            k_window=k_window,
            efermi=efermi,
            interactive=interactive,
            separate=separate,
            save_pickle=save_pickle,
            save_txt=save_txt,
            colours=colours,
            labels=labels,
            legend_fontsize=legend_fontsize,
            axes_labels_fontsize=axes_labels_fontsize,
            title_fontsize=title_fontsize,
            k_points=k_points,
        )


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
        "-kw",
        "--k-window",
        default=None,
        metavar=("min", "max"),
        type=float,
        nargs=2,
        help='K-point window for the plots.',
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
        "-sep",
        "--separate",
        default=False,
        action="store_true",
        help='Whether to plot projected DOS for each atom of the same type separately.',
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
        "--custom",
        default=None,
        metavar="description",
        type=str,
        nargs="*",
        help='Custom PDOS plot.',
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
    parser.add_argument(
        "-kp",
        "--k-points",
        default=None,
        metavar="G 0.23 K 0.64 Y 1.2 ...",
        type=str,
        nargs="*",
        help='Plot coordinates of high symmetry points.',
    )

    return parser
