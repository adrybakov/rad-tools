from argparse import ArgumentParser
import os

import matplotlib.pyplot as plt
import numpy as np
from termcolor import cprint

from radtools import __version__ as version
from radtools.io.internal import read_template
from radtools.io.tb2j import read_tb2j_model
from radtools.magnons.dispersion import MagnonDispersion
from radtools.decorate.stats import logo
from radtools.spinham.constants import TXT_FLAGS
from radtools.decorate.array import print_2d_array
from radtools.decorate.axes import plot_hlines


def manager(
    input_filename,
    template_file,
    output_name="magnon_dispersion",
    spin=None,
    spiral_vector=None,
    rotation_axis=None,
    k_path=None,
    form_model=False,
    R_vector=None,
    max_distance=None,
    min_distance=None,
    save_txt=False,
    interactive=False,
    verbose=False,
    bravais_type=None,
    join_output=False,
    nodmi=False,
    no_anisotropic=False,
):
    r"""
    :ref:`rad-plot-tb2j-magnons` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-plot-tb2j-magnons>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    head, _ = os.path.split(input_filename)
    out_head, out_tail = os.path.split(output_name)
    if len(out_head) == 0:
        out_head = head
    if len(out_tail) == 0:
        out_tail = "magnon_dispersion"

    output_name = os.path.join(out_head, out_tail)

    # Create the output directory if it does not exist
    if out_head != "":
        os.makedirs(out_head, exist_ok=True)

    # Check input parameters for consistency
    if form_model and template_file is None:
        raise ValueError(
            "--form-model option requires a template file (--template-filename)."
        )

    # Translate sequence of numbers to R vectors
    if R_vector is not None:
        R_vector = np.array(R_vector[: len(R_vector) // 3 * 3], dtype=int).reshape(
            (len(R_vector) // 3, 3)
        )
        R_vector = list(map(tuple, R_vector.tolist()))

    # Read the spinham
    spinham = read_tb2j_model(
        input_filename, quiet=not verbose, bravais_type=bravais_type
    )

    cprint(f"{spinham.variation} crystal detected", "green")

    # Get k points of the spinham
    kp = spinham.kpoints
    # Set custom k path
    if k_path is not None:
        kp.path = k_path

    if verbose:
        print("Predefined high symmetry k points:")
        for name in kp.hs_names:
            print(f"  {name} : {kp.hs_coordinates[name]}")

    # Force symmetry of the template
    if template_file is not None:
        template = read_template(template_file)
    else:
        template = None

    # Form the spinham based on the template or filter it
    if form_model:
        spinham.form_model(template=template)
    else:
        spinham.filter(
            min_distance=min_distance,
            max_distance=max_distance,
            R_vector=R_vector,
            template=template,
        )

    # Set the spin of the atoms
    if spin is not None:
        for i in range(len(spin) // 4):
            atom_name = spin[4 * i]
            atom = spinham.get_atom(atom_name)
            atom_spin = list(map(float, spin[4 * i + 1 : 4 * i + 4]))
            atom.spin_vector = atom_spin

    # Get the magnon dispersion
    dispersion = MagnonDispersion(
        spinham,
        Q=spiral_vector,
        n=rotation_axis,
        nodmi=nodmi,
        noaniso=no_anisotropic,
    )

    fig, ax = plt.subplots()

    omegas = dispersion(kp)

    ax.set_xticks(kp.coordinates(), kp.labels, fontsize=15)
    ax.set_ylabel("E, meV", fontsize=15)
    ax.vlines(
        kp.coordinates(),
        0,
        1,
        transform=ax.get_xaxis_transform(),
        colors="black",
        alpha=0.5,
        lw=1,
        ls="dashed",
    )
    for omega in omegas:
        ax.plot(kp.flatten_points(), omega)

    ax.set_xlim(kp.flatten_points()[0], kp.flatten_points()[-1])
    plot_hlines(ax, [0])

    if save_txt:
        main_separator = "=" * 80 + "\n"
        info = [
            main_separator,
            logo(date_time=True, line_length=80),
            f"\nMagnon dispersion is computed based on the file:\n{input_filename}\n",
        ]
        if template_file is not None:
            info.append(f"With template file:\n{template_file}\n")

        info.append(main_separator)
        info.append(TXT_FLAGS["cell"] + "\n")
        info.append(
            print_2d_array(
                spinham.cell,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=["x", "y", "z"],
            )
            + "\n"
        )
        info.append(f"Detected Bravais lattice: {spinham.variation}\n")
        info.append(main_separator)
        info.append(TXT_FLAGS["atoms"] + "\n")
        header_column = []
        header_row = [
            f"Index Name",
            "a1 (rel)",
            "a2 (rel)",
            "a3 (rel)",
            "S_x",
            "S_y",
            "S_z",
        ]
        atom_data = np.zeros((len(spinham.magnetic_atoms), 6))
        for a_i, atom in enumerate(spinham.magnetic_atoms):
            header_column.append(f"{atom.index:<5} {atom.name:<4}")
            atom_data[a_i, :3] = spinham.get_atom_coordinates(atom, relative=True)
            atom_data[a_i, 3:] = atom.spin_vector
        info.append(
            print_2d_array(
                atom_data,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=header_row,
                header_column=header_column,
            )
            + "\n"
        )

        if spiral_vector is not None:
            info.append(f"Spiral vector: {dispersion.Q} (relative: {spiral_vector})\n")
            info.append(f"Rotation axis: {dispersion.n}\n")

        info.append(main_separator)
        info.append(TXT_FLAGS["kpath"] + "\n")
        info.append(f"  {kp.path_string}\n")
        info.append(TXT_FLAGS["kpoints"] + "\n")
        names_data = []
        header_column = []
        for name in kp.hs_names:
            header_column.append(f"{name:<6}")
            names_data.append(kp.hs_coordinates[name])
        info.append(
            print_2d_array(
                names_data,
                borders=False,
                fmt=">.8f",
                print_result=False,
                header_column=header_column,
            )
            + "\n"
        )

        info.append(TXT_FLAGS["klabels"] + "\n")
        labels_data = []
        header_column = []
        for i in range(len(kp.labels)):
            header_column.append(f"{kp.labels[i]:<10}")
            labels_data.append([kp.coordinates()[i]])
        info.append(
            print_2d_array(
                labels_data,
                borders=False,
                fmt=">.8f",
                print_result=False,
                header_column=header_column,
            )
            + "\n"
        )
        info.append(main_separator)
        info.append(TXT_FLAGS["dispersion"] + "\n")

        header = ["Coordinate"]
        for i in range(omegas.shape[0]):
            header.append(f' {f"omega_{i}":>14}')
        header.append(f"    {'b1 (rel)':>10} {'b2 (rel)':>10} {'b3 (rel)':>10}")
        header.append(f"    {'k_x (abs)':>10} {'k_y (abs)':>10} {'k_z (abs)':>10}")
        info.append("".join(header))

        if not join_output:
            filename = f"{output_name}_data.txt"
            info.extend(
                [
                    "\n",
                    TXT_FLAGS["separate"],
                    "\n",
                    os.path.abspath(filename),
                ]
            )
            info = "".join(info)
            with open(f"{output_name}_info.txt", "w") as file:
                file.write(info)
            info = ""
            comments = "#"
        else:
            info = "".join(info)
            filename = f"{output_name}.txt"
            comments = ""

        np.savetxt(
            filename,
            np.concatenate(
                (
                    [kp.flatten_points()],
                    omegas,
                    kp.points(relative=True).T,
                    kp.points(relative=False).T,
                ),
                axis=0,
            ).T,
            fmt="%.8f"
            + " %.8e" * omegas.shape[0]
            + "   "
            + " %.8f" * 3
            + "   "
            + " %.8f" * 3,
            header=info,
            comments=comments,
        )
    if interactive:
        plt.show()
    else:
        plt.savefig(
            f"{output_name}.png",
            bbox_inches="tight",
            dpi=600,
        )
        cprint(
            f"Results are in {os.path.abspath(out_head)}, seedname: {out_tail}.", "blue"
        )


def create_parser():
    parser = ArgumentParser(
        description="Script for magnon dispersion from TB2J results."
    )

    parser.add_argument(
        "-if",
        "--input-filename",
        metavar="filename",
        type=str,
        required=True,
        help="Relative or absolute path to the 'exchange.out' file, "
        + "including the name and extension of the file itself.",
    )
    parser.add_argument(
        "-tf",
        "--template-file",
        metavar="filename",
        type=str,
        default=None,
        help="Relative or absolute path to the template file, "
        + "including the name and extension of the file.",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        metavar="filename",
        type=str,
        default="magnon_dispersion",
        help="Seedname for the output files.",
    )
    parser.add_argument(
        "-s",
        "--spin",
        metavar="Atom S_x S_y S_z",
        type=str,
        nargs="*",
        default=None,
        help="Spin of the atoms in the spinham.",
    )
    parser.add_argument(
        "-Q",
        "--spiral-vector",
        metavar=("Q_x", "Q_y", "Q_z"),
        type=float,
        nargs=3,
        default=None,
        help="Spin spiral vector. relative to the reciprocal cell.",
    )
    parser.add_argument(
        "-ra",
        "--rotation-axis",
        metavar=("n_x", "n_y", "n_z"),
        type=float,
        nargs=3,
        default=None,
        help="Direction of global rotation axis. In absolute coordinates in real space.",
    )
    parser.add_argument(
        "-kp",
        "--k-path",
        metavar="k-path",
        type=str,
        default=None,
        help="Path in reciprocal space for the magnon dispersion.",
    )
    parser.add_argument(
        "-fm",
        "--form-model",
        action="store_true",
        default=False,
        help="Whether to form the spinham based on the template.",
    )
    parser.add_argument(
        "-R",
        "--R-vector",
        metavar="i j k",
        type=int,
        nargs="*",
        default=None,
        help="R vectors for filtering the spinham.",
    )
    parser.add_argument(
        "-maxd",
        "--max-distance",
        metavar="distance",
        type=float,
        default=None,
        help="(<=) Maximum distance.",
    )
    parser.add_argument(
        "-mind",
        "--min-distance",
        metavar="distance",
        type=float,
        default=None,
        help="(>=) Minimum distance.",
    )
    parser.add_argument(
        "-st",
        "--save-txt",
        action="store_true",
        default=False,
        help="Whether to save data to .txt file.",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        default=False,
        help="Interactive plotting.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )
    parser.add_argument(
        "-bt",
        "--bravais-type",
        metavar="type",
        choices=[
            "CUB",
            "FCC",
            "BCC",
            "TET",
            "BCT",
            "ORC",
            "ORCF",
            "ORCI",
            "ORCC",
            "HEX",
            "RHL",
            "MCL",
            "MCLC",
            "TRI",
        ],
        type=str,
        default=None,
        help="Bravais lattice type. If not provided, the type is identified automatically.",
    )
    parser.add_argument(
        "-jo",
        "--join-output",
        action="store_true",
        default=None,
        help="Whether to join the output files into a single file.",
    )
    parser.add_argument(
        "-nodmi",
        action="store_true",
        default=False,
        help="Whether to ignore DMI in the spinham.",
    )
    parser.add_argument(
        "-noa",
        "--no-anisotropic",
        action="store_true",
        default=False,
        help="Whether to ignore anisotropic symmetric exchange in the spinham.",
    )

    return parser
