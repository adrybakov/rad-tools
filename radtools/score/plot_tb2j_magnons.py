from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime
from os import makedirs
from os.path import join, abspath

from termcolor import cprint
import numpy as np
import matplotlib.pyplot as plt

from radtools import __version__ as version
from radtools.io.internal import read_template
from radtools.io.tb2j import read_tb2j_model
from radtools.magnons.magnons import MagnonDispersion


def manager(
    input_filename,
    template_file,
    output_name="magnon_dispersion",
    output_path=".",
    spin=None,
    spiral_vector=None,
    rotation_axis=None,
    path=None,
    force_symmetry=False,
    R_vector=None,
    max_distance=None,
    min_distance=None,
    save_txt=False,
    interactive=False,
    verbose=False,
):
    r"""
    :ref:`rad-plot-tb2j-magnons` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-plot-tb2j-magnons>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    # Create the output directory if it does not exist
    makedirs(output_path, exist_ok=True)

    # Get current date and time
    cd = datetime.now()

    # Check input parameters for consistency
    if force_symmetry and template_file is None:
        raise ValueError("--force-symmetry option requires a template file.")

    # Translate sequence of numbers to R vectors
    if R_vector is not None:
        R_vector = np.array(R_vector[: len(R_vector) // 3 * 3], dtype=int).reshape(
            (len(R_vector) // 3, 3)
        )
        R_vector = list(map(tuple, R_vector.tolist()))

    # Read the model
    model = read_tb2j_model(input_filename, quiet=not verbose)

    cprint(f"{model.variation} crystal detected", "green")

    # Get k points of the model
    kp = model.crystal.get_kpoints()  # Set custom k path
    if path is not None:
        kp.path = path

    if verbose:
        print("Predefined high symmetry k points:")
        for name in kp._hs_points:
            print(f"  {name} : {kp._hs_coordinates[name]}")

    # Force symmetry of the template
    if template_file is not None:
        template = read_template(template_file)
    else:
        template = None
    if force_symmetry:
        model.force_symmetry(template=template)

    # Filter the model
    model.filter(
        min_distance=min_distance,
        max_distance=max_distance,
        R_vector=R_vector,
        template=template,
    )

    if spin is not None:
        for i in range(len(spin) // 4):
            atom_name = spin[4 * i]
            atom = model.crystal.get_atom(atom_name)
            atom_spin = list(map(float, spin[4 * i + 1 : 4 * i + 4]))
            atom.spin_vector = atom_spin

    # Get the magnon dispersion
    dispersion = MagnonDispersion(model, Q=spiral_vector, n=rotation_axis)

    fig, ax = plt.subplots()

    omegas = []
    for point in kp.points:
        omegas.append(dispersion.omega(point))

    omegas = np.array(omegas).T

    ax.set_xticks(kp.coordinates, kp.labels, fontsize=15)
    ax.set_ylabel("E, meV", fontsize=15)
    ax.vlines(
        kp.coordinates,
        0,
        1,
        transform=ax.get_xaxis_transform(),
        colors="black",
        alpha=0.5,
        lw=1,
        ls="dashed",
    )
    for omega in omegas:
        ax.plot(kp.flatten_points, omega)

    ax.set_xlim(kp.flatten_points[0], kp.flatten_points[-1])
    ax.set_ylim(0, None)

    if save_txt:
        header = (
            f"Magnon dispersion is computed based on the file: {input_filename}\n"
            + f"on {cd.day} {month_name[cd.month]} {cd.year}"
            + f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools {version}\n\n"
        )
        header += f"\nSpins:\n"
        for atom in model.magnetic_atoms:
            header += f"  {atom.name} : {atom.spin_vector}\n"

        header += f"\nDetected Bravais lattice: {model.variation}\n"

        if spiral_vector is not None:
            header += f"\nSpiral vector: {dispersion.Q} (relative: {spiral_vector})\n"
            header += f"Rotation axis: {dispersion.n}\n"

        header += f"\nK path: {kp.path_string}\n" + f"\nK points: \n"
        for name in kp._hs_points:
            header += f"  {name} : {kp._hs_coordinates[name]}\n"

        header += "\nLabels: \n"
        for i in range(len(kp.labels)):
            header += f"  {kp.labels[i]} : {kp.coordinates[i]:.8f}\n"

        header += "\nMagnon dispersion header:\ncoordinate"
        for i in range(omegas.shape[0]):
            header += f' {f"omega_{i}":>14}'
        header += f"    {'k_x':>10} {'k_y':>10} {'k_z':>10}"

        with open(join(output_path, f"{output_name}_info.txt"), "w") as file:
            file.write(header)

        np.savetxt(
            join(output_path, f"{output_name}.txt"),
            np.concatenate(([kp.flatten_points], omegas, kp.points.T), axis=0).T,
            fmt="%.8f" + " %.8e" * omegas.shape[0] + "   " + " %.8f" * 3,
        )
    if interactive:
        plt.show()
    else:
        plt.savefig(
            join(output_path, f"{output_name}.png"),
            bbox_inches="tight",
            dpi=600,
        )
        cprint(
            f"Results are in {abspath(output_name)}, seedname: {output_name}.", "blue"
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
        "-op",
        "--output-path",
        metavar="path",
        type=str,
        default=".",
        help="Relative or absolute path to the folder for saving outputs.",
    )
    parser.add_argument(
        "-s",
        "--spin",
        metavar="Atom S_x S_y S_z",
        type=str,
        nargs="*",
        default=None,
        help="Spin of the atoms in the model.",
    )
    parser.add_argument(
        "-Q",
        "--spiral-vector",
        metavar="Q_i Q_j Q_k",
        type=float,
        nargs=3,
        default=None,
        help="Spin spiral vector. relative to the reciprocal cell.",
    )
    parser.add_argument(
        "-ra",
        "--rotation-axis",
        metavar="n_i n_j n_k",
        type=float,
        nargs=3,
        default=None,
        help="Direction of global rotation axis. In absolute coordinates in real space.",
    )
    parser.add_argument(
        "-p",
        "--path",
        metavar="k-path",
        type=str,
        default=None,
        help="Path in reciprocal space for the magnon dispersion.",
    )
    parser.add_argument(
        "-fs",
        "--force-symmetry",
        action="store_true",
        default=False,
        help="Whether to force the symmetry of the template on the model.",
    )
    parser.add_argument(
        "-R",
        "--R-vector",
        metavar="i j k",
        type=int,
        nargs="*",
        default=None,
        help="R vectors for filtering the model.",
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

    return parser
