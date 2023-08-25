from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime
import os

import numpy as np
from termcolor import cprint

from radtools import __version__ as version
from radtools.io.tb2j import load_tb2j_model
from radtools.decorate.stats import logo


def manager(
    output_name="template.txt",
    input_filename=None,
    R_vector=None,
    max_distance=None,
    min_distance=None,
    distance=None,
    verbose=False,
    eps=1e-3,
):
    r"""
    :ref:`rad-make-template` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-make-template>`.
    Parameters of the function directly
    correspond to the arguments of the script.

    Parameters
    ----------
    output_name : str, default "template.txt"
        Name for the template output file.

        See also: :ref:`example <output-notes>`.

        Console argument: ``-on`` / ``--output-name``

        Metavar: "filename"
    input_filename : str, optional
        Relative or absolute path to the "exchange.out" file, including name and extension of the file.

        Console argument: ``-if`` / ``--input-filename``

        Metavar: "filename"

        .. versionchanged:: 0.5.12 Renamed from "tb2j_filename"
    R_vector : list of int, optional
        R vectors for filtering the spin Hamiltonian.

        In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to).
        Atom 1 is always located in (0, 0, 0) unit cell, while atom 2 is located in
        R = (i, j, k) unit cell. This parameter tells the script to keep only the
        bonds for which atom 2 is located in one of specified R supercells.
        Supercells are specified by a set of integers separated by spaces.
        They are grouped by three starting from the left and forms a set
        of R vectors. If the last group contains 1 or 2 integers they are ignored.

        Console argument: ``-R`` / ``--R-vector``

        Metavar: "i1 j1 k1 i2 j2 k2 ..."
    max_distance : float, optional
        (<=) Maximum distance.

        All the bonds with the distance between atom 1 and atom 2
        greater than maximum distance are excluded from the model.

        Console argument: ``-maxd`` / ``--max-distance``

        Metavar: "distance"
    min_distance : float, optional
        (>=) Minimum distance.

        All the bonds with the distance between atom 1 and atom 2
        lower than minimum distance are excluded from the Hamiltonian.

        Console argument: ``-mind`` / ``--min-distance``

        Metavar: "distance"
    distance : float, optional
        (=) Exact distance.

        Only the bonds with the exact distance remains in the model.
        There is no point in specifying maximum or minimum distance when
        this parameter is provided.

        Console argument: ``-d`` / ``--distance``

        Metavar: "distance"
    verbose : bool, default False
        Verbose output, propagates to the called methods.

        Console argument: ``-v`` / ``--verbose``
    eps : float, default 1e-3
        Epsilon for the distance comparison.

        Console argument: ``--eps``
    """

    n_sep = 80

    out_head, out_tail = os.path.split(output_name)
    if len(out_tail) == 0:
        out_tail = "template.txt"

    output_name = os.path.join(out_head, out_tail)

    # Create the output directory if it does not exist
    if out_head != "":
        os.makedirs(out_head, exist_ok=True)

    # Get current date and time
    cd = datetime.now()

    if distance is not None:
        min_distance = distance
        max_distance = distance

    # Translate sequence of numbers to R vectors
    if R_vector is not None:
        R_vector = np.array(R_vector[: len(R_vector) // 3 * 3], dtype=int).reshape(
            (len(R_vector) // 3, 3)
        )
        R_vector = list(map(tuple, R_vector.tolist()))

    # Template draft
    template = (
        "=" * n_sep
        + "\n"
        + "Neighbors template:\n"
        + f"Atom1 Atom2 {'i':>3} {'j':>3} {'k':>3}\n"
        + "-" * n_sep
        + "\n"
        + "J1 $J_1$\n"
        + f"atom1 atom2 {0:3} {0:3} {0:3}\n"
        + f"atom1 atom2 {1:3} {0:3} {0:3}\n"
        + f"atom1 atom1 {-1:3} {0:3} {2:3}\n"
        + "-" * n_sep
        + "\n"
        + "J2\n"
        + f"atom2 atom1 {9:3} {5:3} {-3:3}\n"
        + f"atom1 atom2 {1:3} {4:3} {0:3}\n"
        + f"atom2 atom2 {1:3} {0:3} {2:3}\n"
        + "=" * n_sep
        + "\n"
    )
    with open(output_name, "w") as file:
        file.write("=" * n_sep + "\n" + logo(date_time=True, line_length=n_sep) + "\n")
        # Write the draft
        if input_filename is None:
            file.write(template)
        # Create the template based on the input file
        else:
            file.write(
                f"Template is created based on the file:\n{os.path.abspath(input_filename)}\n"
            )

            # Read and filter the model
            model = load_tb2j_model(input_filename, quiet=not verbose)
            model.filter(
                min_distance=min_distance, max_distance=max_distance, R_vector=R_vector
            )

            # Write header
            file.write(
                "=" * n_sep
                + "\n"
                + "Neighbors template:\n"
                + f"Atom1 Atom2 {'i':>3} {'j':>3} {'k':>3}\n"
                + "-" * n_sep
                + "\n"
            )

            # Get bonds from the model
            data = []
            for atom1, atom2, R, J in model:
                data.append(
                    (
                        atom1.name,
                        atom2.name,
                        R,
                        model.get_distance(atom1, atom2, R),
                    )
                )

            # Sort bonds by distance
            data.sort(key=lambda x: x[3])

            j = 1
            file.write(f"J{j} " + "$J_{" + f"{j}" + "}$\n")
            file.write(
                f"{data[0][0]:5} {data[0][1]:5} "
                + f"{data[0][2][0]:3.0f} {data[0][2][1]:3.0f} {data[0][2][2]:3.0f}\n"
            )
            for i, (atom1, atom2, R, distance) in enumerate(data[1:]):
                # If distance is the same as the previous one, write the bond
                if abs(distance - data[i][3]) < eps:
                    file.write(
                        f"{atom1:5} {atom2:5} "
                        + f"{R[0]:3.0f} {R[1]:3.0f} {R[2]:3.0f}\n"
                    )
                # If distance is different, start a new group and write the bond
                else:
                    j += 1
                    file.write("-" * n_sep + "\n")
                    file.write(f"J{j} " + "$J_{" + f"{j}" + "}$\n")
                    file.write(
                        f"{atom1:5} {atom2:5} "
                        + f"{R[0]:3.0f} {R[1]:3.0f} {R[2]:3.0f}\n"
                    )

            file.write("=" * n_sep + "\n")
    cprint(
        f"Template draft is in "
        + f"{os.path.abspath(output_name)}, grouped by distance",
        "blue",
    )
    cprint(f"Do not forget to correct the template draft to your needs!", "yellow")


def create_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-on",
        "--output-name",
        default="template.txt",
        metavar="filename",
        type=str,
        help="Name for the template output file.",
    )
    parser.add_argument(
        "-if",
        "--input-filename",
        default=None,
        metavar="filename",
        type=str,
        help='Relative or absolute path to the "exchange.out" file, including name and extension of the file.',
    )
    parser.add_argument(
        "-R",
        "--R-vector",
        default=None,
        metavar="i1 j1 k1 i2 j2 k2 ...",
        type=int,
        nargs="*",
        help="R vectors for filtering the spin Hamiltonian.",
    )
    parser.add_argument(
        "-maxd",
        "--max-distance",
        default=None,
        metavar="distance",
        type=float,
        help="(<=) Maximum distance.",
    )
    parser.add_argument(
        "-mind",
        "--min-distance",
        default=None,
        metavar="distance",
        type=float,
        help="(>=) Minimum distance.",
    )
    parser.add_argument(
        "-d",
        "--distance",
        default=None,
        metavar="distance",
        type=float,
        help="(=) Exact distance.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Verbose output, propagates to the called methods.",
    )
    parser.add_argument(
        "--eps",
        default=1e-3,
        type=float,
        help="Epsilon for the distance comparison.",
    )

    return parser
