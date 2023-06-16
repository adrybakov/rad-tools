from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime
from os import makedirs
from os.path import abspath, split

import numpy as np
from termcolor import cprint

from radtools import __version__ as version
from radtools.io.tb2j import read_tb2j_model


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
    """

    # Create the output directory if it does not exist
    makedirs(split(output_name)[0], exist_ok=True)

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
        "=" * 20
        + "\n"
        + "Neighbors template:\n"
        + "i j R_a R_b R_c\n"
        + "-" * 20
        + "\n"
        + "J1 $J_1$\n"
        + "atom1 atom2  0  0  0\n"
        + "atom1 atom2  1  0  0\n"
        + "atom1 atom1 -1  0  2\n"
        + "-" * 20
        + "\n"
        + "J2\n"
        + "atom2 atom1  9  5 -3\n"
        + "atom1 atom2  1  4  0\n"
        + "atom2 atom2  1  0  2\n"
        + "=" * 20
        + "\n"
    )
    with open(output_name, "w") as file:
        # Write the draft
        if input_filename is None:
            file.write(
                f"Template is created "
                + f"on {cd.day} {month_name[cd.month]} {cd.year}"
                + f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools {version}\n\n"
            )

            file.write(template)
        # Create the template based on the input file
        else:
            file.write(
                f"Template is created based on the file: {input_filename}\n"
                + f"on {cd.day} {month_name[cd.month]} {cd.year}"
                + f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools {version}\n\n"
            )

            # Read and filter the model
            model = read_tb2j_model(input_filename, quiet=not verbose)
            model.filter(
                min_distance=min_distance, max_distance=max_distance, R_vector=R_vector
            )

            # Write header
            file.write(
                "=" * 20
                + "\n"
                + "Neighbors template:\n"
                + "i j R_a R_b R_c\n"
                + "-" * 20
                + "\n"
            )

            # Get bonds from the model
            data = []
            for atom1, atom2, R, J in model:
                data.append(
                    (atom1, atom2, R, model.crystal.get_distance(atom1, atom2, R))
                )

            # Sort bonds by distance
            data.sort(key=lambda x: x[3])

            j = 1
            file.write(f"J{j} " + "$J_{" + f"{j}" + "}$\n")
            file.write(
                f"{data[0][0]:4} {data[0][1]:4} "
                + f"{data[0][2][0]:3.0f} {data[0][2][1]:3.0f} {data[0][2][2]:3.0f}\n"
            )
            for i, (atom1, atom2, R, distance) in enumerate(data[1:]):
                # If distance is the same as the previous one, write the bond
                if abs(distance - data[i][3]) < eps:
                    file.write(
                        f"{atom1:4} {atom2:4} "
                        + f"{R[0]:3.0f} {R[1]:3.0f} {R[2]:3.0f}\n"
                    )
                # If distance is different, start a new group and write the bond
                else:
                    j += 1
                    file.write("-" * 20 + "\n")
                    file.write(f"J{j} " + "$J_{" + f"{j}" + "}$\n")
                    file.write(
                        f"{atom1:4} {atom2:4} "
                        + f"{R[0]:3.0f} {R[1]:3.0f} {R[2]:3.0f}\n"
                    )

            file.write("=" * 20 + "\n")
    cprint(
        f"Template draft is in " + f"{abspath(output_name)}, grouped by distance",
        "green",
    )
    cprint(f"Do not forget to correct the template draft to your needs!", "yellow")


def create_parser():
    parser = ArgumentParser(
        description="Script for the creation of template`s template"
    )

    parser.add_argument(
        "-on",
        "--output-name",
        metavar="filename",
        type=str,
        default="template.txt",
        help="Name for the template output file.",
    )
    parser.add_argument(
        "-if",
        "--input-filename",
        metavar="filename",
        type=str,
        default=None,
        help="Relative or absolute path to the 'exchange.out' file, "
        + "including the name and extension of the file itself.",
    )
    parser.add_argument(
        "-R",
        "--R-vector",
        metavar="R1a R1b R1c",
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
        "-d",
        "--distance",
        metavar="distance",
        type=float,
        default=None,
        help="(=) Exact distance.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=1e-3,
        metavar="value",
        help="Epsilon for the distance comparison.",
    )

    return parser
