#! /usr/local/bin/python3

from argparse import ArgumentParser
import os

import numpy as np
from termcolor import cprint


def manager(input_filename, span=0.1, output_name=""):
    r"""
    :ref:`rad-identify-wannier-centres` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-identify-wannier-centres>`.
    Parameters of the function directly
    correspond to the arguments of the script.

    Parameters
    ----------
    input_filename : str
        Relative or absolute path to the "_centres.xyz" file

        Identified Wannier centres are stored in the "filename_identified" file.

        Console argument: ``-if`` / ``--input-filename``

        Metavar: "filename"

        .. versionchanged:: 0.8.0 Renamed from ``filename``
    span : float, default 0.1
        Distance tolerance between centre and atom. (in Angstroms)

        If some centres remain unidentified try to increase the span.

        Console argument: ``-s`` / ``--span``
    output_name : str, default ""
        Seedname for the output files.

        See also: :ref:`example <output-notes>`

        Console argument: ``-on`` / ``--output-name``
    """

    # Get the output path and the output name from the input filename
    # if they are not specified
    head, tail = os.path.split(input_filename)
    out_head, out_tail = os.path.split(output_name)
    if len(out_head) == 0:
        out_head = head
    if len(out_tail) == 0:
        out_tail = tail + "_identified"

    # Create the output directory if it does not exist
    if out_head != "":
        os.makedirs(out_head, exist_ok=True)

    output_name = os.path.join(out_head, out_tail)

    # Set the default separation tolerance for the span
    separation_tolerance = 10e-8

    # Read atoms and centres
    atom_counter = {}
    with open(input_filename, "r") as file:
        # Read the number of atoms and centres
        n = int(file.readline())

        # Read the file stats
        file_stats = file.readline()

        # Read the atoms and centres
        centres = []
        atoms = []
        for _ in range(0, n):
            line = file.readline().split()
            if line[0] == "X":
                centres.append([f"X", np.array(list(map(float, line[1:])))])
            else:
                if line[0] in atom_counter:
                    atom_counter[line[0]] += 1
                else:
                    atom_counter[line[0]] = 1
                atoms.append(
                    [
                        f"{line[0]}{atom_counter[line[0]]}",
                        np.array(list(map(float, line[1:]))),
                    ]
                )

    # Identify the centres by the closest atom
    for centre in centres:
        min_span = 10000
        name = "None"

        # Find the closest atom
        for atom, a_coord in atoms:
            if np.linalg.norm((centre[1] - a_coord)) < min_span:
                min_span = np.linalg.norm((centre[1] - a_coord))
                name = atom

        if min_span - span > separation_tolerance:
            cprint(
                f"Centre {centre} unidentified, " + "try to increase span",
                "yellow",
            )
            print(
                f"    span limit = {span}\n"
                + f"    minimum distance to the atom ({name}) = {min_span:.8f}\n"
            )
            centre[0] = f"None, closest: {name} ({min_span:.8f})"
        else:
            centre[0] = name

    # Write the output
    with open(output_name, "w") as file:
        file.write(f"{len(atoms) + len(centres):6.0f}\n")
        file.write(f"{file_stats}")
        for centre, coordinate in centres:
            file.write(
                f"X      "
                + f"{coordinate[0]:14.8f}   "
                + f"{coordinate[1]:14.8f}   "
                + f"{coordinate[2]:14.8f}"
                + f"   ->   {centre}\n"
            )
        for atom, coordinate in atoms:
            file.write(
                f"{atom.translate(str.maketrans('','','0123456789')):2}     "
                + f"{coordinate[0]:14.8f}   "
                + f"{coordinate[1]:14.8f}   "
                + f"{coordinate[2]:14.8f}"
                + f"   ->   {atom}\n"
            )
    cprint(f"Results are in {os.path.abspath(output_name)}", "blue")


def create_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-if",
        "--input-filename",
        required=True,
        metavar="filename",
        type=str,
        help='Relative or absolute path to the "*_centres.xyz" file',
    )
    parser.add_argument(
        "-s",
        "--span",
        default=0.1,
        type=float,
        help="Distance tolerance between centre and atom. (in Angstroms)",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        default="",
        type=str,
        help="Seedname for the output files.",
    )

    return parser
