#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join, split

import numpy as np

from rad_tools.routines import cprint


def manager(filename, span, output_path, output_name, no_colour=False):
    head, tail = split(filename)
    if output_path is None:
        output_path = head
    if output_name is None:
        output_name = tail + "_identified"
    separation_tolerance = 10e-8

    # Read atoms and centres
    atom_counter = {}
    with open(filename, "r") as file:
        n = int(file.readline())
        file_stats = file.readline()
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

    # Identify centres localization
    for centre in centres:
        min_span = 10000
        name = "None"
        for atom, a_coord in atoms:
            if np.linalg.norm((centre[1] - a_coord)) < min_span:
                min_span = np.linalg.norm((centre[1] - a_coord))
                name = atom
        if min_span - span > separation_tolerance:
            if no_colour:
                print(f"Centre {centre} unidentified, " + "try to increase span")
            else:
                cprint(
                    f"Centre {centre} unidentified, " + "try to increase span",
                    colour="yellow",
                )
            print(
                f"    span limit = {span}\n"
                + f"    minimum distance to the atom ({name}) = {min_span:.8f}\n"
            )
            centre[0] = f"None, closest: {name} ({min_span:.8f})"
        else:
            centre[0] = name

    # Write the output
    with open(join(output_path, output_name), "w") as file:
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


def create_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "filename",
        type=str,
        help="Relative or absolute path to the _centres.xyz file",
    )
    parser.add_argument(
        "-s",
        "--span",
        type=float,
        default=0.1,
        help="Distance tolerance between centre and atom. (in Angstroms)",
    )
    parser.add_argument(
        "-op",
        "--output-path",
        type=str,
        default=None,
        help="Relative or absolute path to the folder for saving outputs.",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        type=str,
        default=None,
        help="Seedname for the output files.",
    )
    parser.add_argument(
        "-nc",
        "--no-colour",
        action="store_true",
        default=False,
        help="Turn off coloured output.",
    )
    return parser
