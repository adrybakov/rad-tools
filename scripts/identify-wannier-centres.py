#! /usr/local/bin/python3

from argparse import ArgumentParser
from math import sqrt
from os.path import split, join

import numpy as np

from rad_tools.routines import strip_digits, WARNING, RESET, \
    search_on_atoms, search_between_atoms


def identify(filename, span, out_dir, out_name, nocolor=False):
    separation_tolerance = 10E-5

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
                centres.append([f"X",
                                np.array(list(map(float, line[1:])))])
            else:
                if line[0] in atom_counter:
                    atom_counter[line[0]] += 1
                else:
                    atom_counter[line[0]] = 1
                atoms.append([f"{line[0]}{atom_counter[line[0]]}",
                              np.array(list(map(float, line[1:])))])

    # Identify centres localization
    for centre in centres:
        min_span_atom, atom = search_on_atoms(
            centre, atoms)
        min_span_pair, pair = search_between_atoms(
            centre, atoms)
        if min_span_atom > span and min_span_pair > span:
            if not nocolor:
                print(f"{WARNING}", end="")
            print(f"Centre {centre} unindentified, " +
                  "try to increase span")
            if not nocolor:
                print(f"{RESET}", end="")
            print(f"    span limit = {span}\n" +
                  f"    minimum distance to the atom = {min_span_atom:.8f} " +
                  f"({atom})\n"
                  f"    minimum distance to the bond`s centre = {min_span_pair:.8f} " +
                  f"({pair})\n")
            if abs(min_span_atom - min_span_pair) < separation_tolerance:
                centre[0] = f"None, closest: {atom}"
            elif min_span_atom < min_span_pair:
                centre[0] = f"None, closest: {atom}"
            else:
                centre[0] = f"None, closest: {pair}"
        elif abs(min_span_atom - min_span_pair) < separation_tolerance:
            centre[0] = f"{atom} or {pair}"
        elif min_span_atom < min_span_pair:
            centre[0] = atom
        else:
            centre[0] = pair

    # Write the output
    with open(join(out_dir, out_name), "w") as file:
        file.write(f"{len(atoms) + len(centres):6.0f}\n")
        file.write(f"{file_stats}")
        for centre, coordinate in centres:
            file.write(f"X      " +
                       f"{coordinate[0]:14.8f}   " +
                       f"{coordinate[1]:14.8f}   " +
                       f"{coordinate[2]:14.8f}" +
                       f"   ->   {centre}\n")
        for atom, coordinate in atoms:
            file.write(f"{strip_digits(atom):2}     " +
                       f"{coordinate[0]:14.8f}   " +
                       f"{coordinate[1]:14.8f}   " +
                       f"{coordinate[2]:14.8f}" +
                       f"   ->   {atom}\n")


if __name__ == "__main__":
    parser = ArgumentParser(
        description=("Identify wannier centres with respect to the atom " +
                     "or to the point between the atom`s pair"))
    parser.add_argument("filename",
                        type=str,
                        help="""
                        Rellative or absolute path to the _centres.xyz file
                        """)
    parser.add_argument("-s", "--span",
                        type=float,
                        default=0.1,
                        help="""
                        Distance tolerance between centre and atom. (in Angstrom)
                        """)
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default=None,
                        help="""
                        Relative or absolute path to the folder
                        for saving outputs.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default=None,
                        help="""
                        Seedname for the output files.
                        """
                        )
    parser.add_argument("-nc", "--no-colour",
                        action="store_true",
                        default=False,
                        help="""
                        Turn off coloured output.
                        """
                        )
    args = parser.parse_args()

    head, tail = split(args.filename)
    if args.output_dir is None:
        args.output_dir = head
    if args.output_name is None:
        args.output_name = tail + "_identified"

    identify(args.filename, args.span,
             out_dir=args.output_dir, out_name=args.output_name,
             nocolor=args.no_colour)
