#! /usr/local/bin/python3

from argparse import ArgumentParser
from math import sqrt

import numpy as np

from rad_tools.routines import strip_digits, WARNING, RESET


def search_on_atoms(centre, atoms, span):
    name = "None"
    min_span = 100
    atom = ""
    for atom, a_coord in atoms:
        if sqrt(np.sum((centre - a_coord)**2)) < span:
            name = atom
        if sqrt(np.sum((centre - a_coord)**2)) < min_span:
            min_span = sqrt(np.sum((centre - a_coord)**2))
            atom = atom
    return name, min_span, atom


def search_between_atoms(centre, atoms, span):
    pairs = []
    for i, atom in enumerate(atoms):
        for j in range(i+1, len(atoms)):
            pair = f"{atom[0]}-{atoms[j][0]}"
            p_coord = (atom[1]+atoms[j][1])/2
            pairs.append((pair, p_coord))
    return search_on_atoms(centre, pairs, span)


def identify(filename, span=0.1):

    # Read atoms and centres
    atom_counter = {}
    with open(filename, "r") as file:
        n = int(file.readline())
        file_stats = file.readline()
        centres = []
        atoms = []
        for i in range(0, n):
            line = file.readline().split()
            if line[0] == "X":
                centres.append(np.array(list(map(float, line[1:]))))
            else:
                if line[0] in atom_counter:
                    atom_counter[line[0]] += 1
                else:
                    atom_counter[line[0]] = 1
                atoms.append((f"{line[0]}{atom_counter[line[0]]}",
                              np.array(list(map(float, line[1:])))))

    # Identify centres localization
    centres_names = []
    for centre in centres:
        name_atom, min_span_atom, atom = search_on_atoms(
            centre, atoms, span)
        name_pair, min_span_pair, pair = search_between_atoms(
            centre, atoms, span)
        if name_atom == "None" and name_pair == "None":
            # print(f"{WARNING}" +
            #       f"Centre {centre} unindentified, " +
            #       "try to increase --span\n" +
            #       f"{RESET}" +
            #       f"    span limit = {span}\n" +
            #       f"    centre`s min span = {min_span_atom:.8f} " +
            #       f"(with {atom} atom)\n"
            #       f"    centre`s min span = {min_span_pair:.8f} " +
            #       f"(with centre point between {pair} atoms)\n")
            pass
        print(name_atom, min_span_atom, atom,
              name_pair, min_span_pair, pair)
        if min_span_atom < min_span_pair:
            centres_names.append(name_atom)
        else:
            centres_names.append(name_pair)

    # Write the output
    with open(f"{filename}_identified", "w") as file:
        file.write(f"{len(atoms) + len(centres):6.0f}\n")
        file.write(f"{file_stats}")
        for c_i, center in enumerate(centres):
            file.write(f"X      " +
                       f"{center[0]:14.8f}   " +
                       f"{center[1]:14.8f}   " +
                       f"{center[2]:14.8f}" +
                       f"   ->   {centres_names[c_i]:4}\n")
        for atom, coord in atoms:
            file.write(f"{strip_digits(atom):2}     " +
                       f"{coord[0]:14.8f}   " +
                       f"{coord[1]:14.8f}   " +
                       f"{coord[2]:14.8f}" +
                       f"   ->   {atom:4}\n")


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
    args = parser.parse_args()

    identify(args.filename, args.span)
