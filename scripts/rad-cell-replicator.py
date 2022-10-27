#! /usr/local/bin/python3

import numpy as np

from rad_tools.routines import OK, RESET, YELLOW, spaces_around


def main():
    print(f"{OK}Provide lattice vector in a matrix form " +
          f"(first line - a-vector, second - b-vector, third - c-vector){RESET}")
    a_vector = np.array(tuple(map(float, input().split())))

    b_vector = np.array(tuple(map(float, input().split())))

    c_vector = np.array(tuple(map(float, input().split())))

    print(f"{OK}Please provide Atoms coordinates and names.\n" +
          f"Keyword for the end: 'end'{RESET}\n")
    goon = True
    atoms = []
    while goon:
        line = input()
        if line == "end":
            goon = False
        else:
            line = line.split()
            atom = line[0]
            coordinates = np.array(list(map(float, line[1:4])))
            atoms.append((atom, coordinates))

    print(f"{OK}Please provide supercell`s dimensions as" +
          f" three integers separated by spaces.{RESET}\n")

    na, nb, nc = tuple(map(int, input().split()))
    A_vector = a_vector * na
    B_vector = b_vector * nb
    C_vector = c_vector * nc

    supercell = []

    for k in range(0, nc):
        for j in range(0, nb):
            for i in range(0, na):
                for atom, coordinates in atoms:
                    rel_coordinates = ((coordinates[0] + i) / na,
                                       (coordinates[1] + j) / nb,
                                       (coordinates[2] + k) / nc)
                    supercell.append(f"{spaces_around(atom, nchars=10)} " +
                                     f"{spaces_around(round(rel_coordinates[0], 10), nchars=15, align='right')} " +
                                     f"{spaces_around(round(rel_coordinates[1], 10), nchars=15, align='right')} " +
                                     f"{spaces_around(round(rel_coordinates[2], 10), nchars=15, align='right')} ")

    print(f"{YELLOW}Here is a new lattice vectors:{OK}")
    print(
        f"    {spaces_around(round(A_vector[0], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(A_vector[1], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(A_vector[2], 10), nchars=15, align='right')}\n" +
        f"    {spaces_around(round(B_vector[0], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(B_vector[1], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(B_vector[2], 10), nchars=15, align='right')}\n" +
        f"    {spaces_around(round(C_vector[0], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(C_vector[1], 10), nchars=15, align='right')}" +
        f"    {spaces_around(round(C_vector[2], 10), nchars=15, align='right')}\n")
    print(f"{RESET}")

    print(f"{YELLOW}Here is an atoms list:{OK}")
    for line in supercell:
        print(line)
    print(f"{RESET}")


if __name__ == "__main__":

    main()
