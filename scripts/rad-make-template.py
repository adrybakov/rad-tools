#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join
from os import makedirs

import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModelTB2J


def provide_template(output_dir=".", output_name="template.txt",
                     tb2j_filename=None,
                     min_distance=None,
                     max_distance=None,
                     R_vector=None):
    template = '=' * 20 + "\n" + \
        "Neighbors template:\n"\
        "i j R_a R_b R_c\n" +\
        '-' * 20 + "\n" + \
        "J1 $J_1$\n"\
        "atom1 atom2 0 0 0\n"\
        "atom1 atom2 1 0 0\n"\
        "atom1 atom1 -1 0 2\n" +\
        '-' * 20 + "\n" + \
        "J2\n"\
        "atom2 atom1 9 5 -3\n"\
        "atom1 atom2 1 4 0\n"\
        "atom2 atom2 1 0 2\n" +\
        '=' * 20 + "\n"

    with open(join(output_dir, output_name), "w") as file:
        if tb2j_filename is None:
            file.write(template)
        else:
            model = ExchangeModelTB2J(tb2j_filename)
            model = model.filter(min_distance=min_distance,
                                 max_distance=max_distance,
                                 R_vector=R_vector)
            file .write('=' * 20 + "\n" +
                        "Neighbors template:\n" +
                        "i j R_a R_b R_c\n" +
                        '-' * 20 + "\n")
            for atom1, atom2, R in model.file_order:
                file.write(f"{atom1} {atom2} {R[0]} {R[1]} {R[2]}\n")
            file .write('=' * 20 + "\n")


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script for the creation of template`s template")

    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder for saving outputs.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='template.txt',
                        help="""
                        Template file name, default "template.txt"
                        """
                        )
    parser.add_argument("-f", "--filename",
                        type=str,
                        default=None,
                        help="""
                        Relative or absulute path to the *exchange.out* file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-R", "--R-vector",
                        type=int,
                        nargs="*",
                        default=None,
                        help="""
                        R vectors for filtering the model.
                        """
                        )
    parser.add_argument("-maxd", "--max-distance",
                        type=float,
                        default=None,
                        help="""
                        (<=) Maximum distance.
                        """
                        )
    parser.add_argument("-mind", "--min-distance",
                        type=float,
                        default=None,
                        help="""
                        (>=) Minimum distance.
                        """
                        )
    parser.add_argument("-d", "--distance",
                        type=float,
                        default=None,
                        help="""
                        (=) Exact distance.
                        """
                        )

    args = parser.parse_args()

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

    if args.distance is not None:
        args.min_distance = args.max_distance = args.distances

    if args.R_vector is not None:
        args.R_vector = np.array(args.R_vector[:len(args.R_vector) // 3 * 3],
                                 dtype=int).reshape((len(args.R_vector)//3, 3))
        args.R_vector = list(map(tuple, args.R_vector.tolist()))

    provide_template(output_dir=args.output_dir,
                     output_name=args.output_name,
                     tb2j_filename=args.filename,
                     min_distance=args.min_distance,
                     max_distance=args.max_distance,
                     R_vector=args.R_vector)
