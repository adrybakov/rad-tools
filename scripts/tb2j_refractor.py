#! /usr/local/bin/python3
from argparse import ArgumentParser
from os.path import join, split, abspath
from math import atan, sqrt

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from rad_tools.tb2j_tools.file_logic import ExchangeModelTB2J
from rad_tools.routines import check_make_dir, OK, RESET


def main(filename, out_dir, out_name, template):

    model = ExchangeModelTB2J(filename)
    model.filter(template=template)
    with open(join(out_dir, out_name)) as out_file:
        out_file.write()

    # model = model.filter(min_distance=min_distance,
    #                      max_distance=max_distance,
    #                      R_vector=R_vector,
    #                      template=template)


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for refractoring of TB2J results",
                            epilog="""
                            #TODO
                            """)

    parser.add_argument("-f", "--file",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the *exchange.out* file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="""
                        Relative or absolute path to the folder for saving outputs.

    If the folder does not exist then it is created from the specified path.
    The creation is applied recursevly to the path, starting from the right
    until the existing folder is reached.
                        """
                        )
    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='exchange',
                        help="""
                        Seedname for the output files.

    Output files will have the following name structure:
    output-name.display_data_type.png
                        """
                        )
    parser.add_argument("-t", "--template",
                        type=str,
                        default=None,
                        help="""
                        Relative or absolute path to the template file, 
                        including the name and extention of the file.
                        """)

    args = parser.parse_args()

    check_make_dir(args.output_dir)

    main(filename=args.file,
         out_dir=args.output_dir,
         out_name=args.output_name,
         template=args.template)
