#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join, abspath
from os import makedirs
from math import sqrt

import numpy as np

from rad_tools.io.internal import read_template
from rad_tools.io.tb2j import read_exchange_model
from rad_tools.routines import OK, RESET


def extract(filename, out_dir, out_name, template, dmi_verbose=False, verbose=False):
    if verbose and dmi_verbose:
        dmi_verbose = False

    model = read_exchange_model(filename)
    template = read_template(template)
    summary_txt = model.summary_as_txt(template=template,
                                       dmi_verbose=dmi_verbose, verbose=verbose)
    if verbose:
        summary_py = model.summary_as_py(template=template)

    if out_name is not None:
        with open(join(out_dir, out_name + ".txt"), "w") as out_file:
            out_file.write(summary_txt)
        print(f"{OK}Extracted exchange info is in " +
              f"{abspath(join(out_dir, out_name + '.txt'))}{RESET}")
        if verbose:
            with open(join(out_dir, out_name + ".py"), "w") as out_python:
                out_python.write(summary_py)
            print(f"{OK}Verbose info in .py format is in " +
                  f"{abspath(join(out_dir, out_name + '.py'))}{RESET}")
    else:
        print(f"{summary_txt}")


if __name__ == '__main__':
    parser = ArgumentParser(
        description="Script for extracting of template-based model from TB2J results.",
        epilog="""
            For the full description of arguments see the docs:
            https://rad-tools.adrybakov.com/en/stable/user-guide/tb2j-extractor.html
            """)

    parser.add_argument("-f", "--filename",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the * exchange.out * file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-tf", "--template",
                        type=str,
                        required=True,
                        help="""
                        Relative or absolute path to the template file,
                        including the name and extention of the file.
                        """)
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default=".",
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
    parser.add_argument("-dmi",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to print each dmi vector
                        for each exchange group separately.
                        """
                        )
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to print each neighbor in the template
                        in a verbose way.
                        """
                        )

    args = parser.parse_args()

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

    extract(filename=args.filename,
            out_dir=args.output_dir,
            out_name=args.output_name,
            template=args.template,
            dmi=args.dmi,
            verbose=args.verbose)
