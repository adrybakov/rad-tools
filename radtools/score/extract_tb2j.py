from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime
from os import makedirs
from os.path import abspath, split

from termcolor import cprint

from radtools import __version__ as version
from radtools.io.internal import read_template
from radtools.io.tb2j import read_tb2j_model


def manager(
    input_filename,
    template_file,
    output_name=None,
    decimals=4,
    force_symmetry=False,
    isotropic=False,
    anisotropic=False,
    matrix=False,
    dmi=False,
    all=False,
    verbose=False,
):
    r"""
    :ref:`rad-extract-tb2j` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-extract-tb2j>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    # Create the output directory if it does not exist
    makedirs(split(output_name)[0], exist_ok=True)

    # Get current date and time
    cd = datetime.now()

    if all:
        isotropic = True
        anisotropic = True
        matrix = True
        dmi = True

    # Read the model and the template
    model = read_tb2j_model(input_filename, quiet=not verbose)
    template = read_template(template_file)

    # Get the summary of the model
    summary_txt = model.summary_as_txt(
        template=template,
        decimals=decimals,
        force_symmetry=force_symmetry,
        isotropic=isotropic,
        anisotropic=anisotropic,
        matrix=matrix,
        dmi=dmi,
    )

    # Write the summary to the file
    if output_name is not None:
        with open(output_name, "w") as out_file:
            out_file.write(
                f"Exchange values are extracted from: {input_filename}\n"
                + f"on {cd.day} {month_name[cd.month]} {cd.year}"
                + f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools {version}\n\n"
            )
            out_file.write(summary_txt)
        cprint(
            f"Extracted exchange info is in " + f"{abspath(output_name)}",
            "green",
        )
    # Print the summary to the terminal
    else:
        print(f"{summary_txt}")


def create_parser():
    parser = ArgumentParser(
        description="Script for extracting of template-based model from TB2J results."
    )

    parser.add_argument(
        "-if",
        "--input-filename",
        metavar="filename",
        type=str,
        required=True,
        help="Relative or absolute path to the 'exchange.out' file, "
        + "including the name and extension of the file itself.",
    )
    parser.add_argument(
        "-tf",
        "--template-file",
        metavar="filename",
        type=str,
        required=True,
        help="Relative or absolute path to the template file, "
        + "including the name and extension of the file.",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        metavar="filename",
        type=str,
        default=None,
        help="Seedname for the output files.",
    )
    parser.add_argument(
        "-d",
        "--decimals",
        metavar="n",
        type=int,
        default=4,
        help="Decimals after the comma for the exchange values.",
    )
    parser.add_argument(
        "-fs",
        "--force-symmetry",
        action="store_true",
        default=False,
        help="Whether to force the symmetry of the template on the model.",
    )
    parser.add_argument(
        "-i",
        "--isotropic",
        action="store_true",
        default=False,
        help="Whether to output isotropic exchange.",
    )
    parser.add_argument(
        "-a",
        "--anisotropic",
        action="store_true",
        default=False,
        help="Whether to output anisotropic exchange.",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        action="store_true",
        default=False,
        help="Whether to output whole matrix exchange.",
    )
    parser.add_argument(
        "-dmi",
        action="store_true",
        default=False,
        help="Whether to output DMI exchange.",
    )
    parser.add_argument(
        "-all",
        action="store_true",
        default=False,
        help="Whether to all types of exchange.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )

    return parser
