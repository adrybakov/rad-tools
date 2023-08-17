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
    template_file=None,
    output_name=None,
    decimals=4,
    form_model=False,
    no_anisotropic=False,
    no_matrix=False,
    nodmi=False,
    verbose=False,
    max_distance=None,
    min_distance=None,
):
    r"""
    :ref:`rad-extract-tb2j` script.

    Full documentation on the behaviour is available in the
    :ref:`User Guide <rad-extract-tb2j>`.
    Parameters of the function directly
    correspond to the arguments of the script.
    """

    if form_model and template_file is None:
        raise ValueError(
            "Template file is required for forming the model from the template."
            + "--form-model implies --template-file"
        )

    # Create the output directory if it does not exist
    if output_name is not None and split(output_name)[0] != "":
        makedirs(split(output_name)[0], exist_ok=True)

    # Get current date and time
    cd = datetime.now()

    # Read the model and the template
    model = read_tb2j_model(input_filename, quiet=not verbose)
    if template_file is not None:
        template = read_template(template_file)
    else:
        template = None

    model.filter(
        max_distance=max_distance,
        min_distance=min_distance,
        template=template,
    )

    # Remove template if no model formation is required
    if not form_model:
        template = None

    # Write the summary to the file or print it
    model.dump_txt(
        filename=output_name,
        anisotropic=not no_anisotropic,
        matrix=not no_matrix,
        dmi=not nodmi,
        template=template,
        decimals=decimals,
        additional_stats=f"Exchange values are extracted from:\n"
        + f"{abspath(input_filename)}\n",
    )
    if output_name is not None:
        cprint(
            f"Extracted exchange info is in " + f"{abspath(output_name)}",
            "green",
        )


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
        "-fm",
        "--form-model",
        action="store_true",
        default=False,
        help="Whether to form the model from the template.",
    )
    parser.add_argument(
        "-noa",
        "--no-anisotropic",
        action="store_true",
        default=False,
        help="Whether to output anisotropic exchange.",
    )
    parser.add_argument(
        "-nom",
        "--no-matrix",
        action="store_true",
        default=False,
        help="Whether to output whole matrix exchange.",
    )
    parser.add_argument(
        "-nodmi",
        action="store_true",
        default=False,
        help="Whether to output DMI exchange.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
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

    return parser
