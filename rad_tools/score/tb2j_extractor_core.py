from argparse import ArgumentParser
from calendar import month_name
from datetime import datetime
from os import makedirs
from os.path import abspath, join

from rad_tools import __version__ as version
from rad_tools.io import read_tb2j_model
from rad_tools.io.internal import read_template
from rad_tools.routines import OK, RESET


def manager(
    input_filename,
    template_file,
    output_path=".",
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
    Main function call of the tb2j-extractor.py script.

    Full documentation on the behaviour is available
    :ref:`here <tb2j-extractor>`. Parameters of the function directly
    corresponds to the arguments of the script.

    If you want to have the behaviour of the tb2j-extractor.py script
    but in a format of a function call use this function.
    """

    try:
        makedirs(output_path)
    except FileExistsError:
        pass

    cd = datetime.now()

    if all:
        isotropic = True
        anisotropic = True
        matrix = True
        dmi = True

    model = read_tb2j_model(input_filename, quiet=not verbose)
    template = read_template(template_file)
    summary_txt = model.summary_as_txt(
        template=template,
        decimals=decimals,
        force_symmetry=force_symmetry,
        isotropic=isotropic,
        anisotropic=anisotropic,
        out_matrix=matrix,
        out_dmi=dmi,
    )

    if output_name is not None:
        with open(join(output_path, output_name + ".txt"), "w") as out_file:
            out_file.write(
                f"Exchange values are extracted from: {input_filename}\n"
                + f"on {cd.day} {month_name[cd.month]} {cd.year}"
                + f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools {version}\n\n"
            )
            out_file.write(summary_txt)
        print(
            f"{OK}Extracted exchange info is in "
            + f"{abspath(join(output_path, output_name + '.txt'))}{RESET}"
        )
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
        "-op",
        "--output-path",
        metavar="path",
        type=str,
        default=".",
        help="Relative or absolute path to the folder for saving outputs.",
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
        help="Whenever to force the symmetry of the template on the model.",
    )
    parser.add_argument(
        "-i",
        "--isotropic",
        action="store_true",
        default=False,
        help="Whenever to output isotropic exchange.",
    )
    parser.add_argument(
        "-a",
        "--anisotropic",
        action="store_true",
        default=False,
        help="Whenever to output anisotropic exchange.",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        action="store_true",
        default=False,
        help="Whenever to output whole matrix exchange.",
    )
    parser.add_argument(
        "-dmi",
        action="store_true",
        default=False,
        help="Whenever to output DMI exchange.",
    )
    parser.add_argument(
        "-all",
        action="store_true",
        default=False,
        help="Whenever to all types of exchange.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )

    return parser
