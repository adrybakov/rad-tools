# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
from argparse import ArgumentParser

from termcolor import cprint

from radtools._osfix import _winwait
from radtools.io.internal import dump_spinham_txt, load_template
from radtools.io.tb2j import load_tb2j_model


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

    Parameters
    ----------
    input_filename : str
        Relative or absolute path to the "exchange.out" file,
        including the name and extension of the file itself.

        Console argument: ``-if`` / ``--input-filename``
    template_file : str, optional
        Relative or absolute path to the template file,
        including the name and extension of the file.

        Console argument: ``-tf`` / ``--template-file``
    output_name : str, optional
        Name of the output files.

        If this parameter is not specified, the result is printed in the
        standard output stream.

        Console argument: ``-on`` / ``--output-name``
    decimals : int, default 4
        Decimals after the comma for the exchange values.

        Console argument: ``-d`` / ``--decimals``

        .. versionchanged:: 0.5.17 Renamed from ``-acc``/``--accuracy``.
    form_model : bool, default False
        Whether to form the model from the template.

        Console argument: ``-fm`` / ``--form-model``

        .. versionchanged:: 0.8.0 Renamed from ``-fs``/``--force-symmetry``.
    no_anisotropic : bool, default False
        Whether to output anisotropic exchange.

        Console argument: ``-noa`` / ``--no-anisotropic``

        .. versionchanged:: 0.8.0 Renamed from ``-a``/``--anisotropic``.
    no_matrix : bool, default False
        Whether to output whole matrix exchange.

        Console argument: ``-nom`` / ``--no-matrix``

        .. versionchanged:: 0.8.0 Renamed from ``-m``/``--matrix``.
    nodmi : bool, default False
        Whether to output DMI exchange.

        Console argument: ``-nodmi``

        .. versionchanged:: 0.8.0 Renamed from ``-dmi``.
    verbose : bool, default False
        Verbose output, propagates to the called methods.

        Console argument: ``-v`` / ``--verbose``
    max_distance : float, optional
        (<=) Maximum distance.

        All the bonds with the distance between atom 1 and atom 2
        greater than maximum distance are excluded from the model.

        Console argument: ``-maxd`` / ``--max-distance``

        .. versionadded:: 0.8.0
    min_distance : float, optional
        (>=) Minimum distance.

        All the bonds with the distance between atom 1 and atom 2
        lower than minimum distance are excluded from the Hamiltonian.

        Console argument: ``-mind`` / ``--min-distance``

        .. versionadded:: 0.8.0

    """

    if form_model and template_file is None:
        raise ValueError(
            "Template file is required for forming the model from the template."
            + "--form-model implies --template-file"
        )

    head, _ = os.path.split(input_filename)
    if output_name is not None:
        out_head, out_tail = os.path.split(output_name)
        if len(out_head) == 0:
            out_head = head
        if len(out_tail) == 0:
            out_tail = "spinham.txt"

        output_name = os.path.join(out_head, out_tail)

        # Create the output directory if it does not exist
        if out_head != "":
            os.makedirs(out_head, exist_ok=True)

    # Read the model and the template
    model = load_tb2j_model(input_filename, quiet=not verbose)
    if template_file is not None:
        template = load_template(template_file)
    else:
        template = None

    model.filter(
        max_distance=max_distance,
        min_distance=min_distance,
        template=template,
    )

    additional_stats = (
        f"Exchange values are extracted from:\n"
        + f"{os.path.abspath(input_filename)}\n"
    )
    if template_file is not None:
        additional_stats += (
            f"Template used for "
            + ("filtering" if not form_model else "model formation")
            + ":\n"
            + f"{os.path.abspath(template_file)}\n"
        )
    # Remove template if no model formation is required
    if not form_model:
        template = None

    # Write the summary to the file or print it
    dump_spinham_txt(
        model,
        filename=output_name,
        anisotropic=not no_anisotropic,
        matrix=not no_matrix,
        dmi=not nodmi,
        template=template,
        decimals=decimals,
        additional_stats=additional_stats,
    )
    if output_name is not None:
        cprint(
            f"Extracted exchange info is in " + f"{os.path.abspath(output_name)}",
            "green",
        )


def create_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-if",
        "--input-filename",
        required=True,
        type=str,
        help='Relative or absolute path to the "exchange.out" file,',
    )
    parser.add_argument(
        "-tf",
        "--template-file",
        default=None,
        type=str,
        help="Relative or absolute path to the template file,",
    )
    parser.add_argument(
        "-on",
        "--output-name",
        default=None,
        type=str,
        help="Name of the output files.",
    )
    parser.add_argument(
        "-d",
        "--decimals",
        default=4,
        type=int,
        help="Decimals after the comma for the exchange values.",
    )
    parser.add_argument(
        "-fm",
        "--form-model",
        default=False,
        action="store_true",
        help="Whether to form the model from the template.",
    )
    parser.add_argument(
        "-noa",
        "--no-anisotropic",
        default=False,
        action="store_true",
        help="Whether to output anisotropic exchange.",
    )
    parser.add_argument(
        "-nom",
        "--no-matrix",
        default=False,
        action="store_true",
        help="Whether to output whole matrix exchange.",
    )
    parser.add_argument(
        "-nodmi",
        default=False,
        action="store_true",
        help="Whether to output DMI exchange.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Verbose output, propagates to the called methods.",
    )
    parser.add_argument(
        "-maxd",
        "--max-distance",
        default=None,
        type=float,
        help="(<=) Maximum distance.",
    )
    parser.add_argument(
        "-mind",
        "--min-distance",
        default=None,
        type=float,
        help="(>=) Minimum distance.",
    )

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
    _winwait()


if __name__ == "__main__":
    main()
