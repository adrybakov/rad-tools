from argparse import ArgumentParser
from os.path import join, abspath
from os import makedirs
from datetime import datetime
from calendar import month_name


from rad_tools.io.internal import read_template
from rad_tools.io import read_tb2j_model
from rad_tools.routines import OK, RESET


def manager(filename,
            template_file,
            output_dir=".",
            output_name=None,
            accuracy=4,
            force_symmetry=False,
            isotropic=False,
            anisotropic=False,
            matrix=False,
            dmi=False,
            all=False):
    r"""
    Main function call of the tb2j-extractor.py script.

    Full documentation on the behaviour is available 
    :ref:`here <tb2j-extractor>`. Parameters of the function directly 
    corresponds to the arguments of the script.

    If you want to have the behaviour of the tb2j-extractor.py script 
    but in a format of a fuction call use this function.
    """

    try:
        makedirs(output_dir)
    except FileExistsError:
        pass

    cd = datetime.now()

    if all:
        isotropic = True
        anisotropic = True
        matrix = True
        dmi = True

    model = read_tb2j_model(filename, quiet=True)
    template = read_template(template_file)
    summary_txt = model.summary_as_txt(template=template,
                                       decimals=accuracy,
                                       force_symmetry=force_symmetry,
                                       isotropic=isotropic,
                                       anisotropic=anisotropic,
                                       out_matrix=matrix,
                                       out_dmi=dmi)

    if output_name is not None:
        with open(join(output_dir, output_name + ".txt"), "w") as out_file:
            out_file.write(
                f"Exchange values are extracted from: {filename}\n" +
                f"on {cd.day} {month_name[cd.month]} {cd.year}" +
                f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools\n\n")
            out_file.write(summary_txt)
        print(f"{OK}Extracted exchange info is in " +
              f"{abspath(join(output_dir, output_name + '.txt'))}{RESET}")
    else:
        print(f"{summary_txt}")


def get_parser():
    parser = ArgumentParser(
        description="Script for extracting of template-based model from TB2J results.")

    parser.add_argument("-f", "--filename",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the exchange.out file,
                        including the name and extention of the file itself.
                        """)
    parser.add_argument("-tf", "--template-file",
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
    parser.add_argument("-acc", "--accuracy",
                        type=int,
                        default=4,
                        help="""
                        Accuracy for the exchange values.
                        """
                        )
    parser.add_argument("-fs", "--force-symmetry",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to force the symmetry of the
                        template on the model.
                        """)
    parser.add_argument("-i", "--isotropic",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to output isotropic exchange.
                        """)
    parser.add_argument("-a", "--anisotropic",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to output anisotropic exchange.
                        """)
    parser.add_argument("-m", "--matrix",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to output whole matrix exchange.
                        """)
    parser.add_argument("-dmi",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to output DMI exchange.
                        """)
    parser.add_argument("-all",
                        action="store_true",
                        default=False,
                        help="""
                        Whenever to all types of exchange.
                        """)

    return parser
