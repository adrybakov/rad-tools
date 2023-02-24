from argparse import ArgumentParser
from os.path import abspath, split
from os import makedirs
from datetime import datetime
from calendar import month_name

import numpy as np

from rad_tools.io.tb2j import read_exchange_model
from rad_tools.routines import OK, RESET, YELLOW


def manager(output_name="template",
            tb2j_filename=None,
            R_vector=None,
            max_distance=None,
            min_distance=None,
            distance=None):
    r"""
    Main function call of the rad-make-template.py script.

    Full documentation on the behaviour is available
    :ref:`here <rad-make-template>`. Parameters of the function directly
    corresponds to the arguments of the script.

    If you want to have the behaviour of the rad-make-template.py script
    but in a format of a fuction call use this function.
    """

    if distance is not None:
        min_distance = distance
        max_distance = distance

    if R_vector is not None:
        R_vector = np.array(R_vector[:len(R_vector) // 3 * 3],
                            dtype=int).reshape((len(R_vector)//3, 3))
        R_vector = list(map(tuple, R_vector.tolist()))

    try:
        makedirs(split(output_name)[0])
    except FileExistsError:
        pass
    except FileNotFoundError:
        pass

    template = ('=' * 20 + "\n" +
                "Neighbors template:\n" +
                "i j R_a R_b R_c\n" +
                '-' * 20 + "\n" +
                "J1 $J_1$\n" +
                "atom1 atom2  0  0  0\n" +
                "atom1 atom2  1  0  0\n" +
                "atom1 atom1 -1  0  2\n" +
                '-' * 20 + "\n" +
                "J2\n" +
                "atom2 atom1  9  5 -3\n" +
                "atom1 atom2  1  4  0\n" +
                "atom2 atom2  1  0  2\n" +
                '=' * 20 + "\n")
    cd = datetime.now()
    with open(f"{output_name}.txt", "w") as file:
        if tb2j_filename is None:
            file.write(
                f"Template is created " +
                f"on {cd.day} {month_name[cd.month]} {cd.year}" +
                f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools\n\n")

            file.write(template)
        else:
            file.write(
                f"Template is created based on the file: {tb2j_filename}\n" +
                f"on {cd.day} {month_name[cd.month]} {cd.year}" +
                f" at {cd.hour}:{cd.minute}:{cd.second} by rad-tools\n\n")

            model = read_exchange_model(tb2j_filename)
            model.filter(min_distance=min_distance,
                         max_distance=max_distance,
                         R_vector=R_vector)
            file.write('=' * 20 + "\n" +
                       "Neighbors template:\n" +
                       "i j R_a R_b R_c\n" +
                       '-' * 20 + "\n" + "Name placeholder" + "\n")
            for atom1, atom2, R in model.bonds:
                file.write(f"{atom1:4} {atom2:4} " +
                           f"{R[0]:3.0f} {R[1]:3.0f} {R[2]:3.0f}\n")
            file.write('=' * 20 + "\n")
    print(f"{OK}Template draft is in " +
          f"{abspath(output_name)}.txt{RESET}")
    print(f"{YELLOW}Do not forget to correct the template draft " +
          f"to your needs!{RESET}")


def get_parser():
    parser = ArgumentParser(
        description="Script for the creation of template`s template")

    parser.add_argument("-on", "--output-name",
                        type=str,
                        default='template',
                        help="""
                        Relative or absolute path to the template output file.
                        """
                        )
    parser.add_argument("-tf", "--tb2j-filename",
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

    return parser
