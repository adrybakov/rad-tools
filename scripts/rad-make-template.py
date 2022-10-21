#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join
from os import makedirs

def provide_template(output_dir=".", output_name="template.txt"):
    template = '=' * 20 + "\n" + \
        "Neighbors template:\n"\
        "i j R_a R_b R_c\n" +\
        '-' * 20 + "\n" + \
        "J1 $J_1$\n"\
        "atom1 atom2 0 0 0\n"\
        "atom1 atom2 1 0 0\n"\
        "atom1 atom1 -1 0 2\n" +\
        '-' * 20 + "\n" + \
        "J2 $J_2$\n"\
        "atom2 atom1 9 5 -3\n"\
        "atom1 atom2 1 4 0\n"\
        "atom2 atom2 1 0 2\n" +\
        '=' * 20 + "\n"

    with open(join(output_dir, output_name), "w") as file:
        file.write(template)
    

if __name__ == "__main__":
    parser = ArgumentParser(description="Script for template creation")

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

    args = parser.parse_args()

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

    provide_template(output_dir=args.output_dir,
         output_name=args.output_name)
