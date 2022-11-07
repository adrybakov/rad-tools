#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join, abspath
from os import makedirs
from math import sqrt

import numpy as np

from rad_tools.exchange.model import ExchangeModelTB2J
from rad_tools.exchange.template import ExchangeTemplate
from rad_tools.routines import OK, RESET, spaces_around


def main(filename, out_dir, out_name, template, dmi=False):

    model = ExchangeModelTB2J(filename)
    template = ExchangeTemplate(template)

    output_line = ""
    for name in template.names:
        J_iso = 0
        J_aniso = np.zeros((3, 3), dtype=float)
        DMI = np.zeros(3, dtype=float)
        abs_DMI = 0
        for bond in template.names[name]:
            atom1 = bond[0]
            atom2 = bond[1]
            R = bond[2]
            J_iso += model.bonds[atom1][atom2][R].iso
            J_aniso += model.bonds[atom1][atom2][R].aniso
            DMI += model.bonds[atom1][atom2][R].dmi
            abs_DMI += sqrt(model.bonds[atom1][atom2][R].dmi[0]**2 +
                            model.bonds[atom1][atom2][R].dmi[1]**2 +
                            model.bonds[atom1][atom2][R].dmi[2]**2)
        J_iso /= len(template.names[name])
        J_aniso /= len(template.names[name])
        DMI /= len(template.names[name])
        abs_DMI /= len(template.names[name])
        output_line += (
            f"{name}\n" +
            f"    Isotropic: {round(J_iso, 4)}\n" +
            f"    Anisotropic:\n" +
            f"        {spaces_around(round(J_aniso[0][0], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[0][1], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[0][2], 4), nchars=7)}\n" +
            f"        {spaces_around(round(J_aniso[1][0], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[1][1], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[1][2], 4), nchars=7)}\n" +
            f"        {spaces_around(round(J_aniso[2][0], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[2][1], 4), nchars=7)}  " +
            f"{spaces_around(round(J_aniso[2][2], 4), nchars=7)}\n" +
            f"    |DMI|: {round(abs_DMI, 4)}\n" +
            f"    |DMI/J| {round(abs(abs_DMI/J_iso), 4)}\n")

        if dmi:
            for bond in template.names[name]:
                atom1 = bond[0]
                atom2 = bond[1]
                R = bond[2]
                DMI = model.bonds[atom1][atom2][R].dmi
                output_line += (
                    f"    DMI: " +
                    f"{spaces_around(round(DMI[0], 4), nchars=7)} " +
                    f"{spaces_around(round(DMI[1], 4), nchars=7)} " +
                    f"{spaces_around(round(DMI[2], 4), nchars=7)} {R}\n")
            output_line += "\n"
        else:
            output_line += (
                f"    DMI: " +
                f"{round(DMI[0], 4)} " +
                f"{round(DMI[1], 4)} " +
                f"{round(DMI[2], 4)}\n\n")

    if out_name is not None:
        with open(join(out_dir, out_name), "w") as out_file:
            out_file.write(output_line)
        print(f"{OK}Extracted exchange info is in " +
              f"{abspath(join(out_dir, out_name))}{RESET}")
    else:
        print(f"{OK}{output_line}{RESET}")


if __name__ == '__main__':
    parser = ArgumentParser(
        description="Script for extracting of template-based model from TB2J results.",
        epilog="""
            For the full description of arguments see the docs:
            https: // rad-tools.adrybakov.com/en/stable/user-guide/tb2j-refractor.html
            """)

    parser.add_argument("-f", "--filename",
                        type=str,
                        required=True,
                        help="""
                        Relative or absulute path to the * exchange.out * file,
                        including the name and extention of the file itself.
                        """
                        )
    parser.add_argument("-t", "--template",
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

    args = parser.parse_args()

    try:
        makedirs(args.output_dir)
    except FileExistsError:
        pass

    main(filename=args.filename,
         out_dir=args.output_dir,
         out_name=args.output_name,
         template=args.template,
         dmi=args.dmi)
