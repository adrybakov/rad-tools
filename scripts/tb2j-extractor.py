#! /usr/local/bin/python3

from argparse import ArgumentParser
from os.path import join, abspath
from os import makedirs
from math import sqrt

import numpy as np

from rad_tools.exchange.model import ExchangeModelTB2J
from rad_tools.exchange.template import ExchangeTemplate
from rad_tools.routines import OK, RESET


def main(filename, out_dir, out_name, template, dmi=False, verbose=False):
    if verbose and dmi:
        dmi = False

    model = ExchangeModelTB2J(filename)
    template = ExchangeTemplate(template)

    output_line = ""
    if verbose:
        output_python_iso = "iso = {\n"
        output_python_aniso = "aniso = {\n"
        output_python_dmi = "dmi = {\n"
        output_python_matrix = "matrix = {\n"
    for name in template.names:

        output_line += (f"{name}\n")
        if verbose:
            output_python_iso += f"    '{name}':\n" + "    {\n"
            output_python_aniso += f"    '{name}':\n" + "    {\n"
            output_python_dmi += f"    '{name}':\n" + "    {\n"
            output_python_matrix += f"    '{name}':\n" + "    {\n"

            for bond in template.names[name]:
                atom1 = bond[0]
                atom2 = bond[1]
                R = bond[2]
                J_iso = model.bonds[atom1][atom2][R].iso
                J_aniso = model.bonds[atom1][atom2][R].aniso
                DMI = model.bonds[atom1][atom2][R].dmi
                abs_DMI = sqrt(model.bonds[atom1][atom2][R].dmi[0]**2 +
                               model.bonds[atom1][atom2][R].dmi[1]**2 +
                               model.bonds[atom1][atom2][R].dmi[2]**2)
                matrix = model.bonds[atom1][atom2][R].matrix
                output_python_iso += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): {J_iso},\n")
                output_python_aniso += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): np.array([" +
                    f"[{J_aniso[0][0]}, {J_aniso[0][1]}, {J_aniso[0][2]}], " +
                    f"[{J_aniso[1][0]}, {J_aniso[1][1]}, {J_aniso[1][2]}], " +
                    f"[{J_aniso[2][0]}, {J_aniso[2][1]}, {J_aniso[2][2]}]]),\n")
                output_python_dmi += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}):" +
                    f" np.array([{DMI[0]}, {DMI[1]}, {DMI[2]}]),\n")
                output_python_matrix += (
                    8 * " " +
                    f"({R[0]}, {R[1]}, {R[2]}): np.array([" +
                    f"[{matrix[0][0]}, {matrix[0][1]}, {matrix[0][2]}], " +
                    f"[{matrix[1][0]}, {matrix[1][1]}, {matrix[1][2]}], " +
                    f"[{matrix[2][0]}, {matrix[2][1]}, {matrix[2][2]}]]),\n")

                output_line += (
                    f"  {atom1:3} {atom2:3} ({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n" +
                    f"    Isotropic: {J_iso:.4f}\n" +
                    f"    Anisotropic:\n" +
                    f"        {J_aniso[0][0]:7.4f}  " +
                    f"{J_aniso[0][1]:7.4f}  " +
                    f"{J_aniso[0][2]:7.4f}\n" +
                    f"        {J_aniso[1][0]:7.4f}  " +
                    f"{J_aniso[1][1]:7.4f}  " +
                    f"{J_aniso[1][2]:7.4f}\n" +
                    f"        {J_aniso[2][0]:7.4f}  " +
                    f"{J_aniso[2][1]:7.4f}  " +
                    f"{J_aniso[2][2]:7.4f}\n" +
                    f"    DMI: {DMI[0]:.4f} " +
                    f"{DMI[1]:.4f} " +
                    f"{DMI[2]:.4f}\n"
                    f"    |DMI|: {abs_DMI:.4f}\n" +
                    f"    |DMI/J| {abs_DMI/J_iso:.4f}\n" +
                    f"    Matrix:\n" +
                    f"        {matrix[0][0]:7.4f}  " +
                    f"{matrix[0][1]:7.4f}  " +
                    f"{matrix[0][2]:7.4f}\n" +
                    f"        {matrix[1][0]:7.4f}  " +
                    f"{matrix[1][1]:7.4f}  " +
                    f"{matrix[1][2]:7.4f}\n" +
                    f"        {matrix[2][0]:7.4f}  " +
                    f"{matrix[2][1]:7.4f}  " +
                    f"{matrix[2][2]:7.4f}\n\n")

            output_python_iso += "    },\n"
            output_python_aniso += "    },\n"
            output_python_dmi += "    },\n"
            output_python_matrix += "    },\n"

        else:
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
                f"    Isotropic: {J_iso:.4f}\n" +
                f"    Anisotropic:\n" +
                f"        {J_aniso[0][0]:7.4f}  " +
                f"{J_aniso[0][1]:7.4f}  " +
                f"{J_aniso[0][2]:7.4f}\n" +
                f"        {J_aniso[1][0]:7.4f}  " +
                f"{J_aniso[1][1]:7.4f}  " +
                f"{J_aniso[1][2]:7.4f}\n" +
                f"        {J_aniso[2][0]:7.4f}  " +
                f"{J_aniso[2][1]:7.4f}  " +
                f"{J_aniso[2][2]:7.4f}\n" +
                f"    |DMI|: {abs_DMI:.4f}\n" +
                f"    |DMI/J| {abs(abs_DMI/J_iso):.4f}\n")

            if dmi:
                for bond in template.names[name]:
                    atom1 = bond[0]
                    atom2 = bond[1]
                    R = bond[2]
                    DMI = model.bonds[atom1][atom2][R].dmi
                    output_line += (
                        f"    DMI: " +
                        f"{DMI[0]:7.4f} " +
                        f"{DMI[1]:7.4f} " +
                        f"{DMI[2]:7.4f} ({R[0]:2.0f}, {R[1]:2.0f}, {R[2]:2.0f})\n")
                output_line += "\n"
            else:
                output_line += (
                    f"    DMI: " +
                    f"{DMI[0]:.4f} " +
                    f"{DMI[1]:.4f} " +
                    f"{DMI[2]:.4f}\n")
        output_line += "\n"

    if out_name is not None:
        with open(join(out_dir, out_name + ".txt"), "w") as out_file:
            out_file.write(output_line)
        print(f"{OK}Extracted exchange info is in " +
              f"{abspath(join(out_dir, out_name + '.txt'))}{RESET}")
        if verbose:
            output_python_iso += "}\n\n"
            output_python_aniso += "}\n\n"
            output_python_dmi += "}\n\n"
            output_python_matrix += "}\n\n"
            with open(join(out_dir, out_name + ".py"), "w") as out_python:
                out_python.write("import numpy as np\n")
                out_python.write(output_python_iso)
                out_python.write(output_python_aniso)
                out_python.write(output_python_dmi)
                out_python.write(output_python_matrix)
            print(f"{OK}Verbose info is in " +
                  f"{abspath(join(out_dir, out_name + '.py'))}{RESET}")

    else:
        print(f"{output_line}")


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

    main(filename=args.filename,
         out_dir=args.output_dir,
         out_name=args.output_name,
         template=args.template,
         dmi=args.dmi,
         verbose=args.verbose)
