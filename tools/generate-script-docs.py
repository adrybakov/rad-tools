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
import re
import sys
from argparse import ArgumentParser

ROOT_DIR = "."
DOCS_DIR = os.path.abspath(
    os.path.join(ROOT_DIR, "docs", "source", "user-guide", "scripts")
)
SCRIPT_DIR = os.path.abspath(os.path.join(ROOT_DIR, "radtools", "score"))


def parse_manager(scriptfile, debug=False):
    if debug:
        print(f"Parsing {scriptfile}")
    scriptfile = open(scriptfile, "r")
    line = scriptfile.readline()

    parameters = {}
    required = {}
    default = {}
    # Find the manager function
    while not line.startswith("def manager"):
        line = scriptfile.readline()

    # Read the parameters of the function
    tmp_parameters = []
    if "):" in line:
        line = line.split("(")[1].split("):")[0]
        tmp_parameters = [i.strip() for i in line.split(",")]
    else:
        line = scriptfile.readline()
        while not line.startswith("):"):
            tmp_parameters.append(line.split(",")[0])
            line = scriptfile.readline()

    if debug:
        print("Parameters in the manager parenthesis:")
        for parameter in tmp_parameters:
            print(f"{parameter}")

    for parameter in tmp_parameters:
        parameter = parameter.replace(" ", "").replace("\n", "")
        parameter = parameter.split("=")
        parameters[parameter[0]] = []
        if len(parameter) == 1:
            required[parameter[0]] = True
        else:
            required[parameter[0]] = False
            default[parameter[0]] = parameter[1]

    if debug:
        print("Data extracted from the manager parenthesis:")
        print(f'    {"Name":20} : {"Required":8} : {"Default":10}')
        for parameter in parameters:
            print(
                f"    {parameter:20} : "
                + f"{str(required[parameter]):8} : "
                + f"{default[parameter] if parameter in default else 'nothing':10}"
            )

    # Find the parameters section
    next_line = line
    while not line.startswith("    Parameters") or not next_line.startswith(
        " " * 4 + "-" * 10
    ):
        line = next_line
        next_line = scriptfile.readline()
    # Skip underline
    line = next_line

    read_parameters = True
    current_parameter = None
    types = {}
    while read_parameters:
        line = scriptfile.readline()
        if line.startswith(" " * 4) and line[4] != " ":
            if ":" in line and line.split(":")[0].replace(" ", "") in parameters:
                current_parameter = line.split(":")[0].replace(" ", "")
                types[current_parameter] = line.split(":")[1].split(",")[0].strip()
            else:
                read_parameters = False
        elif line.startswith(" " * 8) or line == "\n":
            parameters[current_parameter].append(line)
        else:
            read_parameters = False

    if debug:
        print("Types from the parameters section:")
        for parameter in types:
            print(f"    {parameter:20} : {types[parameter]}")
        print("Raw data extracted from the parameters section:")
        for parameter in parameters:
            print(f"{parameter}")
            print("".join(parameters[parameter]), end="")

    console_arguments = {}
    metavars = {}
    choices = {}
    version_notes = {}
    docs = {}
    for parameter in parameters:
        version_notes[parameter] = []
        console_arguments[parameter] = None
        metavars[parameter] = None
        choices[parameter] = None
        docs[parameter] = []
        for line in parameters[parameter]:
            if "Console argument" in line:
                line = line.translate(str.maketrans("", "", "` \n")).split(":")[1]
                console_arguments[parameter] = line.split("/")
            elif "Metavar:" in line:
                line = line.split(":")[1].strip()
                metavars[parameter] = line
            elif "Choices:" in line:
                line = line.split(":")[1].strip().split(",")
                choices[parameter] = line
            elif re.match(" *\.\.", line):
                version_notes[parameter].append(line.replace("\n", "").lstrip())
            else:
                docs[parameter].append(line.lstrip(" "))

    if debug:
        print("=" * 80)
        print("Docs:")
        for parameter in docs:
            print(parameter)
            print("".join(docs[parameter]), end="")
        print("Console arguments:")
        for parameter in console_arguments:
            print(f"    {parameter:20} : {console_arguments[parameter]}")
        print("Version notes:")
        for parameter in version_notes:
            print(f"    {parameter:20} : {version_notes[parameter]}")
        print("Required:")
        for parameter in required:
            print(f"    {parameter:20} : {required[parameter]}")
        print("Default:")
        for parameter in default:
            print(f"    {parameter:20} : {default[parameter]}")
        print("Types:")
        for parameter in types:
            print(f"    {parameter:20} : {types[parameter]}")
        print("Metavars:")
        for parameter in metavars:
            print(f"    {parameter:20} : {metavars[parameter]}")
        print("Choices:")
        for parameter in choices:
            print(f"    {parameter:20} : {choices[parameter]}")

    return (
        docs,
        console_arguments,
        version_notes,
        required,
        default,
        types,
        metavars,
        choices,
    )


def main(script="all", debug=False):
    # Get filenames of the script docs
    docs_files = []
    docs_dirs = []
    for _, dirnames, filenames in os.walk(DOCS_DIR):
        docs_files.extend(filenames)
        docs_dirs.extend(dirnames)
        break

    # Prepare data {script-name : doc-file}
    docfiles = {}
    for filename in docs_files:
        if filename.startswith("rad-"):
            docfiles[filename.split(".rst")[0].split("rad-")[1]] = os.path.join(
                DOCS_DIR, filename
            )
    for dirname in docs_dirs:
        if dirname.startswith("rad-"):
            docfiles[dirname.split("rad-")[1]] = os.path.join(
                DOCS_DIR, dirname, "index.rst"
            )

    if debug:
        print("Found documentation files:")
        for scriptname, docfile in docfiles.items():
            print(f"{scriptname} : {docfile}")

    # Check versus input
    if re.fullmatch("rad\-.*", script):
        script = re.split("rad\-", script)[1]
    if re.fullmatch(".*\.py", script):
        script = re.split("\.py", script)[0]

    if debug:
        print(f"Input script: {script}")

    if script != "all":
        if script in docfiles:
            docfiles = {script: docfiles[script]}
        else:
            sys.tracebacklimit = 0
            raise ValueError(f" Documentation file for script {script} not found")

    for script, docfile in docfiles.items():
        # Parse the scripts
        (
            docs,
            console_arguments,
            version_notes,
            required,
            default,
            types,
            metavars,
            choices,
        ) = parse_manager(
            os.path.join(SCRIPT_DIR, f"{script.replace('-','_')}.py"), debug=debug
        )

        # Read the docs, before Arguments
        docs_header = []
        with open(docfile, "r") as file:
            line = file.readline()
            while (
                line
                and not line.startswith("Arguments")
                and not line.startswith(f".. _rad-{script}_arguments:")
            ):
                docs_header.append(line)
                line = file.readline()

        # Write new docs
        with open(docfile, "w", encoding="utf-8") as file:
            file.write("".join(docs_header))

            # Header for the Arguments section
            file.write(f".. _rad-{script}_arguments:\n\nArguments\n=========\n")
            for parameter in docs:
                # Subsection for each parameter
                file.write(f"\n.. _rad-{script}_{parameter.replace('_','-')}:\n\n")
                if len(console_arguments[parameter]) == 1:
                    argheader = f"{console_arguments[parameter][0]}"
                else:
                    console_arguments[parameter].sort(key=lambda x: len(x))
                    argheader = (
                        f"{console_arguments[parameter][0]}, "
                        + f"{console_arguments[parameter][1]}"
                    )
                file.write(f"{argheader}\n")
                file.write("-" * len(argheader) + "\n")
                file.write(
                    "".join(docs[parameter]).rstrip() + "\n\n.. code-block:: text\n\n"
                )
                if required[parameter]:
                    file.write(f"    required\n")
                elif parameter in default and default[parameter] != "None":
                    file.write(f"    default: {default[parameter]}\n")
                else:
                    file.write(f"    optional\n")
                file.write(f"    type: {types[parameter]}\n\n")
                if len(version_notes[parameter]) > 0:
                    for note in version_notes[parameter]:
                        file.write(f"{note}\n")

        # Read the script implementation, before Arguments
        script_header = []
        with open(
            os.path.join(SCRIPT_DIR, f"{script.replace('-','_')}.py"), "r"
        ) as file:
            line = file.readline()
            while line and not line.startswith("def create_parser():"):
                script_header.append(line)
                line = file.readline()

        # Write new script implementation
        with open(
            os.path.join(SCRIPT_DIR, f"{script.replace('-','_')}.py"),
            "w",
            encoding="utf-8",
        ) as file:
            file.write("".join(script_header))
            file.write("def create_parser():\n")
            file.write("\n    parser = ArgumentParser()\n")
            for parameter in docs:
                # Add new parameter
                file.write("    parser.add_argument(\n")
                file.write(f'        "{console_arguments[parameter][0]}",\n')
                if len(console_arguments[parameter]) > 1:
                    file.write(f'        "{console_arguments[parameter][1]}",\n')

                # Write required or default
                if required[parameter]:
                    file.write(f"        required=True,\n")
                else:
                    file.write(f"        default={default[parameter]},\n")

                # Write metavar
                if metavars[parameter] is not None:
                    file.write(f"        metavar={metavars[parameter]},\n")

                # Write type, nargs and action
                if types[parameter] == "bool":
                    file.write(
                        f'        action="store_true",\n'
                        if default[parameter] == "False"
                        else f'        action="store_false",\n'
                    )
                elif "tuple" in types[parameter] or "list" in types[parameter]:
                    tmp_type = types[parameter].split()
                    if len(tmp_type) == 3:
                        file.write(f"        type={tmp_type[2]},\n")
                        file.write(f'        nargs="*",\n')
                    elif len(tmp_type) == 4:
                        file.write(f"        type={tmp_type[3]},\n")
                        file.write(f"        nargs={tmp_type[2]},\n")
                    else:
                        raise ValueError(
                            f"Unknown type {types[parameter]} for parameter {parameter}"
                        )
                else:
                    file.write(f"        type={types[parameter]},\n")

                # Write choices
                if choices[parameter] is not None:
                    file.write("        choices=[\n")
                    for choice in choices[parameter]:
                        file.write(f"            {choice},\n")
                    file.write("        ],\n")

                # Write help
                file.write(f"        help='{docs[parameter][0].strip()}',\n    )\n")
            file.write("\n    return parser\n")
        print(f"Updated {script}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-s", "--script", type=str, default="all", help="Name of the script"
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Print debug information",
    )
    args = parser.parse_args()
    main(**vars(args))
