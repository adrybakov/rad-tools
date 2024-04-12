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
DOCS_DIR = os.path.join(ROOT_DIR, "docs", "source", "user-guide", "scripts")
EXAMPLE_DIR = os.path.join(ROOT_DIR, "docs", "examples")
ARGUMENTS_TO_EXTEND = [
    "-if",
    "--input-filename",
    "--input-folder",
    "-on",
    "--output-name",
    "-tf",
    "--template-file",
]


def main(script="all", debug=False):
    # Get list of files and directories with docs
    for _, tmp_dirnames, tmp_filenames in os.walk(DOCS_DIR):
        break

    dirnames = []
    filenames = []

    # Filter out non-relevant files and directories
    for name in tmp_dirnames:
        if re.fullmatch("rad\-.*", name):
            dirnames.append(re.split("rad\-", name)[1])

    for name in tmp_filenames:
        if re.fullmatch("rad\-.*\.rst", name):
            filenames.append(re.split("\.rst", re.split("rad\-", name)[1])[0])

    # Check versus input
    if re.fullmatch("rad\-.*", script):
        script = re.split("rad\-", script)[1]
    if script != "all":
        if script in dirnames:
            dirnames = [script]
            filenames = []
        elif script in filenames:
            filenames = [script]
            dirnames = []
        else:
            sys.tracebacklimit = 0
            raise ValueError(f"Script {script} not found")

    # Dict: {script-name : [file, ...]}
    files = dict(zip(filenames, [[f"rad-{name}.rst"] for name in filenames]))

    # Parse folder for the actual files
    dirnames = dict([(i, []) for i in dirnames])
    for name in dirnames:
        if name not in files:
            files[name] = []
        for _, _, tmp_filenames in os.walk(os.path.join(DOCS_DIR, f"rad-{name}")):
            for tmp_filename in tmp_filenames:
                if re.fullmatch(".*\.rst", tmp_filename):
                    files[name].append(os.path.join(f"rad-{name}", tmp_filename))
    if debug:
        for name in files:
            print(f"{name}:")
            for filename in files[name]:
                print(f"  {filename}")

    # Parse code from doc files
    parsed_code = {}
    for script in files:
        # Custom behaviour for rad-plot-dos
        if script == "plot-dos":
            os.system(
                f"rm -r {os.path.join(EXAMPLE_DIR, f'rad-plot-dos', 'style-examples')}"
            )
        for filename in files[script]:
            with open(os.path.join(DOCS_DIR, filename), "r") as file:
                for line in file:
                    if re.match(" *\.\. code\-block\:\: bash", line):
                        N = len(line) - len(line.lstrip(" "))
                        file.readline()
                        code = file.readline()
                        while re.match(" " * (N + 4), code):
                            if script not in parsed_code:
                                parsed_code[script] = []
                            parsed_code[script].append(code)
                            code = file.readline()

    # Update arguments
    code_to_run = []
    for script in parsed_code:
        for code in parsed_code[script]:
            code = code.split()
            for e_i, entry in enumerate(code):
                if entry in ARGUMENTS_TO_EXTEND:
                    code[e_i + 1] = os.path.join(
                        EXAMPLE_DIR, f"rad-{script}", code[e_i + 1]
                    )
            if len(code) == 1:
                code.append("-on")
                code.append(os.path.join(EXAMPLE_DIR, f"rad-{script}", ""))
            code_to_run.append(" ".join(code))

    if debug:
        for c in code_to_run:
            print(c)
        print(len(code_to_run))
    else:
        # Run the code
        for code in code_to_run:
            os.system(code)


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
