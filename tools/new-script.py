from argparse import ArgumentParser
import os
from termcolor import colored
import sys

ROOT_DIR = "."


def main(name: str):
    if name == "None":
        sys.tracebacklimit = 0
        raise ValueError(
            "".join(
                [
                    colored("\nName is undefined\n", "red"),
                    "For the make command use the syntax:\n\n",
                    "    make new-script NAME=name\n",
                ]
            )
        )
    name_py = name.replace("-", "_")
    # Check if the script already exists
    if os.path.isfile(os.path.join(ROOT_DIR, "radtools", "score", f"{name_py}.py")):
        sys.tracebacklimit = 0
        raise FileExistsError(
            "".join(
                [
                    colored("\nScript implementation file already exists\n", "red"),
                    "If you want to rewrite it, delete it manually\n",
                ]
            )
        )
    if os.path.isfile(os.path.join(ROOT_DIR, "scripts", f"rad-{name}.py")):
        sys.tracebacklimit = 0
        raise FileExistsError(
            "".join(
                [
                    colored("\nScript file already exists\n", "red"),
                    "If you want to rewrite it, delete it manually\n",
                ]
            )
        )
    if os.path.isfile(
        os.path.join(
            ROOT_DIR, "docs", "source", "user-guide", "scripts", f"rad-{name}.rst"
        )
    ):
        sys.tracebacklimit = 0
        raise FileExistsError(
            "".join(
                [
                    colored("\nDocumentation file already exists\n", "red"),
                    "If you want to rewrite it, delete it manually\n",
                ]
            )
        )

    # Script interface
    with open(
        os.path.join(ROOT_DIR, "scripts", f"rad-{name}.py"), "w", encoding="utf-8"
    ) as file:
        file.write(
            "#! /usr/local/bin/python3\n"
            + "\n"
            + "from radtools._osfix import _winwait\n"
            + f"from radtools.score.{name_py} import create_parser, manager\n"
            + "\n"
            + 'if __name__ == "__main__":\n'
            + "    parser = create_parser()\n"
            + "    args = parser.parse_args()\n"
            + "    manager(**vars(args))\n"
            + "    _winwait()\n"
        )

    # Script implementation
    with open(
        os.path.join(ROOT_DIR, "radtools", "score", f"{name_py}.py"),
        "w",
        encoding="utf-8",
    ) as file:
        file.write(
            "from argparse import ArgumentParser\n"
            + "\n"
            + "# This function is called by the script\n"
            + "def manager(input_filename):\n"
            + "    pass\n"
            + "\n"
        )

    # Docs
    with open(
        os.path.join(
            ROOT_DIR, "docs", "source", "user-guide", "scripts", f"rad-{name}.rst"
        ),
        "w",
        encoding="utf-8",
    ) as file:
        N = len(f"rad-{name}.py")
        file.write(
            f".. _rad-{name}:\n"
            + "\n"
            + "*" * N
            + "\n"
            + f"rad-{name}.py"
            + "\n"
            + "*" * N
            + "\n"
        )
    # setup.py
    setup_content = []
    with open(os.path.join(ROOT_DIR, "setup.py"), "r") as file:
        for line in file:
            setup_content.append(line)
            if line.startswith("    scripts=["):
                setup_content.append(f'        "scripts/rad-{name}.py",\n')
    with open(os.path.join(ROOT_DIR, "setup.py"), "w", encoding="utf-8") as file:
        file.writelines(setup_content)

    print("Templates are created")
    print(
        f"Write manager() function in \n{os.path.abspath(os.path.join(ROOT_DIR, 'radtools', 'score', f'{name_py}.py'))}"
    )
    print(
        f"Write documentation in \n{os.path.abspath(os.path.join(ROOT_DIR, 'docs', 'source', 'user-guide', 'scripts', f'rad-{name}.rst'))}"
    )
    print(
        "Follow the documentation guidelines for the manager() function. "
        + "When it is ready run\n\n    make generate-script-docs SCRIPT=new-script\n\n"
        + "This command generates Arguments section in docs and "
        + "parser in script implementation"
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=True,
        help='Name of the script with no "rad-" prefix',
    )
    args = parser.parse_args()
    main(args.name)
