import git
import sys
from termcolor import colored, cprint
import os
from calendar import month_name
from datetime import datetime
import re
from argparse import ArgumentParser


class FATAL(Exception):
    pass


ROOT_DIR = os.path.abspath(".")  # Root directory of the project


def envelope(message: str):
    """
    Decorator for printing a message before and "Done" after a function.

    Parameters
    ----------
    message : str
        Message to print before the function.
    """

    def wrapper(func):
        def inner(*args, **kwargs):
            print(f"{f'{message}'} ... ")
            func(*args, **kwargs)
            print(colored("Done", "green"))
            print(f"{'':=^40}")

        return inner

    return wrapper


@envelope(message="Checking git branch")
def check_active_branch(repo: git.Repo):
    """
    Check if the active branch is stable.

    Parameters
    ----------
    repo : git.Repo
        Git repository object.
    """

    if repo.active_branch.name != "stable":
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored("\nYou are not on stable branch\n", "red"),
                    f"You are on '{repo.active_branch.name}' branch.\n",
                    "Please checkout to the stable branch by running\n\n",
                    "    git checkout stable\n",
                ]
            )
        )


@envelope(message="Updating __init__.py")
def update_init(repo: git.Repo, version):
    """
    Update __init__.py file.

    Change the __git_hash__, __release_date__ and __version__ variables.

    Parameters
    ----------
    repo : git.Repo
        Git repository object.
    version : str
        Target version for the release.
    """
    cd = datetime.now()
    sha = repo.head.object.hexsha

    variables = ["__git_hash__", "__release_date__", "__version__"]
    values = [sha, f"{cd.day} {month_name[cd.month]} {cd.year}", version]
    good = [False, False, False]
    good_message = {False: "not updated", True: "updated"}

    # Read __init__.py
    init_file_path = os.path.join(ROOT_DIR, "radtools", "__init__.py")
    init_file = open(init_file_path, "r")
    init_file_content = init_file.readlines()
    init_file.close()

    # Update __init__.py
    for l_i, line in enumerate(init_file_content):
        # Update git hash
        for v_i, variable in enumerate(variables):
            if re.fullmatch(f"{variable}.*\n", line):
                init_file_content[l_i] = f'{variable} = "{values[v_i]}"\n'
                good[v_i] = True

    # Check if all variables were updated
    if not all(good):
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored(
                        "\nFailed to update some variables in '__init__.py':\n", "red"
                    ),
                    *[
                        f"    {variables[i]:20} : {good_message[good[i]]}\n"
                        for i in range(len(variables))
                    ],
                    "Please check the file and try again\n",
                ]
            )
        )

    # Write __init__.py
    init_file = open(init_file_path, "w")
    init_file.writelines(init_file_content)
    init_file.close()


@envelope(message="Checking release notes")
def check_release_notes(version: str):
    """
    Check if the release notes are up to date.

    Parameters
    ----------
    version : str
        Target version for the release.
    """
    path = os.path.join(ROOT_DIR, "docs", "source", "release-notes")
    # (major, minor, rest)
    major, minor, rest = tuple(map(int, version.split(".")[:3]))

    for dirpath, dirnames, filenames in os.walk(path):
        break
    # (major, minor)
    files = []
    for filename in sorted(filenames):
        if re.fullmatch("[0-9]*\.[0-9]*\.rst", filename):
            files.append(tuple(map(int, filename.split(".")[:2])))

    # Check the minor version file
    if (major, minor) not in files:
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored(
                        f"\nRelease notes for the minor version {major}.{minor} are not available\n",
                        "red",
                    ),
                    "Please add the file\n\n",
                    f"    {major}.{minor}.rst\n\n",
                    "to the directory\n\n",
                    f"    docs/source/release-notes/\n",
                ]
            )
        )

    # Update the toctree in the index file
    # Read and process current content
    index_file = open(os.path.join(path, "index.rst"), "r")
    lines = []
    skip_empty = False
    for line in index_file:
        untouched_line = line
        line = line.translate(str.maketrans("", "", " \n"))
        if re.fullmatch(f":maxdepth:1", line):
            lines.append(untouched_line + "\n")
            skip_empty = True
            lines.extend([f"    {major}.{i}\n" for i in range(minor, -1, -1)])
        elif not re.match(f"{major}\.", line) and (not skip_empty or line != ""):
            lines.append(untouched_line)
    index_file.close()
    # Write new content
    index_file = open(os.path.join(path, "index.rst"), "w")
    index_file.writelines(lines)
    index_file.close()

    # Check the version note in the minor version file
    # For the x.x.0 versions there are no note required, but 'Whats new?' section
    file = open(os.path.join(path, f"{major}.{minor}.rst"), "r")
    found_note = False
    for line in file:
        line = line.translate(str.maketrans("", "", " \n"))
        if re.fullmatch(f"v{major}.{minor}.{rest}", line) or (
            re.fullmatch("Whatsnew\?", line) and rest == 0
        ):
            found_note = True
            break
    file.close()
    if not found_note:
        sys.tracebacklimit = 0
        if rest == 0:
            raise FATAL(
                "".join(
                    [
                        colored(
                            f"\n'Whats new?' section for the {major}.{minor}.{rest} version is not available\n",
                            "red",
                        ),
                        f"Please add it to the file\n\n",
                        f"    docs/source/release-notes/{major}.{minor}.rst\n\n",
                        "Follow the style:\n\n",
                        f"    Whats new?\n",
                        f"    {'':-^10}\n",
                        "    Text\n",
                    ]
                )
            )
        else:
            raise FATAL(
                "".join(
                    [
                        colored(
                            f"\nRelease notes for the {major}.{minor}.{rest} version are not available\n",
                            "red",
                        ),
                        f"Please add the note for the {major}.{minor}.{rest} version ",
                        "to the file\n\n",
                        f"    docs/source/release-notes/{major}.{minor}.rst\n\n",
                        "Follow the style:\n\n",
                        f"    {major}.{minor}.{rest}\n",
                        f"    {'':-^{len(version)}}\n",
                        "    Note text\n",
                    ]
                )
            )


@envelope(message="Checking if everything is committed and pushed")
def check_git_status(repo: git.Repo):
    """
    Check if everything is committed and pushed.

    Parameters
    ----------
    repo : git.Repo
        Git repository object.
    """
    status = repo.git.status()
    if "nothing to commit, working tree clean" not in status:
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored("\nThere are uncommitted changes\n", "red"),
                    "Please commit them:\n\n",
                    status,
                ]
            )
        )
    if "Your branch is up to date with" not in status:
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored("\nThere are unpushed changes\n", "red"),
                    "Please push them:\n\n",
                    status,
                ]
            )
        )


def main(version: str):
    if version == "undefined":
        sys.tracebacklimit = 0
        raise FATAL(
            "".join(
                [
                    colored("\nVersion is undefined\n", "red"),
                    "For the make command use the syntax:\n\n",
                    "    make prepare-release VERSION=x.x.x\n",
                ]
            )
        )

    print(f"{'':=^40}\n{f'Preparing {version} release':^40}\n{'':=^40}")
    repo = git.Repo(search_parent_directories=True)

    # the order of checks is important, for example,
    # if the update_init() is called before check_active_branch()
    # it will never pass the check_active_branch() because the
    # update_init() will change the __init__.py file and there
    # will be uncommitted changes

    # check_active_branch(repo)

    check_release_notes(version=version)

    # check_git_status(repo)

    update_init(repo, version=version)

    print(colored(f"{f'{version} ready to deploy':^40}", "green"))
    print(f"{'':=^40}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        metavar="x.x.x",
        type=str,
        required=True,
        help="Version to release",
    )
    args = parser.parse_args()
    main(**vars(args))
