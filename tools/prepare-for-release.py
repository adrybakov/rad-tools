import git
import sys
from termcolor import colored
from os.path import join, abspath
from calendar import month_name
from datetime import datetime
import re

ROOT_DIR = abspath(".")  # Root directory of the project


def check_active_branch(repo: git.Repo):
    if repo.active_branch.name != "stable":
        sys.tracebacklimit = 0
        raise Exception(
            "".join(
                [
                    colored("\nYou are not on stable branch\n", "red"),
                    f"You are on '{repo.active_branch.name}' branch.\n",
                    "Please checkout to the stable branch by running\n\n",
                    "    git checkout stable\n",
                ]
            )
        )


def update_init(repo: git.Repo):
    print("Updating __init__.py...")

    cd = datetime.now()
    sha = repo.head.object.hexsha

    # Read __init__.py
    init_file_path = join(ROOT_DIR, "radtools", "__init__.py")
    init_file = open(init_file_path, "r")
    init_file_content = init_file.readlines()
    init_file.close()

    # Update __init__.py
    for l_i, line in enumerate(init_file_content):
        # Update git hash
        if re.fullmatch("__git_hash__.*\n", line):
            init_file_content[l_i] = f'__git_hash__ = "{sha}"\n'
        # Update release date
        if re.fullmatch("__release_date__.*\n", line):
            init_file_content[
                l_i
            ] = f'__release_date__ = "{cd.day} {month_name[cd.month]} {cd.year}"\n'

    # Write __init__.py
    init_file = open(init_file_path, "w")
    init_file.writelines(init_file_content)
    init_file.close()

    print("Done")


def main():
    repo = git.Repo(search_parent_directories=True)
    # check_active_branch(repo)
    update_init(repo)


if __name__ == "__main__":
    main()
