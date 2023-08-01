from os import walk
from os.path import join, split, abspath
from termcolor import cprint, colored

ROOT_DIR = "."
INTERFACE_DIR = abspath(join(ROOT_DIR, "scripts"))
CORE_DIR = abspath(join(ROOT_DIR, "radtools", "score"))
DOCS_DIR = abspath(join(ROOT_DIR, "docs", "source", "user-guide", "scripts"))


def prepare_files():
    """
    Prepare files for comparison

    Parameters
    ----------
    rootdir : str
        Root directory of the project

    Returns
    -------
    scripts : list
        List of script files with the parser creation.
    docs : list
        List of docs files.
    parse_names : list
        List of names of the script base names.
    error_message : str
        Error message if some of docs are missing. Empty string otherwise.
    """

    error_message = ""

    print(f"Interface: {INTERFACE_DIR}")
    print(f"Implementation: {CORE_DIR}")
    print(f"Documentation: {DOCS_DIR}")

    # Get filenames of the script interface
    interface = []
    for dirpath, dirnames, filenames in walk(INTERFACE_DIR):
        interface.extend(filenames)
        break

    # Get filenames of the script implementation
    implementation = []
    for dirpath, dirnames, filenames in walk(CORE_DIR):
        implementation.extend(filenames)
        break

    # Get filenames of the script docs
    docs_files = []
    docs_dirs = []
    for dirpath, dirnames, filenames in walk(DOCS_DIR):
        docs_files.extend(filenames)
        docs_dirs.extend(dirnames)
        break

    scripts = []
    docs = []
    names = []
    for i, filename in enumerate(interface):
        # for each python file
        if ".py" in filename:
            # Check if corresponding docs page exists
            doc_file_name = f"{filename.split('.py')[0]}"
            # As a single file
            if doc_file_name + ".rst" in docs_files:
                docs.append(join(DOCS_DIR, f"{doc_file_name}.rst"))
            # As a directory
            elif doc_file_name in docs_dirs:
                docs.append(join(DOCS_DIR, doc_file_name, "index.rst"))
            # Throw an error if the script is a rad-... script and the docs are missing
            elif filename.split("-")[0] == "rad":
                error_message += colored(
                    f"Documentation for {filename} is missing.\n\n", "red"
                )
                continue
            else:
                continue

            # Check if corresponding script implementation exists (only for "rad-..." styled names)
            score_file_name = f"{filename.split('.py')[0].replace('-','_')[4:]}.py"
            if score_file_name in implementation:
                scripts.append(join(CORE_DIR, score_file_name))
            # Use the script interface file
            else:
                scripts.append(join(INTERFACE_DIR, filename))

            # Save the base name
            names.append(split(filename)[1].split(".")[0])
    return scripts, docs, names, error_message


def check_arguments(names, script_files, doc_files):
    """
    Check if all arguments are present in the docs and in the scripts.

    Parameters
    ----------
    names : list
        List of the script base names.
    script_files : list
        List of script files with the parser creation.
        Same order as the names list.
    doc_files : list
        List of docs files. Same order as the names list.

    Returns
    -------
    ok_message : str
        Ok message with the comparison summary.
    error_message : str
        Error message if some of docs are missing. Empty string otherwise.
    """

    ok_message = []
    error_message = []
    for n_i, name in enumerate(names):
        # Parse scripts
        script_arguments = []
        script_file = open(script_files[n_i], "r")
        for line in script_file:
            if "parser.add_argument" in line:
                arg_name = [script_file.readline().split(",")[0]]
                long_arg_name = script_file.readline()
                if "=" not in long_arg_name:
                    long_arg_name = long_arg_name.split(",")[0]
                    arg_name.append(long_arg_name)
                script_arguments.append(
                    set(
                        [
                            name.translate(str.maketrans("", "", '" '))
                            for name in arg_name
                        ]
                    )
                )
        script_file.close()

        # Parse docs
        docs_arguments = []
        docs_links = []
        docs_file = open(doc_files[n_i], "r")
        lines = docs_file.readlines()
        docs_file.close()
        look_for_arguments = False
        for i, line in enumerate(lines):
            if "Arguments" in line:
                look_for_arguments = True
            if look_for_arguments:
                # .. _script-name_arg-name:
                if f".. _{names[n_i]}_" in line:
                    # Save the rst link:
                    docs_links.append(line.split("_")[2].split(":")[0])

                    # -shortname, --longname
                    # or
                    # -shortname
                    # or
                    # --longname
                    # or
                    # name
                    arg_names = lines[i + 2].split(",")
                    docs_arguments.append(
                        set(
                            [
                                name.translate(str.maketrans("", "", " \n"))
                                for name in arg_names
                            ]
                        )
                    )

        # Compare
        all_args_fine = True
        message = [
            f"Comparing arguments for script {split(script_files[n_i])[1]} \n",
            f"({len(docs_arguments)} in docs, {len(script_arguments)} in scripts)\n",
        ]
        if len(docs_arguments) != len(script_arguments):
            all_args_fine = False
            error_message.append(message[0])
            error_message.append(colored(message[1], "red"))
        else:
            ok_message.append("".join(message))

        for a_i, argument in enumerate(script_arguments):
            try:
                if argument != docs_arguments[a_i] or not (
                    f"--{docs_links[a_i]}" in argument
                    or f"-{docs_links[a_i]}" in argument
                    or f"{docs_links[a_i]}" in argument
                ):
                    all_args_fine = False
                    error_message.append(
                        colored(
                            f"Problem:\n"
                            + f"    Index: {a_i+1}\n"
                            + f"    Docs: {docs_arguments[a_i]}\n"
                            + f"    Docs link: {docs_links[a_i]}\n"
                            + f"    Script: {argument}\n",
                            "red",
                        )
                    )

            except IndexError:
                all_args_fine = False
                error_message.append(
                    colored(
                        f"Problem:\n"
                        + f"    Index: {a_i+1}\n"
                        + f"    Docs: Index error\n"
                        + f"    Script: {argument}\n",
                        "red",
                    )
                )
        if all_args_fine:
            ok_message.append(colored(f"All arguments are fine\n", "green"))
        else:
            error_message.append("\n")

    return "".join(ok_message), "".join(error_message)


if __name__ == "__main__":
    scripts, docs, names, error_message = prepare_files()
    ok_message, second_error_message = check_arguments(names, scripts, docs)
    error_message += second_error_message

    if len(error_message) != 0:
        import sys

        sys.tracebacklimit = 0
        raise Exception("\n" + error_message)
    else:
        cprint("All scripts have corresponding docs.", "green")
