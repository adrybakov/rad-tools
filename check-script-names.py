from os import walk
from os.path import join, split

GREEN = "\u001b[32m"
YELLOW = "\u001b[33m"
RED = "\u001b[31m"
RESET = "\u001b[0m"

f_scripts = []
for dirpath, dirnames, filenames in walk("scripts"):
    f_scripts.extend(filenames)
    break

f_scripts_core = []
for dirpath, dirnames, filenames in walk(join("radtools", "score")):
    f_scripts_core.extend(filenames)
    break

f_docs = []
for dirpath, dirnames, filenames in walk(join("docs", "source", "user-guide")):
    f_docs.extend(filenames)
    break

f_compare_scripts = []
f_compare_docs = []
parse_names = []
for i, filename in enumerate(f_scripts):
    if ".py" in filename:
        filename = filename.split(".py")[0]
        if f"{filename}.rst" in f_docs:
            f_compare_docs.append(
                join("docs", "source", "user-guide", f"{filename}.rst")
            )
            if f"{filename.replace('-','_')[4:]}.py" in f_scripts_core:
                f_compare_scripts.append(
                    join("radtools", "score", f"{filename.replace('-','_')[4:]}.py")
                )
            else:
                f_compare_scripts.append(join("scripts", f"{filename}.py"))
            parse_names.append(split(filename)[1].split(".")[0])


for f_i, filename in enumerate(f_compare_scripts):
    # Parse scripts
    script_arguments = []
    script_file = open(filename)
    for line in script_file:
        if "parser.add_argument" in line:
            names = script_file.readline().split(",")
            second_name = script_file.readline()
            if "=" not in second_name:
                second_name = second_name.split(",")
                names.extend(second_name)
            script_arguments.append(set())
            for name in names:
                if name != "\n":
                    script_arguments[-1].add(name.replace('"', "").replace(" ", ""))
    script_file.close()

    # Parse docs
    docs_arguments = []
    docs_links = []
    docs_file = open(f_compare_docs[f_i])
    lines = docs_file.readlines()
    docs_file.close()
    look_for_arguments = False
    for i, line in enumerate(lines):
        if "Arguments" in line:
            look_for_arguments = True
        if look_for_arguments:
            if f".. _{parse_names[f_i]}_" in line:
                docs_links.append(line.split("_")[2].split(":")[0])
                names = lines[i + 2].split(",")
                docs_arguments.append(set())
                for name in names:
                    if name != "\n":
                        name = name.split("\n")[0]
                        docs_arguments[-1].add(name.replace('"', "").replace(" ", ""))

    # Compare
    isok = True
    if len(docs_arguments) != len(script_arguments):
        print(RED)
        isok = False
    print(
        f"Comparing arguments for script {split(filename)[1]} \n"
        + f"({len(docs_arguments)} in docs, {len(script_arguments)} in scripts)"
    )
    print(RESET)
    for i, argument in enumerate(script_arguments):
        try:
            if argument != docs_arguments[i] or not (
                f"--{docs_links[i]}" in argument
                or f"-{docs_links[i]}" in argument
                or f"{docs_links[i]}" in argument
            ):
                isok = False
                print(
                    f"{YELLOW}Problem:\n"
                    + f"    Index: {i+1}\n"
                    + f"    Docs: {docs_arguments[i]}\n"
                    + f"    Docs link: {docs_links[i]}\n"
                    + f"    Script: {argument}{RESET}"
                )

        except IndexError:
            isok = False
            print(
                f"{RED}Problem:\n"
                + f"    Index: {i+1}\n"
                + f"    Docs: Index error\n"
                + f"    Script: {argument}{RESET}"
            )
    if isok:
        print(f"{GREEN}All arguments are fine{RESET}")
    print()
