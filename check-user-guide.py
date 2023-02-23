from os import walk
from os.path import join

GREEN = "\u001b[32m"
YELLOW = "\u001b[33m"
RED = "\u001b[31m"
RESET = "\u001b[0m"

f_scripts = []
for (dirpath, dirnames, filenames) in walk(join("rad_tools", "score")):
    f_scripts.extend(filenames)
    break

f_docs = []
for (dirpath, dirnames, filenames) in walk(join("docs", "source", "user-guide")):
    f_docs.extend(filenames)
    break

f_compare = []
for i, filename in enumerate(f_scripts):
    if ".py" in filename:
        filename = filename.split("_core.py")[0]
        filename = filename.replace("_", "-")
        if f"{filename}.rst" in f_docs:
            f_compare.append(filename)

for filename in f_compare:

    # Parse scripts
    script_arguments = []
    script_file = open(join("rad_tools", "score",
                       f"{filename.replace('-', '_')}_core.py"))
    for line in script_file:
        if "parser.add_argument" in line:
            names = line.split("(")[1].split(",")
            script_arguments.append(set())
            for name in names:
                if name != "\n":
                    script_arguments[-1].add(name.replace('"',
                                                          "").replace(" ", ""))
    script_file.close()

    # Parse docs
    docs_arguments = []
    docs_file = open(join("docs", "source", "user-guide", f"{filename}.rst"))
    lines = docs_file.readlines()
    docs_file.close()
    look_for_arguments = False
    for i, line in enumerate(lines):
        if "Arguments" in line:
            look_for_arguments = True
        if look_for_arguments:
            if f".. _{filename}_" in line:
                names = lines[i+2].split(",")
                docs_arguments.append(set())
                for name in names:
                    if name != "\n":
                        name = name.split("\n")[0]
                        docs_arguments[-1].add(name.replace('"',
                                                            "").replace(" ", ""))

    # Compare
    print(f"Comparing arguments for script {filename}.py")
    isok = True
    for i, argument in enumerate(script_arguments):
        try:
            if argument != docs_arguments[i]:
                isok = False
                print(f"{YELLOW}Problem:\n" +
                      f"    Index: {i+1}\n" +
                      f"    Docs: {docs_arguments[i]}\n" +
                      f"    Script: {argument}{RESET}")
        except IndexError:
            isok = False
            print(f"{RED}Problem:\n" +
                  f"    Index: {i+1}\n" +
                  f"    Docs: Index error\n" +
                  f"    Script: {argument}{RESET}")
    if isok:
        print(f"{GREEN}All arguments are fine{RESET}")
    print()
