#! /usr/local/bin/python3

from argparse import ArgumentParser
from rad_tools.routines import GREEN, RESET, RED, YELLOW
from rad_tools.jarvis.head import Action


def manager(no_log=False):
    do_log = not no_log
    if do_log:
        pass
    print(f"{GREEN}JARVIS is here\n{RESET}")
    action = Action("Initiate")
    print(f"{GREEN}JARVIS is out of here\n{RESET}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-nl", "--no-log",
                        default=False,
                        action="store_true",
                        help="Use if you dont want to keep the log.")
    args = parser.parse_args()

    manager(no_log=args.no_log)
