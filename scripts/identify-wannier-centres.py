#! /usr/local/bin/python3

import sys
from time import sleep
from rad_tools.routines import cprint
from rad_tools.score.identify_wannier_centres_core import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
    if sys.platform == "win32":
        cprint("Press Enter to continue", colour="green")
        input()
