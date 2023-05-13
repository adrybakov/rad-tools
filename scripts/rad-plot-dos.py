#! /usr/local/bin/python3

import sys
from time import sleep

from rad_tools.routines import winwait
from rad_tools.score.plot_dos import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
    winwait()
