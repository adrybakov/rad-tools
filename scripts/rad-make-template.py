#! /usr/local/bin/python3

from rad_tools.routines import winwait
from rad_tools.score.make_template import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
    winwait()
