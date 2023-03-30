#! /usr/local/bin/python3

from rad_tools.score.rad_dos_plotter_core import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
