#! /usr/local/bin/python3

from rad_tools.score.rad_plot_dos_core import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
