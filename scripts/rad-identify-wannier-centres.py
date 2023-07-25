#! /usr/local/bin/python3

from radtools.utils import winwait
from radtools.score.identify_wannier_centres import create_parser, manager

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    manager(**vars(args))
    winwait()
