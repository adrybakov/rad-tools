from argparse import ArgumentParser

from rad_tools.tb2j_tools.file_logic import ExchangeModel
from rad_tools.routines import TerminalCoulours


def argument_mistakes_treatment(max_number, max_distance):
    if max_number is not None and max_distance is not None:
        print(f'{TerminalCoulours.WARNING}'
              f'Maximum number of neighbors ({max_number})'
              f' and maximum distance ({max_distance}) are specified '
              f'at the same moment, maximum number will be used.'
              f'{TerminalCoulours.RESET}')


def manager(filename, out_dir):
    model = ExchangeModel(filename)


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for visualisation of "
                            "TB2J results")

    parser.add_argument("-f", '--file',
                        type=str,
                        required=True,
                        help="TB2J *.out file")
    parser.add_argument("-op", '--output-dir',
                        type=str,
                        default='.',
                        help="Directory for saving outputs "
                        "if there will be one. "
                        "Could be non-existing, "
                        "but the parent directory have to exist")
    parser.add_argument("-n", '--max-number',
                        type=int,
                        default=None,
                        help="Maximum number of neighbors to be shown. "
                        "Not an order of neighbors, but an exact amount. "
                        "Counting exactly as in TB2J *.out file, "
                        "sorting by distance")
    parser.add_argument("-d", '--max-distance',
                        type=float,
                        default=None,
                        help="Maximum distance for the neighbors to be shown "
                        "(<=)"
                        )

    args = parser.parse_args()

    argument_mistakes_treatment(args.max_number, args.max_distance)
