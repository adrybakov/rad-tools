from argparse import ArgumentParser

from rad_tools.tb2j_tools.file_logic import ExchangeModel
from rad_tools.routines import TerminalCoulours


def manager(filename, out_dir,
            max_number=None, max_distance=None, template=None):
    model = ExchangeModel(filename)
    ExchangeModel.filter(distance=max_distance,
                         number=max_number,
                         template=template)


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for visualisation of "
                            "TB2J results")

    parser.add_argument("-f", "--file",
                        type=str,
                        required=True,
                        help="TB2J *.out file")
    parser.add_argument("-op", "--output-dir",
                        type=str,
                        default='.',
                        help="Directory for saving outputs "
                        "if there will be one. "
                        "Could be non-existing, "
                        "but the parent directory have to exist")
    parser.add_argument("-n", "--max-number",
                        type=int,
                        default=None,
                        help="Maximum number of neighbors to be shown. "
                        "Not an order of neighbors, but an exact amount. "
                        "Counting exactly as in TB2J *.out file, "
                        "sorting by distance")
    parser.add_argument("-d", "--max-distance",
                        type=float,
                        default=None,
                        help="Maximum distance for the neighbors to be shown "
                        "(<=)"
                        )
    parser.add_argument("-t", "--template",
                        type=str,
                        default=None,
                        help="Template for filtering the Exchange, it have to be a plain text file "
                        "which will passedf to the Template class")

    args = parser.parse_args()
