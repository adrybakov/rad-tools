from argparse import ArgumentParser

from radtools import __version__, __git_hash__, __doclink__, __release_date__
from radtools.decorate.logo import logo


if __name__ == "__main__":
    parser = ArgumentParser(
        description="rad-tools package",
    )
    parser.add_argument(
        "--info",
        action="store_true",
        default=False,
        help="Prints information about the package.",
    )
    args = parser.parse_args()

    if args.info:
        print(
            logo(
                [
                    f"Version: {__version__}",
                    f"Release date: {__release_date__}",
                    f"Git hash: {__git_hash__}",
                    f"Documentation: {__doclink__}",
                ],
            )
        )
