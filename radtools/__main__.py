from argparse import ArgumentParser

from radtools import __version__, __git_hash__, __doclink__, __release_date__

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
            "".join(
                [
                    f"rad-tools package\n",
                    f"version {__version__}\n",
                    f"Release date: {__release_date__}\n",
                    f"Git hash: {__git_hash__}\n",
                    f"Documentation: {__doclink__}",
                ]
            )
        )
