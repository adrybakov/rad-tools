from argparse import ArgumentParser

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
        print(logo())
