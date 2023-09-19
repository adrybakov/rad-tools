from argparse import ArgumentParser

from radtools.decorate.stats import logo, license

if __name__ == "__main__":
    parser = ArgumentParser(
        description="RAD-tools package",
    )
    parser.add_argument(
        "--license",
        action="store_true",
        default=False,
        help="Prints the license of the package.",
    )
    args = parser.parse_args()

    if args.license:
        print(license())
    else:
        print(logo())
