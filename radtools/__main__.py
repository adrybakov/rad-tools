from radtools import __version__, __git_commit__, __doclink__

if __name__ == "__main__":
    print(
        f"rad-tools package, version {__version__}\n"
        + f"Git hash: {__git_commit__}\n"
        + f"Documentation: {__doclink__}"
    )
