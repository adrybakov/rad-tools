from radtools import __version__, __git_hash__, __doclink__, __release_date__

__all__ = ["logo"]


def logo(info=None, line_length=71, flat=False):
    """
    Logo generator for rad-tools package.

    Returns the logo and information about the package.

    Parameters
    ----------
    info : list of str, optional
        Information about the package.
        will be displayed below the logo.
        Each element should not exceed 59 characters.
        by default it displays the version, release date,
        git hash and documentation link. You can pass th empty list
        to display the logo only.
    line_length : int
        Length of the lines to be returned.
        Minimum value is 71.
    flat : bool
        Whether to return a flat logo or not.

    Returns
    -------
    logo_info : str
        Logo and information about the package.
    """
    if info is None:
        info = [
            f"Version: {__version__}",
            f"Release date: {__release_date__}",
            f"Git hash: {__git_hash__}",
            f"Documentation: {__doclink__}",
        ]
    logo = [
        "██████╗  █████╗ ██████╗       ████████╗ █████╗  █████╗ ██╗      ██████╗",
        "██╔══██╗██╔══██╗██╔══██╗      ╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝",
        "██████╔╝███████║██║  ██║█████╗   ██║   ██║  ██║██║  ██║██║     ╚═█████╗",
        "██╔══██╗██╔══██║██║  ██║╚════╝   ██║   ██║  ██║██║  ██║██║       ╚══██║",
        "██║  ██║██║  ██║██████╔╝         ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝",
        "╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝          ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝ ",
    ]
    if flat:
        logo = [
            "██████   █████  ██████        ████████  █████   █████  ██       ██████ ",
            "██   ██ ██   ██ ██   ██          ██    ██   ██ ██   ██ ██      ██      ",
            "██████  ███████ ██   ██ █████    ██    ██   ██ ██   ██ ██        █████ ",
            "██   ██ ██   ██ ██   ██          ██    ██   ██ ██   ██ ██           ██ ",
            "██   ██ ██   ██ ██████           ██     █████   █████  ███████ ██████  ",
        ]
    cat = [
        "▄   ▄     ",
        "█▀█▀█     ",
        "█▄█▄█     ",
        " ███   ▄▄ ",
        " ████ █  █",
        " ████    █",
        " ▀▀▀▀▀▀▀▀ ",
    ]

    N = 71
    n = 12
    if line_length < N:
        line_length = N

    logo_info = [f"{x:^{N}}" for x in logo]
    if len(info) > 0:
        if len(info) <= len(cat):
            before = (len(cat) - len(info)) // 2
            after = len(cat) - len(info) - before
            for i in range(len(cat)):
                if i < before or i >= len(cat) - after:
                    logo_info.append(f"{' ':{N-n}}{cat[i]:^{n}}")
                else:
                    logo_info.append(f"{info[i-before]:^{N-n}}{cat[i]:^{n}}")
        else:
            before = (len(info) - len(cat)) // 2
            after = len(info) - len(cat) - before
            for i in range(len(info)):
                if i < before or i >= len(info) - after:
                    logo_info.append(f"{info[i-before]:^{N-n}}")
                else:
                    logo_info.append(f"{info[i-before]:^{N-n}}{cat[i-before]:^{n}}")

    logo_info = [f"{x:^{line_length}}\n" for x in logo_info]
    return "".join(logo_info)[:-1]
