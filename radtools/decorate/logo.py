__all__ = ["logo"]


def logo(info, line_length=71):
    """
    Logo generator for rad-tools package.

    Returns the logo and information about the package.

    Parameters
    ----------
    info : list of str
        Information about the package.
        will be displayed below the logo.
        Each element should not exceed 59 characters.
    line_length : int
        Length of the lines to be returned.
        Minimum value is 71.

    Returns
    -------
    logo_info : str
        Logo and information about the package.
    """
    logo = [
        "██████╗  █████╗ ██████╗       ████████╗ █████╗  █████╗ ██╗      ██████╗",
        "██╔══██╗██╔══██╗██╔══██╗      ╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝",
        "██████╔╝███████║██║  ██║█████╗   ██║   ██║  ██║██║  ██║██║      ╚█████╗",
        "██╔══██╗██╔══██║██║  ██║╚════╝   ██║   ██║  ██║██║  ██║██║      ╚═══██╗",
        "██║  ██║██║  ██║██████╔╝         ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝",
        "╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝          ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝ ",
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
