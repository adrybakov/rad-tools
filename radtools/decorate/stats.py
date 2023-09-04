from calendar import month_name
from datetime import datetime

from radtools import __doclink__, __git_hash__, __release_date__, __version__

__all__ = ["logo", "stamp_line", "license"]


def logo(info=None, line_length=71, flat=False, date_time=False, comment=None):
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
    date_time : bool, default False
        Whether to include the date and time to the standard info or not.
    comment : str or bool, optional
        Whether to use some character at the end of each string. If bool and
        True, then "# " is used. If str, then this string is used. If None, then
        no character is used.
    Returns
    -------
    logo_info : str
        Logo and information about the package.
    """
    if info is None:
        info = [
            f"Version: {__version__}",
            f"Documentation: {__doclink__}",
            f"Release date: {__release_date__}",
            f"Git hash: {__git_hash__}",
            f"Licence: MIT License",
        ]
        if date_time:
            cd = datetime.now()
            info.append("")
            info.append(
                f"Generated on {cd.day} {month_name[cd.month]} {cd.year}"
                + f" at {cd.hour}:{cd.minute}:{cd.second} "
            )
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
    if isinstance(comment, bool) and comment:
        comment = "# "
    elif comment is not None:
        comment = str(comment)
    else:
        comment = ""

    logo_info = [f"{x:^{N}}" for x in logo]
    if len(info) > 0:
        if len(info) <= len(cat):
            before = (len(cat) - len(info)) // 2 + (len(cat) - len(info)) % 2
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

    logo_info = [f"{comment}{x:^{line_length}}\n" for x in logo_info]
    return "".join(logo_info)[:-1]


def stamp_line(date_time=True, version=True, githash=False, doclink=False):
    """
    Return one-line information about the package.

    Parameters
    ----------
    date_time : bool, default True
        Whether to include the release date or not.
    version : bool, default True
        Whether to include the version number or not.
    githash : bool, default False
        Whether to include the git hash or not.
    doclink : bool, default False
        Whether to include the documentation link or not.

    Returns
    -------
    info : str
        Information about the package.
    """

    line = []
    if date_time:
        cd = datetime.now()
        line.append(
            f"on {cd.day} {month_name[cd.month]} {cd.year}"
            + f" at {cd.hour}:{cd.minute}:{cd.second} "
        )
    line.append("by rad-tools ")
    if version:
        line.append(f"{__version__} ")
    if githash:
        line.append(f"(githash {__git_hash__}) ")
    if doclink:
        line.append(f"Documentation: {__doclink__}")
    return "".join(line)


def license():
    """
    Return the license of the package.

    Returns
    -------
    license : str
        License of the package.
    """
    license = [
        "MIT License\n",
        "\n",
        "Copyright (c) 2022-2023 Andrey Rybakov\n",
        "\n",
        "Permission is hereby granted, free of charge, to any person obtaining a copy\n",
        'of this software and associated documentation files (the "Software"), to deal\n',
        "in the Software without restriction, including without limitation the rights\n",
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n",
        "copies of the Software, and to permit persons to whom the Software is\n",
        "furnished to do so, subject to the following conditions:\n",
        "\n",
        "The above copyright notice and this permission notice shall be included in all\n",
        "copies or substantial portions of the Software.\n",
        "\n",
        'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n',
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n",
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n",
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n",
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n",
        "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n",
        "SOFTWARE.\n",
        "",
    ]

    return "".join(license)
