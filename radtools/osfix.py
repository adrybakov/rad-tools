import sys


def _winwait():
    r"""
    Add "Press Enter to continue" behaviour to Windows.

    Its a hotfix for Window`s terminal behaviour.
    """
    if sys.platform == "win32":
        cprint("Press Enter to continue", "green")
        input()
