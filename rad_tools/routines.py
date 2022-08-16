class TerminalCoulours:

    BLACK = '\u001b[30m'
    RED = '\u001b[31m'
    GREEN = '\u001b[32m'
    YELLOW = '\u001b[33m'
    BLUE = '\u001b[34m'
    MAGENTA = '\u001b[35m'
    CYAN = '\u001b[36m'
    WHITE = '\u001b[37m'
    RESET = ' \u001b[0m'

    WARNING = YELLOW
    OK = GREEN
    ERROR = RED


def get_256_colours(n: int):
    if type(n) != int or not 0 <= n <= 255:
        raise ValueError(f'Integer n have to be in range 0 <= n <= 255. '
                         f'You provided n = {n}, type<n> = {type(n)}')
    return f'\033[38:5:{n}m'
