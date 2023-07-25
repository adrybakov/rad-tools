import numpy as np
from termcolor import colored

__all__ = ["print_2d_array"]


def print_2d_array(
    array, fmt=".2f", highlight=False, print_result=True, borders=True, shift=0
):
    r"""
    Print 1D and 2D arrays in the terminal.

    .. versionadded:: 0.7

    .. versionchanged:: 0.7.11 Renamed from ``print_2D_array``

    Parameters
    ----------
    array : (N,) or (N, M) |array_like|_
        Array to be printed.
    fmt : str
        Format string.
    highlight : bool, default False
        Whether to highlight positive and negative values.
        Only works for real-valued arrays.

        .. versionchanged:: 0.7.11 Renamed from ``posneg``

    print_result : bool, default True
        Whether to print the result or return it as a string.
    borders : bool, default True
        Whether to print borders around the array.

        .. versionadded:: 0.7.11

    shift : int, default 0
        Shifts the array to the right by ``shift`` columns.

        .. versionadded:: 0.7.11

    Returns
    -------
    string : str
        String representation of the array.
        Returned only if ``print_result`` is False.
    """

    array = np.array(array)
    if (len(array.shape) == 1 and array.shape[0] != 0) or (
        len(array.shape) == 2 and array.shape[1] != 0
    ):
        # Convert 1D array to 2D array
        if len(array.shape) == 1:
            array = np.array([array])

        # Array dimensions
        N = len(array)
        M = len(array[0])

        # Find the longest number, used for string formatting
        n = max(
            len(f"{np.amax(array.real):{fmt}}"), len(f"{np.amin(array.real):{fmt}}")
        )
        n = max(
            n, len(f"{np.amax(array.imag):{fmt}}"), len(f"{np.amin(array.imag):{fmt}}")
        )

        # Check if input fmt exceeds the longest number
        try:
            n = max(n, int(fmt.split(".")[0]))
        except ValueError:
            pass

        fmt = f"{n}.{fmt.split('.')[1]}"

        def print_number(number, fmt, highlight=False, condition=None):
            if condition is None:
                condition = number
            string = ""
            # Highlight positive and negative values with colours
            if highlight:
                if condition > 0:
                    string += colored(f"{number:{fmt}}", "red", attrs=["bold"])
                elif condition < 0:
                    string += colored(f"{number:{fmt}}", "blue", attrs=["bold"])
                else:
                    string += colored(f"{number:{fmt}}", "green", attrs=["bold"])
            # Print without colours
            else:
                string += f"{number:{fmt}}"
            return string

        def print_complex(number, fmt, highlight):
            string = ""
            if number.real != 0:
                string += print_number(number.real, fmt, highlight)
            else:
                string += " " * len(print_number(number.real, fmt))

            if number.imag > 0:
                sign = "+"
            else:
                sign = "-"

            if number.imag != 0:
                string += f" {sign} i"
                string += print_number(
                    abs(number.imag), f"<{fmt}", highlight, condition=number.imag
                )
            else:
                string += " " * (
                    len(print_number(abs(number.imag), fmt, condition=number.imag)) + 4
                )
            return string

        def print_border(symbol_start, symbol_middle, symbol_end, n):
            string = symbol_start
            for j in range(0, M):
                # If at least one complex value is present in the column
                if np.iscomplex(array[:, j]).any():
                    string += f"{(2*n + 6)*'─'}"
                else:
                    string += f"{(n + 2)*'─'}"
                if j == M - 1:
                    string += symbol_end + "\n"
                else:
                    string += symbol_middle
            return string

        string = ""
        for i in range(0, N):
            substring = " " * shift
            if borders:
                substring += "│"

            for j in range(0, M):
                # Print complex values
                if np.iscomplex(array[:, j]).any():
                    # Print complex part if it is non-zero
                    substring += " " + print_complex(array[i][j], fmt, highlight)
                # Print real values
                else:
                    # Highlight positive and negative values with colours
                    substring += " " + print_number(array[i][j].real, fmt, highlight)
                if borders:
                    substring += " │"
            if i != N - 1 or borders:
                substring += "\n"

            if borders:
                # Header of the table
                if i == 0:
                    symbol_start = "┌"
                    symbol_middle = "┬"
                    symbol_end = "┐"
                    substring = (
                        " " * shift
                        + print_border(symbol_start, symbol_middle, symbol_end, n)
                        + substring
                    )

                # Footer of the table
                if i == N - 1:
                    symbol_start = "└"
                    symbol_middle = "┴"
                    symbol_end = "┘"
                    substring += (
                        " " * shift
                        + print_border(symbol_start, symbol_middle, symbol_end, n)[:-1]
                    )
                # Middle of the table
                else:
                    symbol_start = "├"
                    symbol_middle = "┼"
                    symbol_end = "┤"
                    substring += " " * shift + print_border(
                        symbol_start, symbol_middle, symbol_end, n
                    )

            string += substring
        if print_result:
            print(string)
        else:
            return string
    else:
        if print_result:
            print(None)
        else:
            return None
