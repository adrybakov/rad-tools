import numpy as np
from termcolor import colored

__all__ = ["print_2d_array"]


def print_2d_array(
    array,
    fmt=">.2f",
    highlight=False,
    print_result=True,
    borders=True,
    shift=0,
    header_row=None,
    footer_row=None,
    header_column=None,
    footer_column=None,
):
    r"""
    Print 1D and 2D arrays in the terminal.

    .. versionadded:: 0.7

    .. versionchanged:: 0.7.11 Renamed from ``print_2D_array``

    Parameters
    ----------
    array : (N,) or (N, M) |array_like|_
        Array to be printed.
    fmt : str, default ">.2f"
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

    header_row : list, optional
        Header of the table (top). Has to have the same length as the number of columns.
        Has to had one more element for each if ``header_column`` or ``footer_column`` is not None.

        .. versionadded:: 0.8.0

    footer_row : list, optional
        Footer of the table (bottom). Has to have the same length as the number of columns.
        Has to had one more element for each if ``header_column`` or ``footer_column`` is not None.

        .. versionadded:: 0.8.0

    header_column : list, optional
        Header of the table (left). Has to have the same length as the number of rows.

        .. versionadded:: 0.8.0

    footer_column : list, optional
        Footer of the table (right). Has to have the same length as the number of rows.

        .. versionadded:: 0.8.0

    Returns
    -------
    string : str
        String representation of the array.
        Returned only if ``print_result`` is False.
    """

    top_border = ("┌", "┬", "┐")
    middle_border = ("├", "┼", "┤")
    bottom_border = ("└", "┴", "┘")
    if borders:
        vert = "│"
        space_vert = f" {vert}"
    else:
        vert = ""
        space_vert = ""

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

        # Check if header_row, footer_row, header_column and footer_column have the correct length
        if header_column is not None and len(header_column) != N:
            raise ValueError(
                f"header_column has to have the same length"
                + f" as the number of rows ({N})."
                + f"It has length of {len(header_column)}. "
            )
        if footer_column is not None and len(footer_column) != N:
            raise ValueError(
                f"footer_column has to have the same length"
                + f" as the number of rows ({N})."
                + f"It has length of {len(footer_column)}. "
            )
        if header_row is not None and len(header_row) != M + int(
            header_column is not None
        ) + int(footer_column is not None):
            raise ValueError(
                f"header_row has to have the same length "
                + f"({M + int(header_column is not None) + int(footer_column is not None)})"
                + f" as the number of columns ({M}) + one more element for each if "
                + "header_column or footer_column is not None "
                + f"(+{int(header_column is not None) + int(footer_column is not None)}). "
                + f"It has length of {len(header_row)}."
            )
        if footer_row is not None and len(footer_row) != M + int(
            header_column is not None
        ) + int(footer_column is not None):
            raise ValueError(
                f"footer_row has to have the same length "
                + f"({M + int(header_column is not None) + int(footer_column is not None)})"
                + f" as the number of columns ({M}) + one more element for each if "
                + "header_column or footer_column is not None "
                + f"(+{int(header_column is not None) + int(footer_column is not None)}). "
                + f"It has length of {len(footer_row)}."
            )

        # Fix issue#3
        array[array == 0] = 0

        # Define functions for printing numbers and borders
        def print_number(number, fmt, highlight=False, condition=None):
            if condition is None:
                condition = number
            string = ""
            if number is None:
                return " " * len(f"{1:{fmt}}")
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

        def print_complex(number, fmt, highlight=False):
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
            result = [symbol_start]
            if header_column is not None:
                result.append(f"{(n[0]+2)*'─'}{symbol_middle}")
                n_index_shift = 1
            else:
                n_index_shift = 0
            for j in range(0, M):
                # If at least one complex value is present in the column
                if np.iscomplex(array[:, j]).any():
                    result.append(f"{(2*n[j + n_index_shift] + 6)*'─'}")
                else:
                    result.append(f"{(n[j + n_index_shift] + 2)*'─'}")
                if j != M - 1 or footer_column is not None:
                    result.append(symbol_middle)

            if footer_column is not None:
                result.append(f"{(n[-1]+2)*'─'}")
            result.append(f"{symbol_end}\n")
            return "".join(result)

        # Get maximum string length for each column of an array
        n = np.zeros(M, dtype=int)
        n_full = np.zeros(M, dtype=int)
        for column in range(M):
            for row in range(N):
                if array[row][column] is None or np.isnan(array[row][column]):
                    pass
                elif np.iscomplex(array[row, column]):
                    n[column] = max(
                        n[column], (len(print_complex(array[row][column], fmt)) - 4) / 2
                    )
                    n_full[column] = max(
                        n_full[column], len(print_complex(array[row][column], fmt))
                    )
                else:
                    n[column] = max(
                        n[column], len(print_number(array[row][column].real, fmt))
                    )

        # Parse fmt
        tmp_fmt, post_fmt = fmt.split(".")
        mid_fmt = "".join(c for c in tmp_fmt if c.isdigit())
        pre_fmt = "".join(c for c in tmp_fmt if not c.isdigit())

        # Force provided fmt
        if len(mid_fmt) != 0:
            n = np.amax([n, [int(mid_fmt) for _ in n]], axis=0)

        # Get format for headers and footers
        if header_column is not None:
            tmp_n = max([len(str(x)) for x in header_column])
            n = np.concatenate(([tmp_n], n))
            n_full = np.concatenate(([tmp_n], n_full))
        if footer_column is not None:
            tmp_n = max([len(str(x)) for x in footer_column])
            n = np.concatenate((n, [tmp_n]))
            n_full = np.concatenate((n_full, [tmp_n]))
        if header_row is not None:
            n = np.amax([n, [len(str(x)) for x in header_row]], axis=0)
        if footer_row is not None:
            n = np.amax([n, [len(str(x)) for x in footer_row]], axis=0)

        n_full = np.amax([n, n_full], axis=0)

        string = []
        # Open borders
        if borders:
            string.append(" " * shift + print_border(*top_border, n))
        # Write header row
        if header_row is not None:
            string.append(" " * shift + vert)
            string.extend(
                [
                    f" {x:{pre_fmt}{n_full[x_i]}}{space_vert}"
                    for x_i, x in enumerate(header_row)
                ]
            )
            string.append("\n")
            if borders:
                string.append(" " * shift + print_border(*middle_border, n))

        # Write array
        for i in range(0, N):
            substring = [" " * shift, vert]
            if header_column is not None:
                substring.append(f" {header_column[i]:{pre_fmt}{n[0]}}{space_vert}")
                n_index_shift = 1
            else:
                n_index_shift = 0

            for j in range(0, M):
                fmt = f"{pre_fmt}{n[j + n_index_shift]}.{post_fmt}"
                substring.append(" ")
                # Print blank for None
                if array[i][j] is None or np.isnan(array[i][j]):
                    # Highlight positive and negative values with colours
                    substring.append(
                        " " * (len(f"{'':{pre_fmt}{n[j + n_index_shift]}}"))
                    )
                # Print complex values
                elif np.iscomplex(array[:, j]).any():
                    # Print complex part if it is non-zero
                    substring.append(print_complex(array[i][j], fmt, highlight))
                # Print real values
                else:
                    # Highlight positive and negative values with colours
                    substring.append(print_number(array[i][j].real, fmt, highlight))
                substring.append(space_vert)

            if footer_column is not None:
                substring.append(f" {footer_column[i]:{pre_fmt}{n[-1]}}{space_vert}")

            substring.append("\n")

            if borders and i != N - 1:
                # Middle of the table
                substring.append(" " * shift + print_border(*middle_border, n))

            string.extend(substring)

        # Write footer row
        if footer_row is not None:
            if borders:
                string.append(" " * shift + print_border(*middle_border, n))
            string.append(" " * shift + vert)
            string.extend(
                [
                    f" {x:{pre_fmt}{n_full[x_i]}}{space_vert}"
                    for x_i, x in enumerate(footer_row)
                ]
            )
            string.append("\n")
        # Close borders
        if borders:
            string.append(" " * shift + print_border(*bottom_border, n))

        # Print or return result
        string = "".join(string)
        if string[-1] == "\n":
            string = string[:-1]
        if print_result:
            print(string)
        else:
            return string
    else:
        if print_result:
            print(None)
        else:
            return None


if __name__ == "__main__":
    array = [[1, 2], [3, 4], [52414345345, 6]]
    print_2d_array(array)
    array = [[0, 1.0, -0.0], [1, 1, 1]]
    print_2d_array(array)
    array = [[1, 2], [3, 4], [5, 6]]
    print_2d_array(array)
    print_2d_array(array, fmt="10.2f")
    print_2d_array(array, borders=False)
    print_2d_array(array, shift=3)
    array = [[1, 2], [3, 4], [52414345345, 6]]
    print_2d_array(array, fmt="10.2E")
    array = [[1, 2 + 1j], [3, 4], [52, 6]]
    print_2d_array(array)
    print_2d_array(array, fmt="4.2E")
    array = [[1, 2 - 1j], [3, 4], [52, 6]]
    print_2d_array(array)
    array = [[1, 1j], [3, 4], [52, 6]]
    print_2d_array(array)
    print_2d_array([])
    print_2d_array([[]])
    array = [[1, 2], [3, 4], [52414345345, 6]]
    print_2d_array(array)

    print_2d_array(array, header_row=["a", "b"])
    print_2d_array(array, header_row=["a", "adfadsasfdb"])
    print_2d_array(array, footer_row=["a", "b"])
    print_2d_array(array, header_column=["a", "b", "c"])
    print_2d_array(array, footer_column=["a", "b", "c"])
    print_2d_array(array, footer_column=["a1234234123", "b", "c"])
    print_2d_array(array, header_row=["", "a", "b"], header_column=["a", "b", "c"])
    print_2d_array(array, footer_row=["a", "b", ""], footer_column=["a", "b", "c"])
    print_2d_array(
        array,
        header_row=["", "a", "b", ""],
        header_column=["a", "b", "c"],
        footer_row=["", "a", "b", ""],
        footer_column=["a", "b", "c"],
        shift=3,
    )
    print_2d_array(
        array,
        header_row=["", "a", "b"],
        header_column=["a", "b", "c"],
        footer_row=["", "a", "b"],
    )

    print_2d_array([1, None], header_row=["a", "b"])
    print_2d_array([None, 1])
    print_2d_array([None, 1], header_row=["a", "b"])
    print_2d_array([1, None, 1], header_row=["a", "b", "c"])
    print_2d_array(
        array,
        header_row=["", "a", "b", ""],
        header_column=["a", "b", "c"],
        footer_row=["", "a", "b", ""],
        footer_column=["a", "b", "c"],
        borders=False,
    )
    array = [[1, 2], [3, 4], [5, 6 + 1j]]
    print_2d_array(array, header_row=["a", "b"], fmt="^.2f")
