r"""
Input-output for the files related with this package.
"""

from rad_tools.exchange.template import ExchangeTemplate


def read_template(filename):
    r"""
    Read template from the template file.

    Parameters
    ----------
    filename : str
        Path to the template file.

    Returns
    -------
    template : :py:class:`.ExchangeTemplate`
        Exchange template.
    """

    # Constants
    major_sep = "=" * 20
    minor_sep = "-" * 20
    neighbors_flag = "Neighbors template:"

    template = ExchangeTemplate()

    with open(filename) as file:
        line = True
        while line:
            line = file.readline()

            # Read the neighbors names
            if line and neighbors_flag in line:
                while line and major_sep not in line:
                    if line and minor_sep in line:
                        line = file.readline()
                        name = line.split()[0]
                        try:
                            latex_name = line.split()[1]
                        except IndexError:
                            latex_name = name
                        template.names[name] = []
                        template.latex_names[name] = latex_name
                        line = file.readline()
                        while line and minor_sep not in line and major_sep not in line:
                            atom1 = line.split()[0]
                            atom2 = line.split()[1]

                            R = tuple(map(int, line.split()[2:]))
                            template.names[name].append((atom1, atom2, R))
                            line = file.readline()
                    if line and major_sep in line:
                        break
                    if line and minor_sep not in line:
                        line = file.readline()
    return template
