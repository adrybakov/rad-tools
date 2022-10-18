class ExchangeTemplate:
    """
    Provide a template for sorting the exchange from TB2J file.

    In addition store the technical details for plotting, orbital decompozition.

    Parameters
    ----------
    filename : str
        Path to the template file.
    """
    # Constants
    _major_sep = '=' * 90 + '\n'
    _minor_sep = '-' * 88 + '\n'
    _neighbors_flag = 'Neighbors template:'

    def __init__(self, filename: str) -> None:

        self.names = {}
        self.latex_names = {}
        with open(filename) as file:
            line = True
            while line:
                line = file.readline()

                # Read the neighbors names
                if line and self._neighbors_flag in line:
                    while line and self._major_sep not in line:
                        if line and self._minor_sep in line:
                            line = file.readline()
                            name = line.split()[0]
                            try:
                                latex_name = line.split()[1]
                            except IndexError:
                                latex_name = name
                            self.names[name] = []
                            self.latex_names[name] = latex_name
                            line = file.readline()
                            while line and\
                                    self._minor_sep not in line and\
                                    self._major_sep not in line:
                                atom1 = line.split()[0]
                                atom2 = line.split()[1]

                                R = tuple(map(int, line.split()[2:]))
                                self.names[name].append((atom1, atom2, R))
                                line = file.readline()
                        if line and self._major_sep in line:
                            break
                        if line and self._minor_sep not in line:
                            line = file.readline()

    def get_list(self):
        """
        Translate bonds present in the template into the list.

        Returns
        -------
        plained_template : list
            List with the bond specifications:
            [(atom_1, atom_2, R), ...]
        """
        plained_template = []
        for name in self.names:
            for atom1, atom2, R in self.names[name]:
                plained_template.append((atom1, atom2, R))
        return plained_template


if __name__ == "__main__":
    tmp = ExchangeTemplate(
        "/Users/rybakov.ad/Projects/rad-tools/utest/tb2j_tools/resourses/template.txt")
    print('here')
