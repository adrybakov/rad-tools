class ExchangeTemplate:
    """
    Provide a template for sorting the exchange from TB2J file.

    In addition store the technical details for plotting, orbital decompozition.

    Parameters
    ----------
    filename : str
        Path to the template file.

    Methods
    -------
    
    _prepare_dict
        Prepare the array of nested dicts.

    _read_neighbors
        Read and format the Neighbors section.

    _plain_template
        Translate ``template`` into the list.
    
    """

    def __init__(self, filename: str) -> None:

        self.template = {}
        self.names_for_plot = {}
        self.plained_template = []

        self.__headers_to_functions = {
            'Neighbors template:\n': self._read_neighbors
        }
        self.__major_sep = '=' * 90 + '\n'
        self.__minor_sep = '-' * 88 + '\n'
        with open(filename, 'r') as file:
            self.__file_data = file.readlines()
        i = 0
        while i < len(self.__file_data):
            if self.__file_data[i] in self.__headers_to_functions:
                i = self.__headers_to_functions[self.__file_data[i]](i)
            else:
                i += 1

    def _prepare_dict(self, dict: dict, atom_1: str, atom_2: str):
        """
        Prepare the array of nested dicts.

        Parameters
        ----------
        dict : dict
            dictionary to be prepared.
        atom_1 : str
            mark of the first-level atom. 
        atom_2 : str
            mark of the second-level atom. 
        """

        if atom_1 not in dict:
            dict[atom_1] = {}
        if atom_2 not in dict[atom_1]:
            dict[atom_1][atom_2] = {}

    def _read_neighbors(self, i: int):
        """
        Read and format the Neighbors section.

        ``__data`` attribute and corresponding public attributes will be
        filled with content here.

        Parameters
        ----------
        i : int
            Index of the line with Neighbor section header
            (see ``__headers_to_functions`` ) in the ``__file_data`` .

        Returns
        -------
        i : int
            Index of next line after the last line of Neighbors section in the
            ``__file_data`` .
        """
        while self.__major_sep not in self.__file_data[i]:
            if self.__minor_sep in self.__file_data[i]:
                i += 1
                key = self.__file_data[i].split()[0]
                self.template[key] = {}
                try:
                    self.names_for_plot[key] = self.__file_data[i].split()[1]
                except IndexError:
                    self.names_for_plot[key] = key
                i += 1
                while (self.__minor_sep not in self.__file_data[i] and
                       self.__major_sep not in self.__file_data[i]):
                    atom_1, atom_2 = tuple(self.__file_data[i].split()[:2])
                    R = tuple(map(int, self.__file_data[i].split()[2:5]))
                    self._prepare_dict(self.template[key], atom_1, atom_2)
                    self.template[key][atom_1][atom_2][R] = True
                    i += 1
            else:
                i += 1
        self._plain_template()
        return i

    def _plain_template(self):
        """
        Translate ``template`` into the list.

        ``plained_template`` = [(atom_1, atom_2, R), ...]
        """
        self.plained_template = []
        for name in self.template:
            for atom_1 in self.template[name]:
                for atom_2 in self.template[name][atom_1]:
                    for R in self.template[name][atom_1][atom_2]:
                        self.plained_template.append((atom_1,
                                                      atom_2,
                                                      R))
