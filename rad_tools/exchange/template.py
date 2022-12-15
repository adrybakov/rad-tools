r"""
Exchange template.
"""


class ExchangeTemplate:
    r"""
    Store a template for sorting the exchange from TB2J file.

    In addition stores the technical details for plotting, 
    orbital decompozition, etc.

    Attributes
    ----------
    names : dict
        Dictionary of neighbours from the template file.

        .. code-block:: python

            {name : [(atom1, atom2, R), ...], ...}

    latex_names : dict
        The dictionary of Letex version of names from `names`.

        .. code-block:: python

            {name : latex_name, ...}
    """

    def __init__(self) -> None:

        self.names = {}
        self.latex_names = {}

    def get_list(self):
        r"""
        Translate bonds present in the template into the list.

        Returns
        -------
        plained_template : list
            List with the bond specifications:

            .. code-block:: python

                [(atom_1, atom_2, R), ...]
        """
        plained_template = []
        for name in self.names:
            for atom1, atom2, R in self.names[name]:
                plained_template.append((atom1, atom2, R))
        return plained_template
