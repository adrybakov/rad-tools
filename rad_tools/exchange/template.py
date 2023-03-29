r"""
Exchange template.

Write a tutorial with docstrings here.
"""


class ExchangeTemplate:
    r"""
    Store a template for sorting the exchange from TB2J file.

    In addition stores the technical details for plotting,
    orbital decomposition, etc.

    Attributes
    ----------
    names : dict
        Dictionary of neighbors from the template file.

        .. code-block:: python

            {name : [(atom1, atom2, R), ...], ...}

    latex_names : dict
        The dictionary of Latex version of names from `names`.

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
        template_list : list
            List with the bond specifications:

            .. code-block:: python

                [(atom_1, atom_2, R), ...]
        """
        template_list = []
        for name in self.names:
            for atom1, atom2, R in self.names[name]:
                template_list.append((atom1, atom2, R))
        return template_list
