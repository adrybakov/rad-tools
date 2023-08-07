__all__ = ["ColpaFailed", "NotationError"]


class ColpaFailed(Exception):
    r"""
    Raised when Diagonalization via Colpa fails.
    """

    def __init__(self):
        self.message = "Diagonalization via Colpa failed."

    def __str__(self):
        return self.message


class NotationError(ValueError):
    r"""
    Raised when the notation (or individual property) is not defined.

    Gives a summary of the notation (or individual property) and how to set it.

    Parameters
    ----------
    name : str
        Name of the corresponding attribute.

    """

    def __init__(self, name):
        self.message = (
            f"\n\nNotation`s interpretation is not set for the property {name}.\n"
            + f"Set the notation first:\n"
            + f"    SpinHamiltonian.{name} = True  "
            + f"or  SpinHamiltonian.{name} = False\n\n"
            + f"Note: When the attribute is set for the first time it sets the interpretation, "
            + "afterwards it change the notation.\n\n"
            + f"If you want to set the interpretation again, use \n"
            + f"    SpinHamiltonian.set_interpretation({name} = True)"
            + "\nor\n"
            + f"    SpinHamiltonian.set_interpretation({name} = False)\n"
        )

    def __str__(self):
        return self.message
