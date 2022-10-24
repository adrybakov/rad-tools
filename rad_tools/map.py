

from tracemalloc import start


class CalculationMap:
    """
    Mapper.

    Attributes
    ----------
    parameters : dict
        Dictionary of the parameters.

        .. code-block:: python

            {name : values, ...}
    """

    def __init__(self) -> None:
        self.parameters = {}

    def add_parameter(self,
                      name,
                      values=None,
                      range=None):
        """
        Add parameter to the class.

        Parameters
        ----------
        name : str 
            The main identificator of the parameter.
        values : array-like
            exact values of the parameter to be considered.
        range : tuple 
            Three numbers: start, stop, step. 
            If `values` are not specified the range will be used to 
            create them. If (`stop` - `start`) != n * `step`, where n - int
            then the last point will be < `stop`.
        """
        if values is None:
            if range is None:
                raise ValueError("Either values or range need to be specified")
            else:
                values = []
                value = range[0]
                while value <= range[1]:
                    values.append(value)
                    value += range[2]

        self.parameters[name] = values
