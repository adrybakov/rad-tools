r"""
Most of the scripts are moved to the library and could be called through 
the corresponding ``manager`` function, which results in the behaviour 
analogous to the command line interface.

.. hint::
    Long names of the arguments have to be used, i.e. ``input_path``, not ``ip``.

.. code-block:: python

    from rad_tools import identify_wannier_centers

    identify_wannier_centers("seedname_centres.xyz", 
            span = 0.1, 
            output_path = ".", 
            output_name="identified centres", 
            no_colour=False)


"""

from .identify_wannier_centres import (
    manager as identify_wannier_centres,
)
from .plot_tb2j import (
    manager as plot_tb2j,
)
from .extract_tb2j import (
    manager as extract_tb2j,
)
from .plot_dos import (
    manager as plot_dos,
)
from .make_template import (
    manager as make_template,
)

__all__ = [
    "identify_wannier_centres",
    "plot_tb2j",
    "extract_tb2j",
    "plot_dos",
    "make_template",
]
