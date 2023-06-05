r"""
Most of the scripts are moved to the library and could be called through 
the corresponding function, which results in the behaviour 
analogous to the command line interface. The common logic for the name: 
If script is called "rad-script-name", then the function is called "script-name".

.. hint::
    Long names of the arguments have to be used, i.e. ``input_path``, not ``ip``.

.. admonition:: Example

    Identifying Wannier centres from the file "seedname_centres.xyz"
    with increased span of 0.2, saving the result in the file 
    "identified_centres" of the current directory:

    .. code-block:: python

        from radtools import identify_wannier_centres
        
        identify_wannier_centres("seedname_centres.xyz", 
            span = 0.2, 
            output_path = ".", 
            output_name="identified_centres",
            no_colour=False)

    Makings template based on the "exchange.out" TB2J file, 
    filtering by maximum distance of 5, saving the output in 
    the file "template.txt" of the current folder:

    .. code-block:: python

        from radtools import make_template
        
        make_template(output_name="template",
            input_filename="exchange.out",
            max_distance=5)
"""

from .extract_tb2j import manager as extract_tb2j
from .identify_wannier_centres import manager as identify_wannier_centres
from .make_template import manager as make_template
from .plot_dos import manager as plot_dos
from .plot_tb2j import manager as plot_tb2j

__all__ = [
    "identify_wannier_centres",
    "plot_tb2j",
    "extract_tb2j",
    "plot_dos",
    "make_template",
]
