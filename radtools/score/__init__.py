r"""
Script interface to the radtools package.

The behaviour of the command line interface is defined by the
functions in this module.  

The functions are called with the same names as the scripts,
but the prefix "rad-" is removed and "-" are substituted by "_". 
Function`s arguments directly correspond to the full names of the
arguments of the script (i.e. the argument :ref:`rad-extract-tb2j_input-filename` 
of the script :ref:`rad-extract-tb2j` is passed to the function :py:func:`.extract_tb2j` as the
argument ``input_filename``). 

Full documentation on the behaviour is available in the 
:ref:`scripts-guide`.


.. admonition:: Example

    Identification of Wannier centres from the file "seedname_centres.xyz"
    with increased span of 0.2, saving the result in the file
    "identified_centres" of the current directory:

    .. code-block:: bash

        rad-identify-wannier-centres.py seedname_centres.xyz -s 0.2 -on identified_centres

    The same result could be achieved by calling the function
    :py:func:`.identify_wannier_centres`:

    .. code-block:: python

        from radtools import identify_wannier_centres  
        identify_wannier_centres("seedname_centres.xyz",
            span = 0.2,
            output_name="identified_centres")
"""

from .extract_tb2j import manager as extract_tb2j
from .identify_wannier_centres import manager as identify_wannier_centres
from .make_template import manager as make_template
from .plot_dos import manager as plot_dos
from .plot_fatbands import manager as plot_fatbands
from .plot_tb2j import manager as plot_tb2j
from .plot_tb2j_magnons import manager as plot_tb2j_magnons

__all__ = [
    "identify_wannier_centres",
    "plot_tb2j",
    "extract_tb2j",
    "plot_dos",
    "make_template",
    "plot_tb2j_magnons",
    "plot_fatbands",
]
