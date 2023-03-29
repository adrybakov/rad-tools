r"""
Module for constructing the k-points paths and high-symmetry k-points [1]_. 

Examples
--------

.. code-block:: python

    from rad_tools.kpoints import HighSymmetryPoints

    hexagonal = HighSymmetryPoints().hex()

References
----------
.. [1] 
    Setyawan, W. and Curtarolo, S., 2010. 
    High-throughput electronic band structure calculations: 
    Challenges and tools. 
    Computational materials science, 49(2), pp.299-312.

    :DOI:`10.1016/j.commatsci.2010.05.010`
"""

from rad_tools.kpoints.high_symmetry_point import HighSymmetryPoints
from rad_tools.kpoints.kpoints import KPoints
