.. _announcement:

*****************
RAD-tools splits!
*****************

**DATE: 11.12.2023**

Part of the RAD-tools package is separated into a standalone package called
`WULFRIC <https://wulfric.org>`_.
The functionality connected with Crystal, Lattice, Atom, Kpoints is moved into WULFRIC.

At the moment the source code of RAD-tools is not changed, but in the future
moved part is going to be removed from it and replaced with the dependency on WULFRIC.

The reason for this change is that the functionality of Crystal, Lattice, Atom, Kpoints is
useful not only for the RAD-tools, but also for other packages.

If you are using RAD-tools only with scripts, you do not need to change anything.

Below is the comparison table, which guides the transition if you are using RAD-tools as
a Python package.

+---------------------------------------+---------------------------------------+
| RAD-tools                             | WULFRIC                               |
+=======================================+=======================================+
| import radtools as rad                | import wulfric as wulf                |
+---------------------------------------+---------------------------------------+
| from radtools import Crystal          | from wulfric import Crystal           |
+---------------------------------------+---------------------------------------+
| from radtools import Lattice          | from wulfric import Lattice           |
+---------------------------------------+---------------------------------------+
| from radtools import Atom             | from wulfric import Atom              |
+---------------------------------------+---------------------------------------+
| from radtools import Kpoints          | from wulfric import Kpoints           |
+---------------------------------------+---------------------------------------+
| from radtools import print_2d_array   | from wulfric import print_2d_array    |
+---------------------------------------+---------------------------------------+
| from radtools import lattice_examples | from wulfric import lattice_examples  |
+---------------------------------------+---------------------------------------+
| from radtools PlotlyBackend           | from wulfric PlotlyBackend            |
+---------------------------------------+---------------------------------------+

This table is not full, but it should give you an idea how to change your code.
