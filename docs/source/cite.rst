.. _rad-tools_cite:

***************
How to cite us?
***************

At the moment the package did not reach the stage of the publication.
However, if you use it, please cite in the following form:

.. code-block::
        
    RAD-tools v<version>, rad-tools.org. Git hash: <hash>.

It is not necessary to include the Git hash but it is recommended for the
reproducibility of the calculations.

Use the following BibTeX entry:

.. code-block:: LaTeX

    @online{radtools,
        author = {Andrey Rybakov},
        title  = {RAD-tools v<version>},
        year   = {2022},
        url    = {rad-tools.org},
        note   = {Git hash: <hash>}
    }

You can consult the version and Git hash in the txt output files or by running
the following command in the terminal (you may need to use ``python3``):

.. code-block:: bash

    python -m radtools
